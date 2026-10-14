module JacobiElliptic
using StaticArrays
using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(TYPEDSIGNATURES)
                                         $(DOCSTRING)
                                         """


export Fukushima, Carlson, AGM
export am, sn, cn, dn, nn, sd, dd, nd, sc, cc, dc, nc, ss, cs, ds, ns # Do not export cd to avoid name clash with Base.cd
export ellipj, ellipke
export E, F, K, Pi, J
export default_algorithm, set_default_algorithm!
include(joinpath(@__DIR__, "Fukushima.jl"))
include(joinpath(@__DIR__, "Carlson.jl"))
include(joinpath(@__DIR__, "AGM.jl"))

"""
Supertype for the algorithms used to evaluate the elliptic integrals and Jacobi elliptic
functions. Pass an instance as the first argument to select one explicitly, e.g.
`F(Carlson(), φ, m)`.
"""
abstract type AbstractAlgorithm end

"""
Toshio Fukushima's algorithms. Available for every function in [`JacobiElliptic`](@ref)
except the complete integrals' AGM forms, and the only algorithm providing `asn` and `acn`.
"""
struct Fukushima <: AbstractAlgorithm end

"""
B. C. Carlson's symmetric-form algorithms, built on the `DRF`, `DRD` and `DRJ` duplication
routines. This is the default for every function other than the complete `K` and `E`.
"""
struct Carlson <: AbstractAlgorithm end

"""
Arithmetic-geometric mean algorithm. Defined only for the complete integrals `K(m)` and
`E(m)`, for which it is the default.
"""
struct AGM <: AbstractAlgorithm end


func_syms = [
    :E,
    :F,
    :K,
    :Pi,
    :J,
    :sn,
    :cn,
    :dn,
    :nn,
    :sd,
    :dd,
    :nd,
    :sc,
    :cc,
    :dc,
    :nc,
    :ss,
    :cs,
    :ds,
    :ns,
    :cd,
    :am,
]

# Creates function func(::alg, args...) = algAlg.sym(args...)
for alg in [:Fukushima, :Carlson]
    for func in func_syms
        f = Symbol(alg, :Alg)
        @eval $func(::$alg, args...) = $f.$func(args...)
    end
end

# AGM implements only the complete integrals; everything else delegates to Carlson. The
# one-argument K and E below are more specific than this fallback, so they win for the
# complete forms while E(::AGM, φ, m) lands here.
for func in func_syms
    @eval $func(::AGM, args...) = CarlsonAlg.$func(args...)
end

K(::AGM, m) = ArithmeticGeometricMeanAlg.K(m)
E(::AGM, m) = ArithmeticGeometricMeanAlg.E(m)

"""
$(TYPEDSIGNATURES)

Return the algorithm currently used when none is passed explicitly. Defaults to
[`AGM`](@ref), which supplies the complete `K` and `E` and delegates everything else to
[`Carlson`](@ref). Change it with [`set_default_algorithm!`](@ref).
"""
default_algorithm() = AGM()

"""
$(TYPEDSIGNATURES)

Set the process-wide default algorithm and return it. Accepts an instance or the type
itself, so `set_default_algorithm!(Fukushima())` and `set_default_algorithm!(Fukushima)`
are equivalent.

Intended for configuration at startup.
"""
function set_default_algorithm!(alg::AbstractAlgorithm)
    @eval default_algorithm() = $alg
    return alg
end
set_default_algorithm!(::Type{A}) where {A<:AbstractAlgorithm} = set_default_algorithm!(A())

_effective_default(@nospecialize(f), alg::AbstractAlgorithm) = alg
_effective_default(@nospecialize(f), ::AGM) = Carlson()
_effective_default(::typeof(K), alg::AGM) = alg
_effective_default(::typeof(E), alg::AGM) = alg

"""
$(TYPEDSIGNATURES)

A callable binding `f` to `alg`, so that `AlgorithmCall(f, alg)(args...)` is
`f(alg, args...)`. Used by GPU backends, which need a concrete function to compile.
"""
struct AlgorithmCall{F,A<:AbstractAlgorithm} <: Function
    f::F
    alg::A
end
(call::AlgorithmCall)(args...) = call.f(call.alg, args...)

# Un-prefixed entry points route through the current default.
for func in func_syms
    @eval $func(args...) = $func(default_algorithm(), args...)
end

asn = FukushimaAlg.asn
acn = FukushimaAlg.acn

# Documentation for the entry points generated above. Each dispatches to the default
# algorithm; pass an `AbstractAlgorithm` first to override, e.g. `K(Fukushima(), m)`.
@doc """
``K(m) = \\int_0^{\\pi/2}\\frac{d\\theta}{\\sqrt{1-m\\sin^2\\theta}}``

Complete elliptic integral of the first kind. Defaults to [`AGM`](@ref).

# Arguments
- `m` : Elliptic parameter
""" K

@doc """
``E(m) = \\int_0^{\\pi/2}\\sqrt{1-m\\sin^2\\theta}\\,d\\theta`` and
``E(\\varphi|m) = \\int_0^{\\varphi}\\sqrt{1-m\\sin^2\\theta}\\,d\\theta``

Complete and incomplete elliptic integrals of the second kind. The complete form defaults
to [`AGM`](@ref), the incomplete form to [`Carlson`](@ref).

# Arguments
- `φ` : Amplitude (incomplete form only)
- `m` : Elliptic parameter
""" E

@doc """
``F(\\varphi|m) = \\int_0^{\\varphi}\\frac{d\\theta}{\\sqrt{1-m\\sin^2\\theta}}``

Incomplete elliptic integral of the first kind. Defaults to [`Carlson`](@ref).

# Arguments
- `φ` : Amplitude
- `m` : Elliptic parameter
""" F

@doc """
``\\Pi(n|m)`` and ``\\Pi(n;\\varphi|m) =
\\int_0^{\\varphi}\\frac{d\\theta}{(1-n\\sin^2\\theta)\\sqrt{1-m\\sin^2\\theta}}``

Complete and incomplete elliptic integrals of the third kind. Defaults to
[`Carlson`](@ref). For ``n > 1`` the result is the Cauchy principal value.

# Arguments
- `n` : Characteristic
- `φ` : Amplitude (incomplete form only)
- `m` : Elliptic parameter
""" Pi

@doc """
``u = F(\\varphi\\,|\\,m)`` inverted for the amplitude ``\\varphi``.

Jacobi amplitude. Defaults to [`Carlson`](@ref).

# Arguments
- `u` : Argument
- `m` : Elliptic parameter
""" am

@doc """
``J(n|m)`` and ``J(n;\\varphi|m)``, the associate elliptic integral of the third kind,
related to ``\\Pi`` by ``J = (\\Pi - F)/n``. Defaults to [`Carlson`](@ref).

# Arguments
- `n` : Characteristic
- `φ` : Amplitude (incomplete form only)
- `m` : Elliptic parameter
""" J

for (p, pn) in ((:s, "\\mathrm{sn}"), (:c, "\\mathrm{cn}"), (:d, "\\mathrm{dn}"), (:n, "1"))
    for (q, qn) in
        ((:s, "\\mathrm{sn}"), (:c, "\\mathrm{cn}"), (:d, "\\mathrm{dn}"), (:n, "1"))
        f = Symbol(p, q)
        f === :am && continue
        body = """
            ``\\mathrm{$(f)}(u|m) = \\dfrac{$(pn)(u|m)}{$(qn)(u|m)}``

            Jacobi elliptic function in Glaisher's notation, where `s`, `c`, `d` and `n`
            stand for ``\\mathrm{sn}``, ``\\mathrm{cn}``, ``\\mathrm{dn}`` and ``1``.
            Defaults to [`Carlson`](@ref).

            # Arguments
            - `u` : Argument
            - `m` : Elliptic parameter
            """
        @eval @doc $body $f
    end
end
end
