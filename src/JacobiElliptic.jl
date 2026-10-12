module JacobiElliptic
using StaticArrays
using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(TYPEDSIGNATURES)
                                         $(DOCSTRING)
                                         """


export Fukushima, Carlson
export am, sn, cn, dn, nn, sd, dd, nd, sc, cc, dc, nc, ss, cs, ds, ns # Do not export cd to avoid name clash with Base.cd
export ellipj, ellipke
export E, F, K, Pi, J
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
]

# Creates function func(::alg, args...) = algAlg.sym(args...)
for alg in [:Fukushima, :Carlson]
    for func in func_syms
        f = Symbol(alg, :Alg)
        @eval $func(::$alg, args...) = $f.$func(args...)
    end
end

K(alg::AGM, m) = ArithmeticGeometricMeanAlg.K(m)
E(alg::AGM, m) = ArithmeticGeometricMeanAlg.E(m)


# Sets default functions to Carlson Alg
# func(args...) = CarlsonAlg.func(args...)
alg = :Carlson
for func in func_syms
    f = Symbol(alg, :Alg)
    @eval $func(args...) = $f.$func(args...)
end

K(m) = ArithmeticGeometricMeanAlg.K(m)
E(m) = ArithmeticGeometricMeanAlg.E(m)

asn = FukushimaAlg.asn
acn = FukushimaAlg.acn
am = CarlsonAlg.am

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
``J(n|m)`` and ``J(n;\\varphi|m)``, the associate elliptic integral of the third kind,
related to ``\\Pi`` by ``J = (\\Pi - F)/n``. Defaults to [`Carlson`](@ref).

# Arguments
- `n` : Characteristic
- `φ` : Amplitude (incomplete form only)
- `m` : Elliptic parameter
""" J

for (p, pn) in ((:s, "\\mathrm{sn}"), (:c, "\\mathrm{cn}"), (:d, "\\mathrm{dn}"), (:n, "1"))
    for (q, qn) in ((:s, "\\mathrm{sn}"), (:c, "\\mathrm{cn}"), (:d, "\\mathrm{dn}"), (:n, "1"))
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
