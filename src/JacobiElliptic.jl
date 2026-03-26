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

abstract type AbstractAlgorithm end
struct Fukushima <: AbstractAlgorithm end
struct Carlson <: AbstractAlgorithm end
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
end
