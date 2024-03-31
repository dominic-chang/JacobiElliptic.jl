module JacobiElliptic
using StaticArrays, Setfield
using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(TYPEDSIGNATURES)
    $(DOCSTRING)
    """


export Fukushima, Carlson
export am, sn, cn, dn, nn, sd, dd, nd, sc, cc, dc, nc, ss, cs, ds, ns # Do not export cd to avoid name clash with Base.cd
export ellipj, ellipke
export E, F, K, Pi, J
include(joinpath(@__DIR__, "Fukushima.jl"))
include(joinpath(@__DIR__, "Carlson.jl"))


abstract type AbstractAlgorithm end
struct Fukushima <: AbstractAlgorithm end
struct Carlson <: AbstractAlgorithm end


func_syms = [:E, :F, :K, :Pi, :J, :sn, :cn, :dn, :nn, :sd, :dd, :nd, :sc, :cc, :dc, :nc, :ss, :cs, :ds, :ns, :cd]
sym_list = []


# Creates function func(::alg, args...) = algAlg.sym(args...)
for alg in [:Fukushima, :Carlson]
    for func in func_syms
        f = Symbol(alg, :Alg)
        @eval $func(::$alg, args...) = $f.$func(args...)
    end
end

default_sym_list = []

# Sets default functions to Carlson Alg
# func(::args...) = CarlsonAlg.func(args...)
alg = :Carlson
for func in func_syms
    f = Symbol(alg, :Alg)
    @eval $func(args...) = $f.$func(args...)
end

asn = FukushimaAlg.asn
acn = FukushimaAlg.asn
am = CarlsonAlg.am

end