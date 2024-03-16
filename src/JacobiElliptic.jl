module JacobiElliptic
using StaticArrays, Setfield
using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(TYPEDSIGNATURES)
    $(DOCSTRING)
    """


export Fukushima, Carlson
export am, sn, cn, dn, nn, sd, cd, dd, nd, sc, cc, dc, nc, ss, cs, ds, ns
export ellipj, ellipke
export E, F, K, Pi, J
include(joinpath(@__DIR__, "Fukushima.jl"))
include(joinpath(@__DIR__, "Carlson.jl"))


abstract type AbstractAlgorithm end
struct Fukushima <: AbstractAlgorithm end
struct Carlson <: AbstractAlgorithm end

alg = Carlson()

func_syms = [:E, :F, :K, :Pi, :J, :sn, :cn, :dn, :nn, :sd, :dd, :nd, :sc, :cc, :dc, :nc, :ss, :cs, :ds, :ns, :am, :cd]
sym_list = []
for alg in [:Fukushima, :Carlson]
    for func in func_syms
        sub_sym = Expr(:., Meta.parse(string(alg)*"Alg"), Meta.parse(":($func)"))
        sym = Expr(:(=), Expr(:call, func, Expr(:(::), alg), :(args...)), Expr(:call, sub_sym, :(args...)))
        push!(sym_list, sym)
    end
end
@eval begin
    $(sym_list...)
end

default_sym_list = []
for func in func_syms
    sub_sym = Expr(:., Meta.parse(string(alg)*"Alg"), Meta.parse(":($func)"))
    sym = Expr(:(=), Expr(:call, func, :(args...)), Expr(:call, func, :(Carlson()), :(args...)))
    push!(default_sym_list, sym)
end

@eval begin
    $(default_sym_list...)
end

asn = FukushimaAlg.asn
acn = FukushimaAlg.asn

end