module JacobiEllipticForwardDiffExt

using JacobiElliptic
using ForwardDiff
import JacobiElliptic: ArithmeticGeometricMeanAlg, CarlsonAlg, FukushimaAlg
include("derivatives.jl")

JacobiElliptic.FukushimaAlg._dup_steps(u::ForwardDiff.Dual) =
    JacobiElliptic.FukushimaAlg._dup_steps(ForwardDiff.value(u))

function JacobiElliptic.CarlsonAlg._zero(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    zero(V)
end

function JacobiElliptic.CarlsonAlg._one(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    one(V)
end

_types = [ForwardDiff.Dual, Real]
for (T1, T2) in zip(_types, _types)
    function JacobiElliptic.CarlsonAlg._isequals(a::T1, b::T2)
        return ForwardDiff.value(a) == ForwardDiff.value(b)
    end
end

# ForwardDiff.DiffRules.@define_diffrule JacobiElliptic.CarlsonAlg._sqrt(x) = :(inv(2 * JacobiElliptic.CarlsonAlg._sqrt($x)))
function JacobiElliptic.CarlsonAlg._sqrt(x::ForwardDiff.Dual{T}) where {T}
    xval = x.value
    fval = JacobiElliptic.CarlsonAlg._sqrt(xval)
    ∂xf = inv(2 * fval)
    return ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
end

for alg in (ArithmeticGeometricMeanAlg, CarlsonAlg, FukushimaAlg)
    # Reduces to the separate kernels if the fused kernels are not defined.
    ke =
        isdefined(alg, :ellipke) ? :(($alg).ellipke) :
        :(m -> (($alg).K(m), ($alg).E(m)))
    fe =
        isdefined(alg, :_F_E) ? :(($alg)._F_E) :
        :((ϕ, m) -> (($alg).F(ϕ, m), ($alg).E(ϕ, m)))

    #----------------------------------------------------------------------------------
    # Complete K(m) and E(m)
    #----------------------------------------------------------------------------------
    for (f, primal, ∂) in ((:K, :K, :∂K), (:E, :E, :∂E))
        name = :(($alg).$f)
        @eval function $name(x::ForwardDiff.Dual{T}) where {T}
            xval = ForwardDiff.value(x)
            K, E = $ke(xval)
            ∂K, ∂E = _∂ellipke(K, E, xval)
            return ForwardDiff.Dual{T}($primal, $∂ * ForwardDiff.partials(x))
        end
    end

    family = (:F, :Pi, :J, :am, :sn, :cn, :dn, :sc, :_scd_am)
    if !all(s -> isdefined(alg, s), family)
        missing_names = filter(s -> !isdefined(alg, s), family)
        @assert !any(s -> isdefined(alg, s), family) "$alg is missing $missing_names; \
            the rules below need the whole incomplete and Jacobi family or none of it"
        continue
    end

    #----------------------------------------------------------------------------------
    # Incomplete F(ϕ, m) and E(ϕ, m)
    #----------------------------------------------------------------------------------
    for (f, primal, ∂ϕ_of, ∂m_of) in ((:F, :F, :∂ϕF, :∂mF), (:E, :E, :∂ϕE, :∂mE))
        name = :(($alg).$f)

        @eval function $name(x::ForwardDiff.Dual{T}, y::U) where {T,U}
            xval = ForwardDiff.value(x)
            fval = $name(xval, y)
            ∂xf = $(f === :F ? :(inv(√(1 - y * sin(xval)^2))) : :(√(1 - y * sin(xval)^2)))
            return ForwardDiff.Dual{T}(fval, ∂xf * ForwardDiff.partials(x))
        end

        @eval function $name(x::U, y::ForwardDiff.Dual{T}) where {T,U}
            yval = ForwardDiff.value(y)
            F, E = $fe(x, yval)
            ∂ϕF, ∂mF, ∂ϕE, ∂mE = _∂F_E(F, E, x, yval)
            return ForwardDiff.Dual{T}($primal, $∂m_of * ForwardDiff.partials(y))
        end

        @eval function $name(
            x::ForwardDiff.Dual{T},
            y::ForwardDiff.Dual{T},
        ) where {T}
            xval = ForwardDiff.value(x)
            yval = ForwardDiff.value(y)
            F, E = $fe(xval, yval)
            ∂ϕF, ∂mF, ∂ϕE, ∂mE = _∂F_E(F, E, xval, yval)
            return ForwardDiff.Dual{T}(
                $primal,
                $∂ϕ_of * ForwardDiff.partials(x) + $∂m_of * ForwardDiff.partials(y),
            )
        end
    end

    #----------------------------------------------------------------------------------
    # Associate complete elliptic integral J(n, m)
    #----------------------------------------------------------------------------------
    for mask ∈ 1:3
        n_is_dual = mask & 1 != 0
        m_is_dual = mask & 2 != 0
        n_type = n_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        m_type = m_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        n_value = n_is_dual ? :(ForwardDiff.value(n)) : :n
        m_value = m_is_dual ? :(ForwardDiff.value(m)) : :m
        tangent_terms = Any[]
        n_is_dual && push!(tangent_terms, :(∂n_j * ForwardDiff.partials(n)))
        m_is_dual && push!(tangent_terms, :(∂m_j * ForwardDiff.partials(m)))
        tangent = reduce((a, b) -> :($a + $b), tangent_terms)

        @eval function ($alg).J(n::$n_type, m::$m_type) where {T}
            n_value = $n_value
            m_value = $m_value
            j, ∂n_j, ∂m_j = _complete_j_derivatives(
                ($alg).E,
                ($alg).K,
                ($alg).Pi,
                ($alg).J,
                n_value,
                m_value,
            )
            return ForwardDiff.Dual{T}(j, $tangent)
        end
    end

    #----------------------------------------------------------------------------------
    # Associate incomplete elliptic integral J(n, φ, m)
    #----------------------------------------------------------------------------------
    for mask ∈ 1:7
        n_is_dual = mask & 1 != 0
        φ_is_dual = mask & 2 != 0
        m_is_dual = mask & 4 != 0
        n_type = n_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        φ_type = φ_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        m_type = m_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        n_value = n_is_dual ? :(ForwardDiff.value(n)) : :n
        φ_value = φ_is_dual ? :(ForwardDiff.value(φ)) : :φ
        m_value = m_is_dual ? :(ForwardDiff.value(m)) : :m
        tangent_terms = Any[]
        n_is_dual && push!(tangent_terms, :(∂n_j * ForwardDiff.partials(n)))
        φ_is_dual && push!(tangent_terms, :(∂φ_j * ForwardDiff.partials(φ)))
        m_is_dual && push!(tangent_terms, :(∂m_j * ForwardDiff.partials(m)))
        tangent = reduce((a, b) -> :($a + $b), tangent_terms)

        @eval function ($alg).J(n::$n_type, φ::$φ_type, m::$m_type) where {T}
            n_value = $n_value
            φ_value = $φ_value
            m_value = $m_value
            j, ∂n_j, ∂φ_j, ∂m_j = _incomplete_j_derivatives(
                ($alg).E,
                ($alg).F,
                ($alg).Pi,
                ($alg).J,
                n_value,
                φ_value,
                m_value,
            )
            return ForwardDiff.Dual{T}(j, $tangent)
        end
    end

    #----------------------------------------------------------------------------------
    # Jacobi sn(u, m), cn(u, m) and dn(u, m)
    #----------------------------------------------------------------------------------
    # Every one of these partials wants the whole (sn, cn, dn) triple -- the u-partials
    # by identity, the m-partials through `_*_parameter_derivative` -- so one
    # `_jacobi_data` call serves all of them. Reaching for sn, cn, dn and am separately
    # evaluates the amplitude four times over.
    for (f, ∂u_of, ∂m_of, primal) in (
        (:sn, :_∂u_sn, :_sn_parameter_derivative, :s),
        (:cn, :_∂u_cn, :_cn_parameter_derivative, :c),
        (:dn, :_∂u_dn, :_dn_parameter_derivative, :d),
    )
        name = :(($alg).$f)
        for mask ∈ 1:3
            u_is_dual = mask & 1 != 0
            m_is_dual = mask & 2 != 0
            u_type = u_is_dual ? :(ForwardDiff.Dual{T}) : :U
            m_type = m_is_dual ? :(ForwardDiff.Dual{T}) : :U
            u_value = u_is_dual ? :(ForwardDiff.value(u)) : :u
            m_value = m_is_dual ? :(ForwardDiff.value(m)) : :m
            terms = Any[]
            u_is_dual &&
                push!(terms, :($∂u_of(s, c, d, mval) * ForwardDiff.partials(u)))
            m_is_dual && push!(
                terms,
                :($∂m_of(e, c / d, uval, mval, s, c, d) * ForwardDiff.partials(m)),
            )
            tangent = reduce((a, b) -> :($a + $b), terms)
            # `U` only appears in the signature when one of the two stays non-Dual, so
            # the method has to be assembled rather than written out literally.
            params = (u_is_dual && m_is_dual) ? (:T,) : (:T, :U)
            signature = Expr(:where, :($name(u::$u_type, m::$m_type)), params...)
            body = quote
                uval = $u_value
                mval = $m_value
                s, c, d, e =
                    _jacobi_data(($alg)._scd_am, ($alg).E, uval, mval)
                return ForwardDiff.Dual{T}($primal, $tangent)
            end

            @eval $(Expr(:function, signature, body))
        end
    end
end

end # module
