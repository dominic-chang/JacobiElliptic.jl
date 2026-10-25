module JacobiEllipticZygoteExt

using JacobiElliptic
using Zygote: @adjoint
import JacobiElliptic: CarlsonAlg, FukushimaAlg
include("derivatives.jl")


@adjoint JacobiElliptic.CarlsonAlg._sqrt(a) = JacobiElliptic.CarlsonAlg._sqrt(a),
(c̄ -> (c̄ * inv(2 * JacobiElliptic.CarlsonAlg._sqrt(a)),))

for alg in (CarlsonAlg, FukushimaAlg)
    ke =
        isdefined(alg, :ellipke) ? :(($alg).ellipke) :
        :(m -> (($alg).K(m), ($alg).E(m)))
    fe =
        isdefined(alg, :_F_E) ? :(($alg)._F_E) :
        :((ϕ, m) -> (($alg).F(ϕ, m), ($alg).E(ϕ, m)))

    @eval begin
        #------------------------------------------------------------------------------
        # Elliptic F(a, b) and E(a, b)
        #------------------------------------------------------------------------------

        @adjoint function ($alg).F(a, b)
            F, E = $fe(a, b)
            ∂ϕF, ∂mF = _∂F_E(F, E, a, b)
            return F, c̄ -> (c̄ * ∂ϕF, c̄ * ∂mF)
        end

        @adjoint function ($alg).E(a, b)
            F, E = $fe(a, b)
            _, _, ∂ϕE, ∂mE = _∂F_E(F, E, a, b)
            return E, c̄ -> (c̄ * ∂ϕE, c̄ * ∂mE)
        end

        #------------------------------------------------------------------------------
        # Elliptic Pi(a, b) and Pi(a, b, c)
        #------------------------------------------------------------------------------

        @adjoint function ($alg).Pi(a, b)
            Pi = ($alg).Pi(a, b)
            K, E = $ke(b)
            ∂n, ∂m = _∂Pi(K, E, Pi, a, b)
            return Pi, c̄ -> (c̄ * ∂n, c̄ * ∂m)
        end

        @adjoint function ($alg).Pi(a, b, c)
            Pi = ($alg).Pi(a, b, c)
            F, E = $fe(b, c)
            ∂n, ∂ϕ, ∂m = _∂Pi(F, E, Pi, a, b, c)
            return Pi, c̄ -> (c̄ * ∂n, c̄ * ∂ϕ, c̄ * ∂m)
        end

        #------------------------------------------------------------------------------
        # Associate complete elliptic integral J(n, m)
        #------------------------------------------------------------------------------

        @adjoint function ($alg).J(n, m)
            j, ∂n_j, ∂m_j =
                _complete_j_derivatives(($alg).E, ($alg).K, ($alg).Pi, ($alg).J, n, m)
            return j, c̄ -> (c̄ * ∂n_j, c̄ * ∂m_j)
        end

        #------------------------------------------------------------------------------
        # Associate incomplete elliptic integral J(n, φ, m)
        #------------------------------------------------------------------------------

        @adjoint function ($alg).J(n, φ, m)
            j, ∂n_j, ∂φ_j, ∂m_j =
                _incomplete_j_derivatives(($alg).E, ($alg).F, ($alg).Pi, ($alg).J, n, φ, m)
            return j, c̄ -> (c̄ * ∂n_j, c̄ * ∂φ_j, c̄ * ∂m_j)
        end
    end

    #----------------------------------------------------------------------------------
    # Jacobi sn(a, b), cn(a, b) and dn(a, b)
    #----------------------------------------------------------------------------------
    # Both partials of each want the whole (sn, cn, dn) triple -- the amplitude partial
    # by identity, the parameter partial through `_*_parameter_derivative` -- so one
    # `_jacobi_data` call serves them. Reaching for sn, cn, dn and am separately
    # evaluates the amplitude four times over.
    for (f, ∂u_of, ∂m_of, primal) in (
        (:sn, :_∂u_sn, :_sn_parameter_derivative, :s),
        (:cn, :_∂u_cn, :_cn_parameter_derivative, :c),
        (:dn, :_∂u_dn, :_dn_parameter_derivative, :d),
    )
        @eval @adjoint function ($alg).$f(a, b)
            s, c, d, e = _jacobi_data(($alg)._scd_am, ($alg).E, a, b)
            ∂u = $∂u_of(s, c, d, b)
            ∂m = $∂m_of(e, c / d, a, b, s, c, d)
            return $primal, c̄ -> (c̄ * ∂u, c̄ * ∂m)
        end
    end
end


end # module
