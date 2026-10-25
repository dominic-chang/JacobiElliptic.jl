module JacobiEllipticEnzymeExt

using JacobiElliptic
using Enzyme
using .EnzymeRules
import JacobiElliptic: ArithmeticGeometricMeanAlg, CarlsonAlg, FukushimaAlg
import .EnzymeRules: forward, reverse, augmented_primal
include("derivatives.jl")

#--------------------------------------------------------------------------------------
# Shared scaffolding
#--------------------------------------------------------------------------------------
# Every rule below evaluates its primal and its partials exactly once, up front, and
# then hands the finished partials to these three helpers. The batch branches call
# `tangent`/`seed` once per lane, so anything computed inside those closures rather
# than before them is paid for `width(config)` times over.

# One lane's seed for an active input: `dval` at width 1, `dval[i]` in a batch.
@inline _lane(config, x, i) = EnzymeRules.width(config) == 1 ? x.dval : x.dval[i]

# This input's contribution to a lane's tangent. A `Const` input contributes an exact
# zero without ever touching `∂`, so a partial that is NaN on a path no one asked for
# -- ∂Π/∂m at n == 0, say -- cannot poison the sum through `NaN * 0`. The annotation is
# a type, so the branch folds away at compile time.
@inline _dx(config, ∂, x, i) = x isa Const ? zero(x.val) : ∂ * _lane(config, x, i)

# Assemble whichever forward return the config asked for, with `tangent(i)` supplying
# lane `i`.
@inline function _fwd(config, primal, tangent)
    if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
        if EnzymeRules.width(config) == 1
            return Duplicated(primal, tangent(1))
        end
        return BatchDuplicated(primal, ntuple(tangent, Val(EnzymeRules.width(config))))
    elseif EnzymeRules.needs_shadow(config)
        if EnzymeRules.width(config) == 1
            return tangent(1)
        end
        return ntuple(tangent, Val(EnzymeRules.width(config)))
    elseif EnzymeRules.needs_primal(config)
        return primal
    end
    return nothing
end

# Scatter the return adjoint onto one input. `seed(d)` turns the adjoint `d` into this
# input's contribution and closes over partials that are already computed.
@inline function _rev(config, dret, x, seed)
    x isa Const && return nothing
    if EnzymeRules.width(config) == 1
        return dret isa Type{<:Const} ? zero(x.val) : seed(dret.val)
    end
    return dret isa Type{<:Const} ?
           ntuple(i -> zero(x.val), Val(EnzymeRules.width(config))) :
           ntuple(i -> seed(dret.val[i]), Val(EnzymeRules.width(config)))
end

for alg in (ArithmeticGeometricMeanAlg, CarlsonAlg, FukushimaAlg)
    ke = isdefined(alg, :ellipke) ? :(($alg).ellipke(m.val)) :
         :((($alg).K(m.val), ($alg).E(m.val)))
    fe = isdefined(alg, :_F_E) ? :(($alg)._F_E(ϕ.val, m.val)) :
         :((($alg).F(ϕ.val, m.val), ($alg).E(ϕ.val, m.val)))

    if isdefined(alg, :ellipke)
        #------------------------------------------------------------------------------
        # ellipke(m) -> (K(m), E(m))
        #------------------------------------------------------------------------------
        @eval function forward(
            config::EnzymeRules.FwdConfig,
            func::Const{typeof(($alg).ellipke)},
            RT,
            m::Annotation{<:Real},
        )
            K, E = func.val(m.val)
            ∂K, ∂E = _∂ellipke(K, E, m.val)
            function tangent(i)
                m isa Const && return (zero(K), zero(E))
                dm = _lane(config, m, i)
                return (∂K * dm, ∂E * dm)
            end
            return _fwd(config, (K, E), tangent)
        end

        @eval function augmented_primal(
            config::RevConfigWidth,
            func::Const{typeof(($alg).ellipke)},
            ::Type,
            m::Annotation{<:Real},
        )
            K, E = func.val(m.val)
            primal = EnzymeRules.needs_primal(config) ? (K, E) : nothing
            m isa Const && return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
            return EnzymeRules.AugmentedReturn(
                primal,
                nothing,
                _∂ellipke(K, E, m.val),
            )
        end

        @eval function reverse(
            config::RevConfigWidth,
            func::Const{typeof(($alg).ellipke)},
            dret,
            tape,
            m::Annotation{T},
        ) where {T}
            m isa Const && return (nothing,)
            ∂K, ∂E = tape
            # the adjoint arrives as the pair (K̄, Ē); both paths feed the one input
            return (_rev(config, dret, m, d -> ∂K * d[1] + ∂E * d[2]),)
        end
    end

    #----------------------------------------------------------------------------------
    # Elliptic K(m)
    #----------------------------------------------------------------------------------
    # https://enzymead.github.io/Enzyme.jl/stable/# Forward-mode
    # Of note, when we seed both arguments at once the tangent return is the sum of both.
    @eval function forward(
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).K)},
        RT,
        m::Annotation{<:Real},
    )
        K, E = $ke
        ∂K = _∂ellipke(K, E, m.val)[1]
        return _fwd(config, K, i -> _dx(config, ∂K, m, i))
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).K)},
        ::Type,
        m::Annotation{<:Real},
    )
        if m isa Const
            primal = EnzymeRules.needs_primal(config) ? ($alg).K(m.val) : nothing
            return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
        end
        K, E = $ke
        primal = EnzymeRules.needs_primal(config) ? K : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, _∂ellipke(K, E, m.val)[1])
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).K)},
        dret,
        tape,
        m::Annotation{T},
    ) where {T}
        m isa Const && return (nothing,)
        ∂K = tape
        return (_rev(config, dret, m, d -> ∂K * d),)
    end

    #----------------------------------------------------------------------------------
    # Elliptic E(m)
    #----------------------------------------------------------------------------------
    @eval function forward(
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).E)},
        RT,
        m::Annotation{<:Real},
    )
        K, E = $ke
        ∂E = _∂ellipke(K, E, m.val)[2]
        return _fwd(config, E, i -> _dx(config, ∂E, m, i))
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).E)},
        ::Type,
        m::Annotation{<:Real},
    )
        if m isa Const
            primal = EnzymeRules.needs_primal(config) ? ($alg).E(m.val) : nothing
            return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
        end
        K, E = $ke
        primal = EnzymeRules.needs_primal(config) ? E : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, _∂ellipke(K, E, m.val)[2])
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).E)},
        dret,
        tape,
        m::Annotation{T},
    ) where {T}
        m isa Const && return (nothing,)
        ∂E = tape
        return (_rev(config, dret, m, d -> ∂E * d),)
    end

    family = (:F, :Pi, :J, :am, :sn, :cn, :dn, :sc, :_scd_am)
    if !all(s -> isdefined(alg, s), family)
        missing_names = filter(s -> !isdefined(alg, s), family)
        @assert !any(s -> isdefined(alg, s), family) "$alg is missing $missing_names; \
            the rules below need the whole incomplete and Jacobi family or none of it"
        continue
    end

    if isdefined(alg, :_F_E)
        #------------------------------------------------------------------------------
        # _F_E(ϕ, m) -> (F(ϕ|m), E(ϕ|m))
        #------------------------------------------------------------------------------
        @eval function forward(
            config::EnzymeRules.FwdConfig,
            func::Const{typeof(($alg)._F_E)},
            RT,
            ϕ::Annotation{<:Real},
            m::Annotation{<:Real},
        )
            F, E = func.val(ϕ.val, m.val)
            ∂ϕF, ∂mF, ∂ϕE, ∂mE = _∂F_E(F, E, ϕ.val, m.val)
            tangent(i) = (
                _dx(config, ∂ϕF, ϕ, i) + _dx(config, ∂mF, m, i),
                _dx(config, ∂ϕE, ϕ, i) + _dx(config, ∂mE, m, i),
            )
            return _fwd(config, (F, E), tangent)
        end

        @eval function augmented_primal(
            config::RevConfigWidth,
            func::Const{typeof(($alg)._F_E)},
            ::Type,
            ϕ::Annotation{<:Real},
            m::Annotation{<:Real},
        )
            F, E = func.val(ϕ.val, m.val)
            primal = EnzymeRules.needs_primal(config) ? (F, E) : nothing
            if ϕ isa Const && m isa Const
                return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
            end
            return EnzymeRules.AugmentedReturn(
                primal,
                nothing,
                _∂F_E(F, E, ϕ.val, m.val),
            )
        end

        @eval function reverse(
            config::RevConfigWidth,
            func::Const{typeof(($alg)._F_E)},
            dret,
            tape,
            ϕ::Annotation{<:Real},
            m::Annotation{<:Real},
        )
            ϕ isa Const && m isa Const && return (nothing, nothing)
            ∂ϕF, ∂mF, ∂ϕE, ∂mE = tape
            # the adjoint arrives as the pair (F̄, Ē); each input collects both paths
            return (
                _rev(config, dret, ϕ, d -> ∂ϕF * d[1] + ∂ϕE * d[2]),
                _rev(config, dret, m, d -> ∂mF * d[1] + ∂mE * d[2]),
            )
        end
    end

    #----------------------------------------------------------------------------------
    # Elliptic F(ϕ, m)
    #----------------------------------------------------------------------------------
    @eval function forward(
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).F)},
        RT,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        F, E = $fe
        ∂ϕF, ∂mF = _∂F_E(F, E, ϕ.val, m.val)
        return _fwd(
            config,
            F,
            i -> _dx(config, ∂ϕF, ϕ, i) + _dx(config, ∂mF, m, i),
        )
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).F)},
        ::Type,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        if ϕ isa Const && m isa Const
            primal = EnzymeRules.needs_primal(config) ? ($alg).F(ϕ.val, m.val) : nothing
            return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
        end
        F, E = $fe
        primal = EnzymeRules.needs_primal(config) ? F : nothing
        ∂ϕF, ∂mF = _∂F_E(F, E, ϕ.val, m.val)
        return EnzymeRules.AugmentedReturn(primal, nothing, (∂ϕF, ∂mF))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).F)},
        dret,
        tape,
        ϕ::Annotation{T},
        m::Annotation{T},
    ) where {T}
        ϕ isa Const && m isa Const && return (nothing, nothing)
        ∂ϕF, ∂mF = tape
        return (
            _rev(config, dret, ϕ, d -> ∂ϕF * d),
            _rev(config, dret, m, d -> ∂mF * d),
        )
    end

    #----------------------------------------------------------------------------------
    # Elliptic E(ϕ, m)
    #----------------------------------------------------------------------------------
    @eval function forward(
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).E)},
        RT,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        F, E = $fe
        _, _, ∂ϕE, ∂mE = _∂F_E(F, E, ϕ.val, m.val)
        return _fwd(
            config,
            E,
            i -> _dx(config, ∂ϕE, ϕ, i) + _dx(config, ∂mE, m, i),
        )
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).E)},
        ::Type,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        if ϕ isa Const && m isa Const
            primal = EnzymeRules.needs_primal(config) ? ($alg).E(ϕ.val, m.val) : nothing
            return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
        end
        F, E = $fe
        primal = EnzymeRules.needs_primal(config) ? E : nothing
        _, _, ∂ϕE, ∂mE = _∂F_E(F, E, ϕ.val, m.val)
        return EnzymeRules.AugmentedReturn(primal, nothing, (∂ϕE, ∂mE))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).E)},
        dret,
        tape,
        ϕ::Annotation{T},
        m::Annotation{T},
    ) where {T}
        ϕ isa Const && m isa Const && return (nothing, nothing)
        ∂ϕE, ∂mE = tape
        return (
            _rev(config, dret, ϕ, d -> ∂ϕE * d),
            _rev(config, dret, m, d -> ∂mE * d),
        )
    end

    #----------------------------------------------------------------------------------
    # Elliptic Pi(n, m)
    #----------------------------------------------------------------------------------
    @eval function forward(
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).Pi)},
        RT,
        n::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        Pi = func.val(n.val, m.val)
        K, E = $ke
        ∂n, ∂m = _∂Pi(K, E, Pi, n.val, m.val)
        return _fwd(
            config,
            Pi,
            i -> _dx(config, ∂n, n, i) + _dx(config, ∂m, m, i),
        )
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).Pi)},
        ::Type,
        n::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        Pi = func.val(n.val, m.val)
        primal = EnzymeRules.needs_primal(config) ? Pi : nothing
        if n isa Const && m isa Const
            return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
        end
        K, E = $ke
        return EnzymeRules.AugmentedReturn(
            primal,
            nothing,
            _∂Pi(K, E, Pi, n.val, m.val),
        )
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).Pi)},
        dret,
        tape,
        n::Annotation{T},
        m::Annotation{T},
    ) where {T}
        n isa Const && m isa Const && return (nothing, nothing)
        ∂n, ∂m = tape
        return (
            _rev(config, dret, n, d -> ∂n * d),
            _rev(config, dret, m, d -> ∂m * d),
        )
    end

    #----------------------------------------------------------------------------------
    # Elliptic Pi(n, ϕ, m)
    #----------------------------------------------------------------------------------
    @eval function forward(
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).Pi)},
        RT,
        n::Annotation{<:Real},
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        Pi = func.val(n.val, ϕ.val, m.val)
        F, E = $fe
        ∂n, ∂ϕ, ∂m = _∂Pi(F, E, Pi, n.val, ϕ.val, m.val)
        return _fwd(
            config,
            Pi,
            i ->
                _dx(config, ∂n, n, i) +
                _dx(config, ∂ϕ, ϕ, i) +
                _dx(config, ∂m, m, i),
        )
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).Pi)},
        ::Type,
        n::Annotation{<:Real},
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        Pi = func.val(n.val, ϕ.val, m.val)
        primal = EnzymeRules.needs_primal(config) ? Pi : nothing
        if n isa Const && ϕ isa Const && m isa Const
            return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
        end
        F, E = $fe
        return EnzymeRules.AugmentedReturn(
            primal,
            nothing,
            _∂Pi(F, E, Pi, n.val, ϕ.val, m.val),
        )
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).Pi)},
        dret,
        tape,
        n::Annotation{T},
        ϕ::Annotation{T},
        m::Annotation{T},
    ) where {T}
        n isa Const && ϕ isa Const && m isa Const && return (nothing, nothing, nothing)
        ∂n, ∂ϕ, ∂m = tape
        return (
            _rev(config, dret, n, d -> ∂n * d),
            _rev(config, dret, ϕ, d -> ∂ϕ * d),
            _rev(config, dret, m, d -> ∂m * d),
        )
    end

    #----------------------------------------------------------------------------------
    # Associate complete elliptic integral J(n, m)
    #----------------------------------------------------------------------------------
    @eval function forward(
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).J)},
        RT,
        n::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        j, ∂n_j, ∂m_j =
            _complete_j_derivatives(($alg).E, ($alg).K, ($alg).Pi, func.val, n.val, m.val)
        return _fwd(
            config,
            j,
            i -> _dx(config, ∂n_j, n, i) + _dx(config, ∂m_j, m, i),
        )
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).J)},
        ::Type,
        n::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        if n isa Const && m isa Const
            primal = EnzymeRules.needs_primal(config) ? func.val(n.val, m.val) : nothing
            return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
        end
        j, ∂n_j, ∂m_j =
            _complete_j_derivatives(($alg).E, ($alg).K, ($alg).Pi, func.val, n.val, m.val)
        primal = EnzymeRules.needs_primal(config) ? j : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, (∂n_j, ∂m_j))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).J)},
        dret,
        tape,
        n::Annotation{T},
        m::Annotation{T},
    ) where {T}
        n isa Const && m isa Const && return (nothing, nothing)
        ∂n_j, ∂m_j = tape
        return (
            _rev(config, dret, n, d -> ∂n_j * d),
            _rev(config, dret, m, d -> ∂m_j * d),
        )
    end

    #----------------------------------------------------------------------------------
    # Associate incomplete elliptic integral J(n, φ, m)
    #----------------------------------------------------------------------------------
    @eval function forward(
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).J)},
        RT,
        n::Annotation{<:Real},
        φ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        j, ∂n_j, ∂φ_j, ∂m_j = _incomplete_j_derivatives(
            ($alg).E,
            ($alg).F,
            ($alg).Pi,
            func.val,
            n.val,
            φ.val,
            m.val,
        )
        return _fwd(
            config,
            j,
            i ->
                _dx(config, ∂n_j, n, i) +
                _dx(config, ∂φ_j, φ, i) +
                _dx(config, ∂m_j, m, i),
        )
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).J)},
        ::Type,
        n::Annotation{<:Real},
        φ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        if n isa Const && φ isa Const && m isa Const
            primal =
                EnzymeRules.needs_primal(config) ? func.val(n.val, φ.val, m.val) : nothing
            return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
        end
        j, ∂n_j, ∂φ_j, ∂m_j = _incomplete_j_derivatives(
            ($alg).E,
            ($alg).F,
            ($alg).Pi,
            func.val,
            n.val,
            φ.val,
            m.val,
        )
        primal = EnzymeRules.needs_primal(config) ? j : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, (∂n_j, ∂φ_j, ∂m_j))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).J)},
        dret,
        tape,
        n::Annotation{T},
        φ::Annotation{T},
        m::Annotation{T},
    ) where {T}
        n isa Const && φ isa Const && m isa Const && return (nothing, nothing, nothing)
        ∂n_j, ∂φ_j, ∂m_j = tape
        return (
            _rev(config, dret, n, d -> ∂n_j * d),
            _rev(config, dret, φ, d -> ∂φ_j * d),
            _rev(config, dret, m, d -> ∂m_j * d),
        )
    end

    #----------------------------------------------------------------------------------
    # Jacobi sn(ϕ, m), cn(ϕ, m), dn(ϕ, m) and sc(ϕ, m)
    #----------------------------------------------------------------------------------
    # All four share one `_jacobi_data` call: every partial here needs the whole
    # (sn, cn, dn) triple plus E at the amplitude, so there is nothing to gain from
    # splitting them and three amplitude evaluations to lose. sc is sn/cn, and the
    # triple gives it for free.
    for (f, ∂u_of, ∂m_of, primal_of) in (
        (:sn, :_∂u_sn, :_sn_parameter_derivative, :s),
        (:cn, :_∂u_cn, :_cn_parameter_derivative, :c),
        (:dn, :_∂u_dn, :_dn_parameter_derivative, :d),
        (:sc, :_∂u_sc, :_sc_parameter_derivative, :(s / c)),
    )
        fn = getproperty(alg, f)

        @eval function forward(
            config::EnzymeRules.FwdConfig,
            func::Const{typeof($fn)},
            RT,
            ϕ::Annotation{<:Real},
            m::Annotation{<:Real},
        )
            s, c, d, e = _jacobi_data(($alg)._scd_am, ($alg).E, ϕ.val, m.val)
            ∂ϕ = $∂u_of(s, c, d, m.val)
            ∂m = $∂m_of(e, c / d, ϕ.val, m.val, s, c, d)
            return _fwd(
                config,
                $primal_of,
                i -> _dx(config, ∂ϕ, ϕ, i) + _dx(config, ∂m, m, i),
            )
        end

        @eval function augmented_primal(
            config::RevConfigWidth,
            func::Const{typeof($fn)},
            ::Type,
            ϕ::Annotation{<:Real},
            m::Annotation{<:Real},
        )
            if ϕ isa Const && m isa Const
                primal =
                    EnzymeRules.needs_primal(config) ? func.val(ϕ.val, m.val) : nothing
                return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
            end
            s, c, d, e = _jacobi_data(($alg)._scd_am, ($alg).E, ϕ.val, m.val)
            primal = EnzymeRules.needs_primal(config) ? $primal_of : nothing
            return EnzymeRules.AugmentedReturn(
                primal,
                nothing,
                ($∂u_of(s, c, d, m.val), $∂m_of(e, c / d, ϕ.val, m.val, s, c, d)),
            )
        end

        @eval function reverse(
            config::RevConfigWidth,
            func::Const{typeof($fn)},
            dret,
            tape,
            ϕ::Annotation{T},
            m::Annotation{T},
        ) where {T}
            ϕ isa Const && m isa Const && return (nothing, nothing)
            ∂ϕ, ∂m = tape
            return (
                _rev(config, dret, ϕ, d -> ∂ϕ * d),
                _rev(config, dret, m, d -> ∂m * d),
            )
        end
    end
end

end
