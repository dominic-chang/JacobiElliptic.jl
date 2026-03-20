module JacobiEllipticEnzymeExt# Should be same name as the file (just like a normal package)

using JacobiElliptic, Enzyme
import .EnzymeRules: forward, reverse, augmented_primal
using .EnzymeRules

for alg in [JacobiElliptic.CarlsonAlg, JacobiElliptic.FukushimaAlg]
    #----------------------------------------------------------------------------------------
    # Elliptic K(m)
    #----------------------------------------------------------------------------------------
    @eval function forward(
        # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
        # Of note, when we seed both arguments at once the tangent return is the sum of both.
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).K)},
        RT,
        m::Annotation{<:Real},
    )
        ∂K_∂m(m) = ($alg).E(m) / (2 * m * (1 - m)) - ($alg).K(m) / 2 / m

        if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return Duplicated(
                    func.val(m.val),
                    (m isa Const ? zero(m.val) : ∂K_∂m(m.val) * m.dval),
                )
            else
                return BatchDuplicated(
                    func.val(m.val),
                    ntuple(
                        i -> (m isa Const ? zero(m.val) : ∂K_∂m(m.val) * m.dval[i]),
                        Val(EnzymeRules.width(config)),
                    ),
                )
            end
        elseif EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return (m isa Const ? zero(m.val) : ∂K_∂m(m.val) * m.dval)
            else
                return ntuple(
                    i -> (m isa Const ? zero(m.val) : ∂K_∂m(m.val) * m.dval[i]),
                    Val(EnzymeRules.width(config)),
                )
            end
        elseif EnzymeRules.needs_primal(config)
            return func.val(m.val)
        else
            return nothing
        end
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).K)},
        ::Type,
        m::Annotation{<:Real},
    )
        K = func.val(m.val)
        E = ($alg).E(m.val)
        primal = EnzymeRules.needs_primal(config) ? K : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, (E, K))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).K)},
        dret,
        tape,
        m::Annotation{T},
    ) where {T}
        E, K = tape
        ∂K_∂m(m) = E / (2 * m * (1 - m)) - K / 2 / m

        dm = if m isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(m.val)
            else
                ∂K_∂m(m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(m.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂K_∂m(m.val) * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end
        return (dm,)
    end

    #----------------------------------------------------------------------------------------
    # Elliptic Pi(n, m)
    #----------------------------------------------------------------------------------------

    @eval function forward(
        # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
        # Of note, when we seed both arguments at once the tangent return is the sum of both.
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).Pi)},
        RT,
        n::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        E, K, Pi = ($alg).E(m.val), ($alg).K(m.val), func.val(n.val, m.val)
        ∂Pi_∂n(n, m) = begin
            if iszero(n)
                if iszero(m)
                    π / 4
                end
                (K - E) / m
            else
                (E + (m - n) * K / n + (n^2 - m) * Pi / n) / (2 * (m - n) * (n - 1))
            end
        end

        ∂Pi_∂m(n, m) = (E / (m - 1) + Pi) / (2 * (n - m))

        if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return Duplicated(
                    func.val(n.val, m.val),
                    (n isa Const ? zero(n.val) : ∂Pi_∂n(n.val, m.val) * n.dval) +
                    (m isa Const ? zero(m.val) : ∂Pi_∂m(n.val, m.val) * m.dval),
                )
            else
                return BatchDuplicated(
                    func.val(n.val, m.val),
                    ntuple(
                        i ->
                            (n isa Const ? zero(n.val) : ∂Pi_∂n(n.val, m.val) * n.dval[i]) +
                            (m isa Const ? zero(m.val) : ∂Pi_∂m(n.val, m.val) * m.dval[i]),
                        Val(EnzymeRules.width(config)),
                    ),
                )
            end
        elseif EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return (n isa Const ? zero(n.val) : ∂Pi_∂n(n.val, m.val) * n.dval) +
                       (m isa Const ? zero(m.val) : ∂Pi_∂m(n.val, m.val) * m.dval)
            else
                return ntuple(
                    i ->
                        (n isa Const ? zero(n.val) : ∂Pi_∂n(n.val, m.val) * n.dval[i]) +
                        (m isa Const ? zero(m.val) : ∂Pi_∂m(n.val, m.val) * m.dval[i]),
                    Val(EnzymeRules.width(config)),
                )
            end
        elseif EnzymeRules.needs_primal(config)
            return func.val(n.val, m.val)
        else
            return nothing
        end
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).Pi)},
        ::Type,
        n::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        Pi = func.val(n.val, m.val)
        E = ($alg).E(m.val)
        K = ($alg).K(m.val)
        primal = EnzymeRules.needs_primal(config) ? Pi : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, (E, K, Pi))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).Pi)},
        dret,
        tape,
        n::Annotation{T},
        m::Annotation{T},
    ) where {T}
        E, K, Pi = tape
        ∂Pi_∂n(n, m) = begin
            if iszero(n)
                if iszero(m)
                    π / 4
                end
                (K - E) / m
            else
                (E + (m - n) * K / n + (n^2 - m) * Pi / n) / (2 * (m - n) * (n - 1))
            end
        end

        ∂Pi_∂m(n, m) = (E / (m - 1) + Pi) / (2 * (n - m))

        dn = if n isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(n.val)
            else
                ∂Pi_∂n(n.val, m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(n.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂Pi_∂n(n.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end

        dm = if m isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(m.val)
            else
                ∂Pi_∂m(n.val, m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(m.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂Pi_∂m(n.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end
        return (dn, dm)
    end

    #----------------------------------------------------------------------------------------
    # Elliptic F(ϕ, m)
    #----------------------------------------------------------------------------------------

    @eval function forward(
        # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
        # Of note, when we seed both arguments at once the tangent return is the sum of both.
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).F)},
        RT,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        ∂F_∂m(ϕ, m) =
            ($alg).E(ϕ, m) / (2 * m * (1 - m)) - ($alg).F(ϕ, m) / 2 / m -
            sin(2 * ϕ) / (4 * (1 - m) * √(1 - m * sin(ϕ)^2))

        ∂F_∂ϕ(ϕ, m) = 1 / √(1 - m * sin(ϕ)^2)

        if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return Duplicated(
                    func.val(ϕ.val, m.val),
                    (ϕ isa Const ? zero(ϕ.val) : ∂F_∂ϕ(ϕ.val, m.val) * ϕ.dval) +
                    (m isa Const ? zero(m.val) : ∂F_∂m(ϕ.val, m.val) * m.dval),
                )
            else
                return BatchDuplicated(
                    func.val(ϕ.val, m.val),
                    ntuple(
                        i ->
                            (ϕ isa Const ? zero(ϕ.val) : ∂F_∂ϕ(ϕ.val, m.val) * ϕ.dval[i]) +
                            (m isa Const ? zero(m.val) : ∂F_∂m(ϕ.val, m.val) * m.dval[i]),
                        Val(EnzymeRules.width(config)),
                    ),
                )
            end
        elseif EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return (ϕ isa Const ? zero(ϕ.val) : ∂F_∂ϕ(ϕ.val, m.val) * ϕ.dval) +
                       (m isa Const ? zero(m.val) : ∂F_∂m(ϕ.val, m.val) * m.dval)
            else
                return ntuple(
                    i ->
                        (ϕ isa Const ? zero(ϕ.val) : ∂F_∂ϕ(ϕ.val, m.val) * ϕ.dval[i]) +
                        (m isa Const ? zero(m.val) : ∂F_∂m(ϕ.val, m.val) * m.dval[i]),
                    Val(EnzymeRules.width(config)),
                )
            end
        elseif EnzymeRules.needs_primal(config)
            return func.val(ϕ.val, m.val)
        else
            return nothing
        end
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).F)},
        ::Type,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        F = func.val(ϕ.val, m.val)
        E = ($alg).E(ϕ.val, m.val)
        sin_2ϕ = sin(2 * ϕ.val)
        sqrt_term = √(1 - m.val * sin(ϕ.val)^2)
        primal = EnzymeRules.needs_primal(config) ? F : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, (E, sin_2ϕ, sqrt_term, F))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).F)},
        dret,
        tape,
        ϕ::Annotation{T},
        m::Annotation{T},
    ) where {T}
        E, sin_2ϕ, sqrt_term, F = tape
        ∂F_∂m(ϕ, m) = E / (2 * m * (1 - m)) - F / 2 / m - sin_2ϕ / (4 * (1 - m) * sqrt_term)
        ∂F_∂ϕ(ϕ, m) = 1 / sqrt_term

        dϕ = if ϕ isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(ϕ.val)
            else
                ∂F_∂ϕ(ϕ.val, m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂F_∂ϕ(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end

        dm = if m isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(m.val)
            else
                ∂F_∂m(ϕ.val, m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(m.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂F_∂m(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end
        return (dϕ, dm)
    end


    #----------------------------------------------------------------------------------------
    # Elliptic E(ϕ, m)
    #----------------------------------------------------------------------------------------

    @eval function forward(
        # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
        # Of note, when we seed both arguments at once the tangent return is the sum of both.
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).E)},
        RT,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        ∂E_∂m(ϕ, m) =
            iszero(m) ? (sin(2 * ϕ) - 2 * ϕ) / 8 : (($alg).E(ϕ, m) - ($alg).F(ϕ, m)) / (2m)

        ∂E_∂ϕ(ϕ, m) = √(1 - m * sin(ϕ)^2)

        if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return Duplicated(
                    func.val(ϕ.val, m.val),
                    (ϕ isa Const ? zero(ϕ.val) : ∂E_∂ϕ(ϕ.val, m.val) * ϕ.dval) +
                    (m isa Const ? zero(m.val) : ∂E_∂m(ϕ.val, m.val) * m.dval),
                )
            else
                return BatchDuplicated(
                    func.val(ϕ.val, m.val),
                    ntuple(
                        i ->
                            (ϕ isa Const ? zero(ϕ.val) : ∂E_∂ϕ(ϕ.val, m.val) * ϕ.dval[i]) +
                            (m isa Const ? zero(m.val) : ∂E_∂m(ϕ.val, m.val) * m.dval[i]),
                        Val(EnzymeRules.width(config)),
                    ),
                )
            end
        elseif EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return (ϕ isa Const ? zero(ϕ.val) : ∂E_∂ϕ(ϕ.val, m.val) * ϕ.dval) +
                       (m isa Const ? zero(m.val) : ∂E_∂m(ϕ.val, m.val) * m.dval)

            else
                return ntuple(
                    i ->
                        (ϕ isa Const ? zero(ϕ.val) : ∂E_∂ϕ(ϕ.val, m.val) * ϕ.dval[i]) +
                        (m isa Const ? zero(m.val) : ∂E_∂m(ϕ.val, m.val) * m.dval[i]),
                    Val(EnzymeRules.width(config)),
                )

            end
        elseif EnzymeRules.needs_primal(config)
            return func.val(ϕ.val, m.val)
        else
            return nothing
        end
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).E)},
        ::Type,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        E = func.val(ϕ.val, m.val)
        F = ($alg).F(ϕ.val, m.val)
        primal = EnzymeRules.needs_primal(config) ? E : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, (E, F))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).E)},
        dret,
        tape,
        ϕ::Annotation{T},
        m::Annotation{T},
    ) where {T}
        E, F = tape
        ∂E_∂m(ϕ, m) = iszero(m) ? (sin(2 * ϕ) - 2 * ϕ) / 8 : (E - F) / (2m)
        ∂E_∂ϕ(ϕ, m) = √(1 - m * sin(ϕ)^2)

        dϕ = if ϕ isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(ϕ.val)
            else
                ∂E_∂ϕ(ϕ.val, m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂E_∂ϕ(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end

        dm = if m isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(ϕ.val)
            else
                ∂E_∂m(ϕ.val, m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂E_∂m(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end
        return (dϕ, dm)
    end


    #----------------------------------------------------------------------------------------
    # Elliptic Pi(n, ϕ, m)
    #----------------------------------------------------------------------------------------
    @eval function forward(
        # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
        # Of note, when we seed both arguments at once the tangent return is the sum of both.
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).Pi)},
        RT,
        n::Annotation{<:Real},
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        E, F, Pi =
            ($alg).E(ϕ.val, m.val), ($alg).F(ϕ.val, m.val), func.val(n.val, ϕ.val, m.val)
        ∂Pi_∂n(n, ϕ, m) = begin
            if iszero(n)
                if iszero(m)
                    (2 * ϕ - sin(2 * ϕ)) / 4
                end
                (F - E) / m
            else
                (
                    E + (m - n) * F / n + (n^2 - m) * Pi / n -
                    n * √(1 - m * sin(ϕ)^2) * sin(2ϕ) / (2(1 - n * sin(ϕ)^2))
                ) / (2 * (m - n) * (n - 1))
            end
        end

        ∂Pi_∂m(n, ϕ, m) =
            (E / (m - 1) + Pi - m * sin(2 * ϕ) / (2 * (m - 1) * √(1 - m * sin(ϕ)^2))) /
            (2 * (n - m))

        ∂Pi_∂ϕ(n, ϕ, m) = 1 / (√(1 - m * sin(ϕ)^2) * (1 - n * sin(ϕ)^2))

        if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return Duplicated(
                    func.val(n.val, ϕ.val, m.val),
                    (n isa Const ? zero(n.val) : ∂Pi_∂n(n.val, ϕ.val, m.val) * n.dval) +
                    (ϕ isa Const ? zero(ϕ.val) : ∂Pi_∂ϕ(n.val, ϕ.val, m.val) * ϕ.dval) +
                    (m isa Const ? zero(m.val) : ∂Pi_∂m(n.val, ϕ.val, m.val) * m.dval),
                )
            else
                return BatchDuplicated(
                    func.val(n.val, ϕ.val, m.val),
                    ntuple(
                        i ->
                            (
                                n isa Const ? zero(n.val) :
                                ∂Pi_∂n(n.val, ϕ.val, m.val) * n.dval[i]
                            ) +
                            (
                                ϕ isa Const ? zero(ϕ.val) :
                                ∂Pi_∂ϕ(n.val, ϕ.val, m.val) * ϕ.dval[i]
                            ) +
                            (
                                m isa Const ? zero(m.val) :
                                ∂Pi_∂m(n.val, ϕ.val, m.val) * m.dval[i]
                            ),
                        Val(EnzymeRules.width(config)),
                    ),
                )
            end
        elseif EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return (n isa Const ? zero(n.val) : ∂Pi_∂n(n.val, ϕ.val, m.val) * n.dval) +
                       (ϕ isa Const ? zero(ϕ.val) : ∂Pi_∂ϕ(n.val, ϕ.val, m.val) * ϕ.dval) +
                       (m isa Const ? zero(m.val) : ∂Pi_∂m(n.val, ϕ.val, m.val) * m.dval)
            else
                return ntuple(
                    i ->
                        (
                            n isa Const ? zero(n.val) :
                            ∂Pi_∂n(n.val, ϕ.val, m.val) * n.dval[i]
                        ) +
                        (
                            ϕ isa Const ? zero(ϕ.val) :
                            ∂Pi_∂ϕ(n.val, ϕ.val, m.val) * ϕ.dval[i]
                        ) +
                        (
                            m isa Const ? zero(m.val) :
                            ∂Pi_∂m(n.val, ϕ.val, m.val) * m.dval[i]
                        ),
                    Val(EnzymeRules.width(config)),
                )
            end
        elseif EnzymeRules.needs_primal(config)
            return func.val(n.val, ϕ.val, m.val)
        else
            return nothing
        end
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
        E = ($alg).E(ϕ.val, m.val)
        F = ($alg).F(ϕ.val, m.val)
        sin_2ϕ = sin(2 * ϕ.val)
        sin_ϕ = sin(ϕ.val)
        sqrt_term = √(1 - m.val * sin_ϕ^2)
        primal = EnzymeRules.needs_primal(config) ? Pi : nothing
        return EnzymeRules.AugmentedReturn(
            primal,
            nothing,
            (E, F, sin_2ϕ, sin_ϕ, sqrt_term, Pi),
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
        E, F, sin_2ϕ, sin_ϕ, sqrt_term, Pi = tape
        ∂Pi_∂n(n, ϕ, m) = begin
            if iszero(n)
                if iszero(m)
                    (2 * ϕ - sin_2ϕ) / 4
                end
                (F - E) / m
            else
                (
                    E + (m - n) * F / n + (n^2 - m) * Pi / n -
                    n * sqrt_term * sin_2ϕ / (2 * (1 - n * sin_ϕ^2))
                ) / (2 * (m - n) * (n - 1))
            end
        end

        ∂Pi_∂m(n, ϕ, m) =
            (E / (m - 1) + Pi - m * sin_2ϕ / (2 * (m - 1) * sqrt_term)) / (2 * (n - m))

        ∂Pi_∂ϕ(n, ϕ, m) = 1 / (sqrt_term * (1 - n * sin_ϕ^2))

        dn = if n isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(n.val)
            else
                ∂Pi_∂n(n.val, ϕ.val, m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(n.val), Val(EnzymeRules.width(config)))
            else
                ntuple(
                    i -> ∂Pi_∂n(n.val, ϕ.val, m.val) * dret.val[i],
                    Val(EnzymeRules.width(config)),
                )
            end
        end

        dϕ = if ϕ isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(ϕ.val)
            else
                ∂Pi_∂ϕ(n.val, ϕ.val, m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
            else
                ntuple(
                    i -> ∂Pi_∂ϕ(n.val, ϕ.val, m.val) * dret.val[i],
                    Val(EnzymeRules.width(config)),
                )
            end
        end

        dm = if m isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(m.val)
            else
                ∂Pi_∂m(n.val, ϕ.val, m.val) * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(m.val), Val(EnzymeRules.width(config)))
            else
                ntuple(
                    i -> ∂Pi_∂m(n.val, ϕ.val, m.val) * dret.val[i],
                    Val(EnzymeRules.width(config)),
                )
            end
        end
        return (dn, dϕ, dm)
    end

    #----------------------------------------------------------------------------------------
    # Jacobi CN(ϕ, m)
    #----------------------------------------------------------------------------------------

    @eval function forward(
        # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
        # Of note, when we seed both arguments at once the tangent return is the sum of both.
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).cn)},
        RT,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        ∂cn_∂m(ϕ, m) = begin
            s = ($alg).sn(ϕ, m)
            a = ($alg).am(ϕ, m)
            d = ($alg).dn(ϕ, m)
            c = ($alg).cn(ϕ, m)
            e = ($alg).E(a, m)
            inv(2m * (1 - m)) * d * s * ((m - 1) * ϕ - m * (c / d) * s + e)

        end

        ∂cn_∂ϕ(ϕ, m) = -($alg).dn(ϕ, m) * ($alg).sn(ϕ, m)

        if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return Duplicated(
                    func.val(ϕ.val, m.val),
                    (ϕ isa Const ? zero(ϕ.val) : ∂cn_∂ϕ(ϕ.val, m.val) * ϕ.dval) +
                    (m isa Const ? zero(m.val) : ∂cn_∂m(ϕ.val, m.val) * m.dval),
                )
            else
                return BatchDuplicated(
                    func.val(ϕ.val, m.val),
                    ntuple(
                        i ->
                            (ϕ isa Const ? zero(ϕ.val) : ∂cn_∂ϕ(ϕ.val, m.val) * ϕ.dval[i]) +
                            (m isa Const ? zero(m.val) : ∂cn_∂m(ϕ.val, m.val) * m.dval[i]),
                        Val(EnzymeRules.width(config)),
                    ),
                )
            end
        elseif EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return (ϕ isa Const ? zero(ϕ.val) : ∂cn_∂ϕ(ϕ.val, m.val) * ϕ.dval) +
                       (m isa Const ? zero(m.val) : ∂cn_∂m(ϕ.val, m.val) * m.dval)

            else
                return ntuple(
                    i ->
                        (ϕ isa Const ? zero(ϕ.val) : ∂cn_∂ϕ(ϕ.val, m.val) * ϕ.dval[i]) +
                        (m isa Const ? zero(m.val) : ∂cn_∂m(ϕ.val, m.val) * m.dval[i]),
                    Val(EnzymeRules.width(config)),
                )

            end
        elseif EnzymeRules.needs_primal(config)
            return func.val(ϕ.val, m.val)
        else
            return nothing
        end
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).cn)},
        ::Type,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        c = func.val(ϕ.val, m.val)
        s = ($alg).sn(ϕ.val, m.val)
        d = ($alg).dn(ϕ.val, m.val)
        e = ($alg).E(($alg).am(ϕ.val, m.val), m.val)
        primal = EnzymeRules.needs_primal(config) ? c : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, (s, d, c, e))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).cn)},
        dret,
        tape,
        ϕ::Annotation{T},
        m::Annotation{T},
    ) where {T}
        s, d, c, e = tape
        ∂cn_∂m =
            inv(2 * m.val * (1 - m.val)) *
            d *
            s *
            ((m.val - 1) * ϕ.val - m.val * (c / d) * s + e)
        ∂cn_∂ϕ = -d * s

        dϕ = if ϕ isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(ϕ.val)
            else
                ∂cn_∂ϕ * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂cn_∂ϕ * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end

        dm = if m isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(m.val)
            else
                ∂cn_∂m * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(m.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂cn_∂m * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end
        return (dϕ, dm)
    end


    #----------------------------------------------------------------------------------------
    # Jacobi SN(ϕ, m)
    #----------------------------------------------------------------------------------------

    @eval function forward(
        # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
        # Of note, when we seed both arguments at once the tangent return is the sum of both.
        config::EnzymeRules.FwdConfig,
        func::Const{typeof(($alg).sn)},
        RT,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        ∂sn_∂m(ϕ, m) = begin
            s = ($alg).sn(ϕ, m)
            a = ($alg).am(ϕ, m)
            d = ($alg).dn(ϕ, m)
            c = ($alg).cn(ϕ, m)
            e = ($alg).E(a, m)
            inv(2m * (1 - m)) * d * c * ((1 - m) * ϕ + m * (c / d) * s - e)

        end

        ∂sn_∂ϕ(ϕ, m) = ($alg).dn(ϕ, m) * ($alg).cn(ϕ, m)


        if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return Duplicated(
                    func.val(ϕ.val, m.val),
                    (ϕ isa Const ? zero(ϕ.val) : ∂sn_∂ϕ(ϕ.val, m.val) * ϕ.dval) +
                    (m isa Const ? zero(m.val) : ∂sn_∂m(ϕ.val, m.val) * m.dval),
                )
            else
                return BatchDuplicated(
                    func.val(ϕ.val, m.val),
                    ntuple(
                        i ->
                            (ϕ isa Const ? zero(ϕ.val) : ∂sn_∂ϕ(ϕ.val, m.val) * ϕ.dval[i]) +
                            (m isa Const ? zero(m.val) : ∂sn_∂m(ϕ.val, m.val) * m.dval[i]),
                        Val(EnzymeRules.width(config)),
                    ),
                )
            end
        elseif EnzymeRules.needs_shadow(config)
            if EnzymeRules.width(config) == 1
                return (ϕ isa Const ? zero(ϕ.val) : ∂sn_∂ϕ(ϕ.val, m.val) * ϕ.dval) +
                       (m isa Const ? zero(m.val) : ∂sn_∂m(ϕ.val, m.val) * m.dval)

            else
                return ntuple(
                    i ->
                        (ϕ isa Const ? zero(ϕ.val) : ∂sn_∂ϕ(ϕ.val, m.val) * ϕ.dval[i]) +
                        (m isa Const ? zero(m.val) : ∂sn_∂m(ϕ.val, m.val) * m.dval[i]),
                    Val(EnzymeRules.width(config)),
                )

            end
        elseif EnzymeRules.needs_primal(config)
            return func.val(ϕ.val, m.val)
        else
            return nothing
        end
    end

    @eval function augmented_primal(
        config::RevConfigWidth,
        func::Const{typeof(($alg).sn)},
        ::Type,
        ϕ::Annotation{<:Real},
        m::Annotation{<:Real},
    )
        s = func.val(ϕ.val, m.val)
        d = ($alg).dn(ϕ.val, m.val)
        c = ($alg).cn(ϕ.val, m.val)
        e = ($alg).E(($alg).am(ϕ.val, m.val), m.val)
        primal = EnzymeRules.needs_primal(config) ? s : nothing
        return EnzymeRules.AugmentedReturn(primal, nothing, (s, d, c, e))
    end

    @eval function reverse(
        config::RevConfigWidth,
        func::Const{typeof(($alg).sn)},
        dret,
        tape,
        ϕ::Annotation{T},
        m::Annotation{T},
    ) where {T}
        s, d, c, e = tape
        ∂sn_∂m =
            inv(2 * m.val * (1 - m.val)) *
            d *
            c *
            ((1 - m.val) * ϕ.val + m.val * (c / d) * s - e)
        ∂sn_∂ϕ = d * c


        dϕ = if ϕ isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(ϕ.val)
            else
                ∂sn_∂ϕ * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂sn_∂ϕ * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end

        dm = if m isa Const
            nothing
        elseif EnzymeRules.width(config) == 1
            if dret isa Type{<:Const}
                zero(m.val)
            else
                ∂sn_∂m * dret.val
            end
        else
            if dret isa Type{<:Const}
                ntuple(i -> zero(m.val), Val(EnzymeRules.width(config)))
            else
                ntuple(i -> ∂sn_∂m * dret.val[i], Val(EnzymeRules.width(config)))
            end
        end
        return (dϕ, dm)
    end

end
end
