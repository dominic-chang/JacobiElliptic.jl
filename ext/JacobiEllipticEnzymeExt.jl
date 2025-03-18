module JacobiEllipticEnzymeExt# Should be same name as the file (just like a normal package)

using JacobiElliptic, Enzyme
import .EnzymeRules: forward, reverse, augmented_primal
using .EnzymeRules

#----------------------------------------------------------------------------------------
# Elliptic F(ϕ, m)
#----------------------------------------------------------------------------------------
function ∂F_∂m(ϕ, m)
    return JacobiElliptic.CarlsonAlg.E(ϕ, m) / (2 * m * (1 - m)) -
           JacobiElliptic.CarlsonAlg.F(ϕ, m) / 2 / m -
           sin(2 * ϕ) / (4 * (1 - m) * √(1 - m * sin(ϕ)^2))
end

function ∂F_∂ϕ(ϕ, m)
    return 1 / √(1 - m * sin(ϕ)^2)
end


function forward(
    # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
    # Of note, when we seed both arguments at once the tangent return is the sum of both.
    config::EnzymeRules.FwdConfig,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.F)},
    RT,
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
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

function augmented_primal(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.F)},
    ::Type,
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
    primal = EnzymeRules.needs_primal(config) ? func.val(ϕ.val, m.val) : nothing

    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function reverse(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.F)},
    dret,
    tape,
    ϕ::Annotation{T},
    m::Annotation{T},
) where {T}
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
            zero(ϕ.val)
        else
            ∂F_∂m(ϕ.val, m.val) * dret.val
        end
    else
        if dret isa Type{<:Const}
            ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
        else
            ntuple(i -> ∂F_∂m(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
        end
    end
    return (dϕ, dm)
end


#----------------------------------------------------------------------------------------
# Elliptic E(ϕ, m)
#----------------------------------------------------------------------------------------
function ∂E_∂m(ϕ, m)
    return iszero(m) ? -π / 8 :
           (JacobiElliptic.CarlsonAlg.E(ϕ, m) - JacobiElliptic.CarlsonAlg.F(ϕ, m)) / (2m)
end

function ∂E_∂ϕ(ϕ, m)
    return √(1 - m * sin(ϕ)^2)
end

function forward(
    # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
    # Of note, when we seed both arguments at once the tangent return is the sum of both.
    config::EnzymeRules.FwdConfig,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.E)},
    RT,
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
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

function augmented_primal(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.E)},
    ::Type,
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
    primal = EnzymeRules.needs_primal(config) ? func.val(ϕ.val, m.val) : nothing

    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function reverse(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.E)},
    dret,
    tape,
    ϕ::Annotation{T},
    m::Annotation{T},
) where {T}
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
# Elliptic Pi(n, m)
#----------------------------------------------------------------------------------------
function ∂Pi_∂n(n, m)
    return (
        JacobiElliptic.CarlsonAlg.E(m) +
        (m - n) * JacobiElliptic.CarlsonAlg.K(m) / n +
        (n^2 - m) * JacobiElliptic.CarlsonAlg.Pi(n, m) / n
    ) / (2 * (m - n) * (n - 1))
end

function ∂Pi_∂m(n, m)
    return (
        JacobiElliptic.CarlsonAlg.E(m) / (m - 1) +
        JacobiElliptic.CarlsonAlg.Pi(n, m) 
    ) / (2 * (n - m))
end

function forward(
    # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
    # Of note, when we seed both arguments at once the tangent return is the sum of both.
    config::EnzymeRules.FwdConfig,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.Pi)},
    RT,
    n::Annotation{<:Real},
    m::Annotation{<:Real},
)
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
                        (
                            n isa Const ? zero(n.val) :
                            ∂Pi_∂n(n.val, m.val) * n.dval[i]
                        ) +
                        (
                            m isa Const ? zero(m.val) :
                            ∂Pi_∂m(n.val, m.val) * m.dval[i]
                        ),
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

function augmented_primal(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.Pi)},
    ::Type,
    n::Annotation{<:Real},
    m::Annotation{<:Real},
)
    primal = EnzymeRules.needs_primal(config) ? func.val(n.val, m.val) : nothing

    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function reverse(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.Pi)},
    dret,
    tape,
    n::Annotation{T},
    m::Annotation{T},
) where {T}
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
            ntuple(
                i -> ∂Pi_∂n(n.val, m.val) * dret.val[i],
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
            ∂Pi_∂m(n.val, m.val) * dret.val
        end
    else
        if dret isa Type{<:Const}
            ntuple(i -> zero(m.val), Val(EnzymeRules.width(config)))
        else
            ntuple(
                i -> ∂Pi_∂m(n.val, m.val) * dret.val[i],
                Val(EnzymeRules.width(config)),
            )
        end
    end
    return (dn, dm)
end



#----------------------------------------------------------------------------------------
# Elliptic Pi(n, ϕ, m)
#----------------------------------------------------------------------------------------
function ∂Pi_∂n(n, ϕ, m)
    return (
        JacobiElliptic.CarlsonAlg.E(ϕ, m) +
        (m - n) * JacobiElliptic.CarlsonAlg.F(ϕ, m) / n +
        (n^2 - m) * JacobiElliptic.CarlsonAlg.Pi(n, ϕ, m) / n -
        n * √(1 - m * sin(ϕ)^2) * sin(2ϕ) / (2(1 - n * sin(ϕ)^2))
    ) / (2 * (m - n) * (n - 1))
end

function ∂Pi_∂m(n, ϕ, m)
    return (
        JacobiElliptic.CarlsonAlg.E(ϕ, m) / (m - 1) +
        JacobiElliptic.CarlsonAlg.Pi(n, ϕ, m) -
        m * sin(2 * ϕ) / (2 * (m - 1) * √(1 - m * sin(ϕ)^2))
    ) / (2 * (n - m))
end

function ∂Pi_∂ϕ(n, ϕ, m)
    return 1 / (√(1 - m * sin(ϕ)^2) * (1 - n * sin(ϕ)^2))
end


function forward(
    # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
    # Of note, when we seed both arguments at once the tangent return is the sum of both.
    config::EnzymeRules.FwdConfig,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.Pi)},
    RT,
    n::Annotation{<:Real},
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
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
                    (n isa Const ? zero(n.val) : ∂Pi_∂n(n.val, ϕ.val, m.val) * n.dval[i]) +
                    (ϕ isa Const ? zero(ϕ.val) : ∂Pi_∂ϕ(n.val, ϕ.val, m.val) * ϕ.dval[i]) +
                    (m isa Const ? zero(m.val) : ∂Pi_∂m(n.val, ϕ.val, m.val) * m.dval[i]),
                Val(EnzymeRules.width(config)),
            )
        end
    elseif EnzymeRules.needs_primal(config)
        return func.val(n.val, ϕ.val, m.val)
    else
        return nothing
    end
end

function augmented_primal(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.Pi)},
    ::Type,
    n::Annotation{<:Real},
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
    primal = EnzymeRules.needs_primal(config) ? func.val(n.val, ϕ.val, m.val) : nothing

    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function reverse(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.Pi)},
    dret,
    tape,
    n::Annotation{T},
    ϕ::Annotation{T},
    m::Annotation{T},
) where {T}
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

function ∂cn_∂m(ϕ, m)
    s = JacobiElliptic.sn(ϕ, m)
    a = JacobiElliptic.CarlsonAlg.am(ϕ, m)
    d = JacobiElliptic.dn(ϕ, m)
    c = JacobiElliptic.cn(ϕ, m)
    e = JacobiElliptic.CarlsonAlg.E(a, m)
    return inv(2m * (1 - m)) * d * s * ((m - 1) * ϕ - m * (c / d) * s + e)

end

function ∂cn_∂ϕ(ϕ, m)
    return -JacobiElliptic.dn(ϕ, m) * JacobiElliptic.sn(ϕ, m)
end

function forward(
    # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
    # Of note, when we seed both arguments at once the tangent return is the sum of both.
    config::EnzymeRules.FwdConfig,
    func::Const{typeof(JacobiElliptic.cn)},
    RT,
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
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

function augmented_primal(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.cn)},
    ::Type,
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
    primal = EnzymeRules.needs_primal(config) ? func.val(ϕ.val, m.val) : nothing

    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function reverse(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.cn)},
    dret,
    tape,
    ϕ::Annotation{T},
    m::Annotation{T},
) where {T}
    dϕ = if ϕ isa Const
        nothing
    elseif EnzymeRules.width(config) == 1
        if dret isa Type{<:Const}
            zero(ϕ.val)
        else
            ∂cn_∂ϕ(ϕ.val, m.val) * dret.val
        end
    else
        if dret isa Type{<:Const}
            ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
        else
            ntuple(i -> ∂cn_∂ϕ(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
        end
    end

    dm = if m isa Const
        nothing
    elseif EnzymeRules.width(config) == 1
        if dret isa Type{<:Const}
            zero(ϕ.val)
        else
            ∂cn_∂m(ϕ.val, m.val) * dret.val
        end
    else
        if dret isa Type{<:Const}
            ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
        else
            ntuple(i -> ∂cn_∂m(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
        end
    end
    return (dϕ, dm)
end


#----------------------------------------------------------------------------------------
# Jacobi SN(ϕ, m)
#----------------------------------------------------------------------------------------

function ∂sn_∂m(ϕ, m)
    s = JacobiElliptic.sn(ϕ, m)
    a = JacobiElliptic.CarlsonAlg.am(ϕ, m)
    d = JacobiElliptic.dn(ϕ, m)
    c = JacobiElliptic.cn(ϕ, m)
    e = JacobiElliptic.CarlsonAlg.E(a, m)
    return inv(2m * (1 - m)) * d * c * ((1 - m) * ϕ + m * (c / d) * s - e)

end

function ∂sn_∂ϕ(ϕ, m)
    return JacobiElliptic.dn(ϕ, m) * JacobiElliptic.cn(ϕ, m)
end

function forward(
    # https://enzymead.github.io/Enzyme.jl/stable/#Forward-mode
    # Of note, when we seed both arguments at once the tangent return is the sum of both.
    config::EnzymeRules.FwdConfig,
    func::Const{typeof(JacobiElliptic.sn)},
    RT,
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
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

function augmented_primal(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.sn)},
    ::Type,
    ϕ::Annotation{<:Real},
    m::Annotation{<:Real},
)
    primal = EnzymeRules.needs_primal(config) ? func.val(ϕ.val, m.val) : nothing

    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function reverse(
    config::RevConfigWidth,
    func::Const{typeof(JacobiElliptic.sn)},
    dret,
    tape,
    ϕ::Annotation{T},
    m::Annotation{T},
) where {T}
    dϕ = if ϕ isa Const
        nothing
    elseif EnzymeRules.width(config) == 1
        if dret isa Type{<:Const}
            zero(ϕ.val)
        else
            ∂sn_∂ϕ(ϕ.val, m.val) * dret.val
        end
    else
        if dret isa Type{<:Const}
            ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
        else
            ntuple(i -> ∂sn_∂ϕ(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
        end
    end

    dm = if m isa Const
        nothing
    elseif EnzymeRules.width(config) == 1
        if dret isa Type{<:Const}
            zero(ϕ.val)
        else
            ∂sn_∂m(ϕ.val, m.val) * dret.val
        end
    else
        if dret isa Type{<:Const}
            ntuple(i -> zero(ϕ.val), Val(EnzymeRules.width(config)))
        else
            ntuple(i -> ∂sn_∂m(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
        end
    end
    return (dϕ, dm)
end

end
