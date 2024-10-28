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
    sin(2*ϕ) / (4 * (1 - m) * √(1 - m * sin(ϕ)^2)) 
end

function ∂F_∂ϕ(ϕ, m)
    return 1 / √(1 - m*sin(ϕ)^2)
end

function forward(
    config::EnzymeRules.FwdConfig,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, 
    RT, 
    ϕ::Annotation{<:Real}, 
    m::Annotation{<:Real}
) 
    if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
        if EnzymeRules.width(config) == 1
            return Duplicated(
                func.val(ϕ.val, m.val), 
                ϕ isa Const ? zero(ϕ.val) : ∂F_∂ϕ(ϕ.val, m.val)*ϕ.dval, 
                m isa Const ? zero(m.val) : ∂F_∂m(ϕ.val, m.val)*m.dval
            )
        else
            return BatchDuplicated(
                func.val(ϕ.val, m.val), 
                ntuple(
                    i -> ϕ isa Const ? zero(ϕ.val) : ∂F_∂ϕ(ϕ.val, m.val)*ϕ.dval[i], Val(EnzymeRules.width(config))
                ),
                ntuple(
                    i -> m isa Const ? zero(m.val) : ∂F_∂m(ϕ.val, m.val)*m.dval[i], Val(EnzymeRules.width(config))
                ),
            )
d        end
    elseif EnzymeRules.needs_shadow(config)
        if EnzymeRules.width(config) == 1
            return DuplicatedNoNeed(
                ϕ isa Const ? zero(ϕ.val) : ∂F_∂ϕ(ϕ.val, m.val)*ϕ.dval, 
                m isa Const ? zero(m.val) : ∂F_∂m(ϕ.val, m.val)*m.dval
            )
        else
            return BatchDuplicatedNoNeed(
                ntuple(i -> ϕ isa Const ? zero(ϕ.val) : ∂F_∂ϕ(ϕ.val, m.val)*ϕ.dval[i], Val(EnzymeRules.width(config))),
                ntuple(i -> m isa Const ? zero(m.val) : ∂F_∂m(ϕ.val, m.val)*m.dval[i], Val(EnzymeRules.width(config)))
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
    m::Annotation{<:Real}
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
    m::Annotation{T}
) where T
    dϕ = if ϕ isa Const
        nothing
    elseif EnzymeRules.width(config) == 1
        ∂F_∂ϕ(ϕ.val, m.val) * dret.val
    else
        ntuple(i -> ∂F_∂ϕ(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
    end

    dm = if m isa Const
        nothing
    elseif EnzymeRules.width(config) == 1
        ∂F_∂m(ϕ.val, m.val) * dret.val
    else
        ntuple(i -> ∂F_∂m(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
    end
    return (dϕ, dm)
end


#----------------------------------------------------------------------------------------
# Elliptic E(ϕ, m)
#----------------------------------------------------------------------------------------
function ∂E_∂m(ϕ, m)
    return iszero(m) ? -π/8 : (JacobiElliptic.CarlsonAlg.E(ϕ, m)  - JacobiElliptic.CarlsonAlg.F(ϕ, m)) / (2m) 
end

function ∂E_∂ϕ(ϕ, m)
    return √(1 - m*sin(ϕ)^2)
end

function forward(
    config::EnzymeRules.FwdConfig,
    func::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, 
    RT, 
    ϕ::Annotation{<:Real}, 
    m::Annotation{<:Real}
) 
    if EnzymeRules.needs_primal(config) && EnzymeRules.needs_shadow(config)
        if EnzymeRules.width(config) == 1
            return Duplicated(
                func.val(ϕ.val, m.val), 
                ϕ isa Const ? zero(ϕ.val) : ∂E_∂ϕ(ϕ.val, m.val)*ϕ.dval, 
                m isa Const ? zero(m.val) : ∂E_∂m(ϕ.val, m.val)*m.dval
            )
        else
            return BatchDuplicated(
                func.val(ϕ.val, m.val), 
                ntuple(
                    i -> ϕ isa Const ? zero(ϕ.val) : ∂E_∂ϕ(ϕ.val, m.val)*ϕ.dval[i], Val(EnzymeRules.width(config))
                ),
                ntuple(
                    i -> m isa Const ? zero(m.val) : ∂E_∂m(ϕ.val, m.val)*m.dval[i], Val(EnzymeRules.width(config))
                ),
            )
        end
    elseif EnzymeRules.needs_shadow(config)
        if EnzymeRules.width(config) == 1
            return DuplicatedNoNeed(
                ϕ isa Const ? zero(ϕ.val) : ∂E_∂ϕ(ϕ.val, m.val)*ϕ.dval, 
                m isa Const ? zero(m.val) : ∂E_∂m(ϕ.val, m.val)*m.dval
            )
        else
            return BatchDuplicatedNoNeed(
                ntuple(i -> ϕ isa Const ? zero(ϕ.val) : ∂E_∂ϕ(ϕ.val, m.val)*ϕ.dval[i], Val(EnzymeRules.width(config))),
                ntuple(i -> m isa Const ? zero(m.val) : ∂E_∂m(ϕ.val, m.val)*m.dval[i], Val(EnzymeRules.width(config)))
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
    m::Annotation{<:Real}
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
    m::Annotation{T}
) where T
    dϕ = if ϕ isa Const
        nothing
    elseif EnzymeRules.width(config) == 1
        ∂E_∂ϕ(ϕ.val, m.val) * dret.val
    else
        ntuple(i -> ∂E_∂ϕ(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
    end

    dm = if m isa Const
        nothing
    elseif EnzymeRules.width(config) == 1
        ∂E_∂m(ϕ.val, m.val) * dret.val
    else
        ntuple(i -> ∂E_∂m(ϕ.val, m.val) * dret.val[i], Val(EnzymeRules.width(config)))
    end
    return (dϕ, dm)
end
end