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

function forward(func::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, ::Type{<:Duplicated}, ϕ::Const, m::Duplicated) 
    return Duplicated(func.val(ϕ.val, m.val), ∂F_∂m(ϕ.val, m.val)*m.dval)
end

function forward(func::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, ::Type{<:Duplicated}, ϕ::Duplicated, m::Const) 
    return Duplicated(func.val(ϕ.val, m.val), ∂F_∂ϕ(ϕ.val, m.val)*ϕ.dval)
end

function forward(func::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, ::Type{<:Duplicated}, ϕ::Duplicated, m::Duplicated) 
    return Duplicated(func.val(ϕ.val, m.val), ∂F_∂m(ϕ.val, m.val)*m.dval + ∂F_∂ϕ(ϕ.val, m.val)*ϕ.dval)
end

function forward(func::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, ::Type{<:Const}, ϕ, m) 
    return zero(promote_type(typeof(ϕ.val), typeof(m.val)))#    func.val(ϕ.val, m.val)
end

function augmented_primal(
    config::RevConfigWidth{N},
    func::Const{typeof(JacobiElliptic.CarlsonAlg.F)},
    ::Union{Type{<:Const}, Type{<:Active}},
    ϕ,
    m
 ) where {N}
    primal = EnzymeRules.needs_primal(config) ? func.val(ϕ.val, m.val) : nothing

    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function reverse(::RevConfigWidth{1}, func::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, dret::Active, tape, ϕ::Const, m::Active) 
    dm = ∂F_∂m(ϕ.val, m.val) * dret.val
    return (nothing, dm)
end

function reverse(::RevConfigWidth{1}, ::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, dret::Active, tape, ϕ::Active, m::Const) 
    dϕ = ∂F_∂ϕ(ϕ.val, m.val) * dret.val
    return (dϕ, nothing)
end

function reverse(::RevConfigWidth{N}, ::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, dret::Active, tape, ϕ::Active, m::Active) where N
    ϕval = ϕ.val
    mval = m.val
    dm = ∂F_∂m(ϕval, mval) * dret.val
    dϕ = ∂F_∂ϕ(ϕval, mval) * dret.val
    return (dϕ, dm)
end

function reverse(::RevConfigWidth, ::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, ::Type{<:Const}, tape, ϕ::Union{Const, Duplicated}, m::Active)
    return (nothing, zero(m.val))
end

function reverse(::RevConfigWidth, ::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, ::Type{<:Const}, tape, ϕ::Active, m::Union{Const, Duplicated}) 
    return (zero(ϕ.val), nothing)
end
function reverse(::RevConfigWidth, ::Const{typeof(JacobiElliptic.CarlsonAlg.F)}, ::Type{<:Const}, tape, ϕ::Active, m::Active)
    return (zero(ϕ.val), zero(m.val))
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

function forward(func::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, ::Type{<:Duplicated}, ϕ::Const, m::Duplicated) 
    return Duplicated(func.val(ϕ.val, m.val), ∂E_∂m(ϕ.val, m.val)*m.dval)
end

function forward(func::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, ::Type{<:Duplicated}, ϕ::Duplicated, m::Const) 
    return Duplicated(func.val(ϕ.val, m.val), ∂E_∂ϕ(ϕ.val, m.val)*ϕ.dval)
end

function forward(func::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, ::Type{<:Duplicated}, ϕ::Duplicated, m::Duplicated) 
    return Duplicated(func.val(ϕ.val, m.val), ∂E_∂m(ϕ.val, m.val)*m.dval + ∂E_∂ϕ(ϕ.val, m.val)*ϕ.dval)
end

function forward(func::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, ::Type{<:Const}, ϕ, m) 
    return zero(promote_type(typeof(ϕ.val), typeof(m.val)))#    func.val(ϕ.val, m.val)
end

function augmented_primal(
    config::RevConfigWidth{N},
    func::Const{typeof(JacobiElliptic.CarlsonAlg.E)},
    ::Union{Type{<:Const}, Type{<:Active}},
    ϕ,
    m
 ) where {N}
    primal = EnzymeRules.needs_primal(config) ? func.val(ϕ.val, m.val) : nothing

    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function reverse(::RevConfigWidth{1}, func::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, dret::Active, tape, ϕ::Const, m::Active) 
    dm = ∂E_∂m(ϕ.val, m.val) * dret.val
    return (nothing, dm)
end

function reverse(::RevConfigWidth{1}, ::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, dret::Active, tape, ϕ::Active, m::Const) 
    dϕ = ∂E_∂ϕ(ϕ.val, m.val) * dret.val
    return (dϕ, nothing)
end

function reverse(::RevConfigWidth{N}, ::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, dret::Active, tape, ϕ::Active, m::Active) where N
    ϕval = ϕ.val
    mval = m.val
    dm = ∂E_∂m(ϕval, mval) * dret.val
    dϕ = ∂E_∂ϕ(ϕval, mval) * dret.val
    return (dϕ, dm)
end

function reverse(::RevConfigWidth, ::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, ::Type{<:Const}, tape, ϕ::Union{Const, Duplicated}, m::Active)
    return (nothing, zero(m.val))
end

function reverse(::RevConfigWidth, ::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, ::Type{<:Const}, tape, ϕ::Active, m::Union{Const, Duplicated}) 
    return (zero(ϕ.val), nothing)
end
function reverse(::RevConfigWidth, ::Const{typeof(JacobiElliptic.CarlsonAlg.E)}, ::Type{<:Const}, tape, ϕ::Active, m::Active)
    return (zero(ϕ.val), zero(m.val))
end



end # module