module JacobiEllipticEnzymeExt# Should be same name as the file (just like a normal package)

using JacobiElliptic, Enzyme
import .EnzymeRules: forward, reverse, augmented_primal
using .EnzymeRules

#----------------------------------------------------------------------------------------
# Elliptic F(ϕ, m)
#----------------------------------------------------------------------------------------
function ∂F_∂m(ϕ, m)
    return JacobiElliptic.E(ϕ, m) / (2 * m * (1 - m)) -
    JacobiElliptic.F(ϕ, m) / 2 / m -
    sin(2*ϕ) / (4 * (1 - m) * √(1 - m * sin(ϕ)^2)) 
end

function ∂F_∂ϕ(ϕ, m)
    return 1 / √(1 - m*sin(ϕ)^2)
end

function forward(
    ::Const{typeof(JacobiElliptic.F)},
    ::Type,
    ϕ::Const,
    m::Duplicated
) 
    println("In custom forward rule.")

    return ∂F_∂m(ϕ.val, m.val)*m.dval
end

function forward(
    ::Const{typeof(JacobiElliptic.F)},
    ::Type,
    ϕ::Duplicated,
    m::Const
) 
    println("In custom forward rule.")

    return ∂F_∂ϕ(ϕ.val, m.val)*ϕ.dval
end

function forward(
    ::Const{typeof(JacobiElliptic.F)},
    ::Type,
    ϕ::Duplicated,
    m::Duplicated
) 
    println("In custom forward rule.")

    return ∂F_∂m(ϕ.val, m.val)*m.dval + ∂F_∂ϕ(ϕ.val, m.val)*ϕ.dval
end

function augmented_primal(
    config,
    func::Const{typeof(JacobiElliptic.F)},
    ::Type{<:Active},
    tape,
    ϕ,
    m
 ) 
    println("In custom augmented primal rule.")
    # Save x in tape if x will be overwritten
    primal = EnzymeRules.needs_primal(config) ? func(ϕ.val, m.val) : nothing

    return EnzymeRules.AugmentedReturn(primal, nothing, tape)
end


function reverse(config::ConfigWidth{N}, ::Const{JacobiElliptic.F}, dret::Active, tape, ϕ::Const, m::Active) where {N}  
    println("In custom reverse rule.")
    # retrieve x value, either from original x or from tape if x may have been overwritten.
    mval = m.val
    dm = ∂F_∂m(ϕ.val, mval) * dret.val
    return (dm, )
end

function reverse(config::ConfigWidth{N}, ::Const{JacobiElliptic.F}, dret::Active, tape, ϕ::Active, m::Const) where {N}  
    println("In custom reverse rule.")
    # retrieve x value, either from original x or from tape if x may have been overwritten.
    ϕval = EnzymeRules.overwritten(config)[2] ? tape : ϕ.val
    dϕ = ∂F_∂ϕ(ϕval, m.val) * dret.val
    return (dϕ, )
end

function reverse(config::ConfigWidth{N}, ::Const{JacobiElliptic.F}, dret::Active, tape, ϕ::Active, m::Active) where {N}  
    println("In custom reverse rule.")
    # retrieve x value, either from original x or from tape if x may have been overwritten.
    ϕval = ϕ.val
    mval = m.val
    dm = ∂F_∂m(ϕval, mval) * dret.val
    dϕ = ∂F_∂ϕ(ϕval, mval) * dret.val
    return (dϕ, dm)
end


end # module