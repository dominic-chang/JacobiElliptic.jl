module JacobiEllipticEnzymeExt# Should be same name as the file (just like a normal package)

using JacobiElliptic, Enzyme
import .EnzymeRules: forward, reverse, augmented_primal
using .EnzymeRules

function forward(func::Const{typeof(JacobiElliptic.F)}, ::Type{<:Duplicated}, ϕ::Type{<:Duplicated}, k::Type{<:Const}) 
    println("forward call")
    m = k.val
    dret = JacobiElliptic.E(ϕ.val, m) / (2 * m * (1 - m)) -
    JacobiElliptic.F(ϕ.val, m) / 2 / m -
    sin(2*ϕ.val) / (4 * (1 - m) * √(1 - m * sin(ϕ.val)^2)) 

    ret = func.val(ϕ.val, m)
    return Duplicated(ret, dret* ϕ.dval)
end

function augmented_primal(config::ConfigWidth{1}, func::Const{typeof(sqrt)}, ::Type{<:Active}, x::Active{<:Complex{T}}) where {T<:Real}
    println("In custom augmented primal rule.")
    if needs_primal(config)
        primal = func.val(x.val)
    else
        primal = nothing
    end

    # Save x in tape if x will be overwritten
    if overwritten(config)[2]
        tape = copy(x.val)
    else
        tape = nothing
    end

    # Return an AugmentedReturn object with shadow = nothing
    return AugmentedReturn(primal, nothing, tape)
end

function reverse(config::ConfigWidth{1}, ::Const{typeof(sqrt)}, dret::Active, tape, x::Active{<:Complex{T}}) where {T<:Real}  
    println("In custom reverse rule.")
    # retrieve x value, either from original x or from tape if x may have been overwritten.
    xval = overwritten(config)[2] ? tape : x.val
    dx = inv(2*sqrt(xval))' * dret.val
    return (dx, )
end


end # module