module JacobiEllipticForwardDiffExt

using JacobiElliptic, ForwardDiff

function JacobiElliptic.CarlsonAlg._zero(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    zero(V)
end
function JacobiElliptic.CarlsonAlg._one(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    one(V)
end

#----------------------------------------------------------------------------------------
# Elliptic E(ϕ, m)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.E(x::ForwardDiff.Dual{T}, y::U) where {T, U} 
    xval = x.value

    ForwardDiff.Dual{T}(JacobiElliptic.CarlsonAlg.E(xval, y), (sqrt(1-y*sin(xval)^2)))
end

function JacobiElliptic.CarlsonAlg.E(x::U, y::ForwardDiff.Dual{T}) where {T, U}
    yval = y.value
    ForwardDiff.Dual{T}(JacobiElliptic.CarlsonAlg.E(x, yval), (JacobiElliptic.CarlsonAlg.E(x, yval)-JacobiElliptic.CarlsonAlg.F(x,yval))/(2yval))
end

function JacobiElliptic.CarlsonAlg.E(x::ForwardDiff.Dual{T}, y::ForwardDiff.Dual{T}) where T
    xval = x.value
    yval = y.value

    ∂yE = iszero(yval) ? -π/8 : (JacobiElliptic.CarlsonAlg.E(xval, yval)-JacobiElliptic.CarlsonAlg.F(xval, yval))/(2yval)
    ForwardDiff.Dual{T}(JacobiElliptic.CarlsonAlg.E(xval, yval), ForwardDiff.Partials(((sqrt(1-yval*sin(xval)^2)), ∂yE)))
end


end # module
