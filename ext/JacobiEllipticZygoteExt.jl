module JacobiEllipticZygoteExt

using JacobiElliptic
using Zygote: @adjoint

#----------------------------------------------------------------------------------------
# Elliptic E(ϕ, m)
#----------------------------------------------------------------------------------------

@adjoint JacobiElliptic.CarlsonAlg.E(a, b) = JacobiElliptic.CarlsonAlg.E(a, b), (c̄ -> c̄ .* ((sqrt(1-b*sin(a)^2)), iszero(b) ? -π/8 : (JacobiElliptic.CarlsonAlg.E(a, b)-JacobiElliptic.CarlsonAlg.F(a,b))/(2b)))

end # module