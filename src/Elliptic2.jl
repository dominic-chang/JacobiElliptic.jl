module Elliptic2
# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi

# jacobi elliptic functions
export Jacobi

# matlab compatible
export ellipj, ellipke

import ForwardDiff
include("jacobi.jl")
include("slatec.jl")

_zero(T)  = zero(T)
_one(T) = one(T)
function _zero(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}  zero(V) end
function _one(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}  one(V) end

function E(phi::A, m::B) where {A,B}
    T = promote_type(A,B)
    (isnan(phi) || isnan(m)) && return T(NaN)

    #if !(0 ≤ m ≤ 1) throw(DomainError(m, "argument m not in [0,1]")) end
    if 2abs(phi) > T(π)
        phi2 = phi + T(π/2)
        return 2fld(phi2,T(π))*E(m) - _E(cos(mod(phi2,T(π))), m)
    end
    _E(sin(phi), m)
end

function _E(sinphi::A, m::B) where {A,B}
    T = promote_type(A,B)
    sinphi2 = sinphi^2
    cosphi2 = one(T) - sinphi2
    y = one(T) - m*sinphi2
    drf,ierr1 = SLATEC.DRF(cosphi2, y, _one(T))
    drd,ierr2 = SLATEC.DRD(cosphi2, y, _one(T))
    if ierr1 == ierr2 == 0
        return sinphi*(drf - m*sinphi2*drd/3)
    elseif ierr1 == ierr2 == 2
        # 2 - (1+m)*sinphi2 < tol
        return sinphi
    else
        return T(NaN)
    end
end

"""
`ellipke(m::Real)`
returns `(K(m), E(m))` for scalar `0 ≤ m ≤ 1`
"""
function ellipke(m::T) where T
    if m < one(T)
        y = 1 - m
        drf,ierr1 = SLATEC.DRF(_zero(T), y, _one(T))
        drd,ierr2 = SLATEC.DRD(_zero(T), y, _one(T))
        @assert ierr1 == 0 && ierr2 == 0
        return (drf, drf - m*drd/3)
    elseif m == 1.
        return (T(Inf), one(T))
    elseif isnan(m)
        return (T(NaN), T(NaN))
    else
        return (T(NaN), T(NaN))
        #throw(DomainError(m, "argument m not <= 1"))
    end
end

E(m) = ellipke(m)[2]

# assumes 0 ≤ m ≤ 1
function rawF(sinphi::A, m::B) where {A,B}
    T = promote_type(A,B)
    (abs(sinphi) == one(T) && m == one(T))  && return sign(sinphi)*T(Inf)
    sinphi2 = sinphi^2
    drf,ierr = SLATEC.DRF(_one(T) - sinphi2, _one(T) - m*sinphi2, _one(T))
    @assert ierr == 0
    sinphi*drf
end

function F(phi::A, m::B) where {A,B}
    T = promote_type(A,B)

    (isnan(phi) || isnan(m)) && return T(NaN) 
    (m < 0 || m > 1) && throw(DomainError(m, "argument m not in [0,1]")) 
    if abs(phi) > T(π/2)
        # Abramowitz & Stegun (17.4.3)
        phi2 = phi + T(π/2)
        return 2*fld(phi2,T(π))*K(m) - rawF(cos(mod(phi2,T(π))), m)
    end
    return rawF(sin(phi), m)
end

function K(m::T) where {T}
    if m < 1
        drf,ierr = SLATEC.DRF(_zero(T), 1 - m, _one(T))
        @assert ierr == 0
        return drf
    elseif m == 1
        return T(Inf)
    elseif isnan(m)
        return T(NaN)
    else
        #throw(DomainError(m, "argument m not <= 1"))
    #end
        return T(NaN)
    end
end

function Pi(n, phi, m)
    if isnan(n) || isnan(phi) || isnan(m) return NaN end
    if m < 0. || m > 1. throw(DomainError(m, "argument m not in [0,1]")) end
    sinphi = sin(phi)
    sinphi2 = sinphi^2
    cosphi2 = 1. - sinphi2
    y = 1. - m*sinphi2
    drf,ierr1 = SLATEC.DRF(cosphi2, y, 1.)
    drj,ierr2 = SLATEC.DRJ(cosphi2, y, 1., 1. - n*sinphi2)
    if ierr1 == 0 && ierr2 == 0
        return sinphi*(drf + n*sinphi2*drj/3)
    elseif ierr1 == 2 && ierr2 == 2
        # 2 - (1+m)*sinphi2 < tol
        return Inf
    elseif ierr1 == 0 && ierr2 == 2
        # 1 - n*sinphi2 < tol
        return Inf
    end
    NaN
end
Π = Pi

function ellipj(u, m, tol)
    phi = Jacobi.am(u, m, tol)
    s = sin(phi)
    c = cos(phi)
    d = sqrt(1. - m*s^2)
    s, c, d
end
ellipj(u, m) = ellipj(u, m, eps(Float64))

end # module
