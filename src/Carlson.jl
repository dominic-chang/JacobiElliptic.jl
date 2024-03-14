module CarlsonAlg
# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi

# jacobi elliptic functions
export Jacobi

# matlab compatible
export ellipj, ellipke

import ForwardDiff
include("jacobi.jl")
include("slatec.jl")

_zero(T) = zero(T)
_one(T) = one(T)
function _zero(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    zero(V)
end
function _one(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    one(V)
end

function E(phi::A, m::B) where {A,B}
    T = promote_type(A, B)
    (isnan(phi) || isnan(m)) && return T(NaN)
    #if !(0 ≤ m ≤ 1) throw(DomainError(m, "argument m not in [0,1]")) end
    if 2abs(phi) > T(π)
        phi2 = phi + T(π / 2)
        return 2fld(phi2, T(π)) * E(m) - _E(cos(mod(phi2, T(π))), m)
    end
    _E(sin(phi), m)
end

function _E(sinphi::A, m::B) where {A,B}
    T = promote_type(A, B)
    if m > one(T)
        minv = inv(m)
        m12 = √m
        sinalpha = m12 * sinphi
        sinalpha2 = sinalpha^2
        cosalpha2 = one(T) - sinalpha2
        yalpha = one(T) - minv * sinalpha2

        drf, ierr1 = DRF(cosalpha2, yalpha, _one(T))
        drd, ierr2 = DRD(cosalpha2, yalpha, _one(T))

        if ierr1 == ierr2 == 0
            return sinalpha*(drf -sinalpha2*drd/3)/m12
        elseif ierr1 == ierr2 == 2
            # 2 - (1+m)*sinphi2 < tol
            return sinphi
        else
            return T(NaN)
        end
    else
        sinphi2 = sinphi^2
        cosphi2 = one(T) - sinphi2
        y = one(T) - m * sinphi2
        drf, ierr1 = DRF(cosphi2, y, _one(T))
        drd, ierr2 = DRD(cosphi2, y, _one(T))

        if ierr1 == ierr2 == 0
            return sinphi * (drf - m * sinphi2 * drd / 3)
        elseif ierr1 == ierr2 == 2
            # 2 - (1+m)*sinphi2 < tol
            return sinphi
        else
            return T(NaN)
        end
    end
end

"""
`ellipke(m::Real)`
returns `(K(m), E(m))` for scalar `0 ≤ m ≤ 1`
"""
function ellipke(m::T) where {T}
    if m < one(T)
        y = 1 - m
        drf, ierr1 = DRF(_zero(T), y, _one(T))
        drd, ierr2 = DRD(_zero(T), y, _one(T))
        @assert ierr1 == 0 && ierr2 == 0
        return (drf, drf - m * drd / 3)
    elseif m == 1.0
        return (T(Inf), one(T))
    elseif isnan(m)
        return (T(NaN), T(NaN))
    else
        throw(DomainError(m, "argument m not <= 1"))
    end
end

E(m) = ellipke(m)[2]

# assumes 0 ≤ m ≤ 1
function rawF(sinphi::A, m::B) where {A,B}
    T = promote_type(A, B)
    (abs(sinphi) == one(T) && m == one(T)) && return sign(sinphi) * T(Inf)
    sinphi2 = sinphi^2
    drf, ierr = DRF(_one(T) - sinphi2, _one(T) - m * sinphi2, _one(T))
    @assert ierr == 0
    sinphi * drf
end

function _F(phi::A, m::B) where {A,B}
    T = promote_type(A, B)

    (isnan(phi) || isnan(m)) && return T(NaN)
    (m < 0 || m > 1) && throw(DomainError(m, "argument m not in [0,1]"))
    if abs(phi) > T(π / 2)
        # Abramowitz & Stegun (17.4.3)
        phi2 = phi + T(π / 2)
        return 2 * fld(phi2, T(π)) * K(m) - rawF(cos(mod(phi2, T(π))), m)
    end
    return rawF(sin(phi), m)
end

function F(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    if m > 1
        ## Abramowitz & Stegum*(17.4.15)
        m12 = sqrt(m)
        theta = asin(m12 * sin(φ))
        signθ = sign(theta)
        absθ = abs(theta)
        return signθ / m12 * _F(absθ, inv(m))
    elseif m < 0
        # Abramowitz & Stegum*(17.4.17)
        n = -m
        m12 = inv(sqrt(1 + n))
        m1m = n / (1 + n)
        newφ = T(π / 2) - φ
        signφ = sign(newφ)
        absφ = abs(newφ)
        return (m12 * K(m1m) - signφ * m12 * _F(absφ, m1m))
    end
    absφ = abs(φ)
    signφ = sign(φ)
    return signφ * _F(absφ, m)
end

function K(m::T) where {T}
    if m < 1
        drf, ierr = DRF(_zero(T), 1 - m, _one(T))
        @assert ierr == 0
        return drf
    elseif m == 1
        return T(Inf)
    elseif isnan(m)
        return T(NaN)
    else
        throw(DomainError(m, "argument m not <= 1"))
    end
end

function Pi(n::A, phi::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    (isnan(n) || isnan(phi) || isnan(m)) && return T(NaN)
    !(0 ≤ m ≤ 1) && throw(DomainError(m, "argument m not in [0,1]"))
    sinphi = sin(phi)
    sinphi2 = sinphi^2
    cosphi2 = 1 - sinphi2
    y = 1 - m * sinphi2
    drf, ierr1 = DRF(cosphi2, y, _one(T))
    drj, ierr2 = DRJ(cosphi2, y, _one(T), 1 - n * sinphi2)
    if ierr1 == 0 && ierr2 == 0
        return sinphi * (drf + n * sinphi2 * drj / 3)
    elseif ierr1 == 2 && ierr2 == 2
        # 2 - (1+m)*sinphi2 < tol
        return T(Inf)
    elseif ierr1 == 0 && ierr2 == 2
        # 1 - n*sinphi2 < tol
        return T(Inf)
    end
    return T(NaN)
end
Π = Pi

function ellipj(u::A, m::B, tol::C) where {A,B,C}

    phi = am(u, m, tol)
    s = sin(phi)
    c = cos(phi)
    d = sqrt(1 - m * s^2)
    return (s, c, d)
end

function ellipj(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    ellipj(u, m, eps(T))
end

end # module
