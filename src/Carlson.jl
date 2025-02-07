module CarlsonAlg
# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi

# matlab compatible
export ellipj, ellipke

include("jacobi.jl")
include("slatec.jl")

_zero(T) = zero(T)
_one(T) = one(T)

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
    #(m < 0 || m > 1) && throw(DomainError(m, "argument m not in [0,1]"))
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
        throw(DomainError("argument m not <= 1"))
    end
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
            return sinalpha * (drf - sinalpha2 * drd / 3) / m12
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
    elseif m == 1
        return (T(Inf), one(T))
    elseif isnan(m)
        return (T(NaN), T(NaN))
    else
        throw(DomainError("argument m not <= 1"))
    end
end

E(m) = ellipke(m)[2]



function _Pi(n::A, sinphi::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    (isnan(n) || isnan(sinphi) || isnan(m)) && return T(NaN)
    #!(0 ≤ m ≤ 1) && throw(DomainError(m, "argument m not in [0,1]"))
    #sinphi = sin(phi)
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
function custom_atanh(a::T) where {T}

    arg1 = abs(one(T) + a)
    arg2 = abs(one(T) - a)

    ans = (log(arg1 / arg2)) / 2
    return ans
end

function FukushimaT(t::A, h::B) where {A,B}
    T = promote_type(A, B)
    if h > zero(T)
        return atan(t * √h) / √(h)
    elseif h == zero(T)
        return t
    else
        arg = t * √(-h)
        ans = abs(arg) < one(T) ? atanh(arg) : custom_atanh(arg)
        return ans / √(-h)
    end
end

function Pi(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)

    if m < zero(T) # Imaginary modulus transformation https://dlmf.nist.gov/19.7#iii
        mc = one(T) - m
        imc = inv(mc)
        mN = -m * imc
        φN = asin(sqrt(mc / (one(T) − m * sin(φ)^2)) * sin(φ))

        nN = (n - m) * imc

        return sqrt(imc) / nN * (mN * F(φN, mN) + imc * n * Pi(nN, φN, mN))
    end # https://link.springer.com/book/10.1007/978-3-642-65138-0 117.01
    if n > one(T)
        nc = one(T) - n
        t1 = tan(φ) / sqrt(one(T) − m * sin(φ)^2)
        h1 = nc * (n − m) / n
        n1 = m / n
        return (FukushimaT(t1, h1) - Pi(n1, φ, m) + F(φ, m))
    end
    if 2abs(φ) > T(π)
        #j = floor(φ/T(π))
        #phi2= φ - j*T(π)
        #signφ = sign(phi2)
        #if abs(phi2) > T(π/2)
        #    j += signφ*one(T)
        #    phi2= phi2 - signφ*T(π)
        #end
        #signφ = sign(phi2)
        phi2 = φ + T(π / 2)
        return 2fld(phi2, T(π)) * Pi(n, m) - _Pi(n, cos(mod(phi2, T(π))), m)
        #return 2j * Pi(n, m) - signφ*Pi(n, abs(phi2), m)
    end

    return _Pi(n, sin(φ), m)
end

#https://doi.org/T(10).1016/j.cam.2011.1107
function Pi(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    n > one(T) && return K(m) - Pi(m / n, m)
    n == zero(T) && return K(m)
    m == zero(T) || m == one(T) && return T(Inf) #atanh(√(-1 + n)*tan(θ))/√(-1 + n)
    kc = √(one(T) - m)
    nc = one(T) - n
    return cel(kc, nc, one(T), one(T))
end

Π = Pi

#https://link-springer-com.ezp-prod1.hul.harvard.edu/article/T(10).1007/BF02165405
function cel(kc::A, p::B, a::C, b::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    #ca = T(1e-6)
    ca = eps(T)
    kc = abs(kc)
    e = kc
    m = one(T)

    f, g, q = T(0), T(0), T(0)
    if p > T(0)
        p = √p
        b = b / p
    else
        f = kc^2
        q = one(T) - f
        g = one(T) - p
        f = f - p
        q = (b - a * p) * q
        p = √(f / g)
        a = (a - b) / g
        b = -q * (g^2 * p) + a * p
    end
    while true
        f = a
        a = b / p + a
        g = e / p
        b = f * g + b
        b = b + b
        p = g + p
        g = m
        m = kc + m
        if abs(g - kc) < g * ca
            break
        end
        kc = √e
        kc = kc + kc
        e = kc * m
    end
    return T(π / 2) * (a * m + b) / (m * (m + p))
end



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
