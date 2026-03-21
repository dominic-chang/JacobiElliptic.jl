module CarlsonAlg
# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi

# matlab compatible
export ellipj, ellipke

@inline function _sqrt(x::Union{Float32,Float64})
    return Core.Intrinsics.sqrt_llvm(x)
end

include("jacobi.jl")
include("slatec.jl")

_zero(T) = zero(T)
_one(T) = one(T)
function _isequals(A, B)
    return A == B;
end

# assumes 0 ≤ m ≤ 1
@inline function rawF(sinphi::A, m::B) where {A,B}
    T = promote_type(A, B)
    _isequals(abs(sinphi), 1) && _isequals(m, 1) && return sign(sinphi) * T(Inf)
    oneT = _one(T)
    sinphi2 = sinphi * sinphi
    drf, ierr = DRF(oneT - sinphi2, muladd(-m, sinphi2, oneT), oneT)
    @assert ierr == 0
    sinphi * drf
end

@inline function _F(phi::A, m::B) where {A,B}
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
    oneT = one(T)
    if m > 1
        ## Abramowitz & Stegum*(17.4.15)
        m12 = _sqrt(m)
        theta = asin(m12 * sin(φ))
        signθ = sign(theta)
        absθ = abs(theta)
        return signθ / m12 * _F(absθ, inv(m))
    elseif m < 0
        # Abramowitz & Stegum*(17.4.17)
        n = -m
        one_plus_n = oneT + n
        m12 = inv(_sqrt(one_plus_n))
        m1m = n / one_plus_n
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
    oneT = one(T)
    if isnan(m)
        return T(NaN)
    elseif m < zero(T)
        # Transformation for m < 0
        m_abs = -m
        one_plus_m_abs = m_abs + oneT
        m_transformed = m_abs / one_plus_m_abs
        k_trans = K(m_transformed)
        sqrt_factor = _sqrt(one_plus_m_abs)
        k_result = k_trans / sqrt_factor
        return k_result
    elseif m > oneT
        # Reciprocal modulus transformation for m > 1
        k = _sqrt(m)
        k_inv = inv(k)
        m_inv = inv(m)
        k_complete = K(m_inv)
        k_transformed = k_inv * k_complete
        return k_transformed
    elseif m == oneT
        return T(Inf)
    else
        # 0 ≤ m < 1
        drf, ierr = DRF(_zero(T), oneT - m, _one(T))
        @assert ierr == 0
        return drf
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

@inline function _E(sinphi::A, m::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)
    inv3 = T(1 / 3)
    if m > one(T)
        minv = inv(m)
        m12 = _sqrt(m)
        sinalpha = m12 * sinphi
        sinalpha2 = sinalpha * sinalpha
        cosalpha2 = oneT - sinalpha2
        yalpha = muladd(-minv, sinalpha2, oneT)

        drf, ierr1 = DRF(cosalpha2, yalpha, _one(T))
        drd, ierr2 = DRD(cosalpha2, yalpha, _one(T))

        if ierr1 == ierr2 == 0
            return sinalpha * muladd(-sinalpha2 * inv3, drd, drf) / m12
        elseif ierr1 == ierr2 == 2
            # 2 - (1+m)*sinphi2 < tol
            return sinphi
        else
            return T(NaN)
        end
    else
        sinphi2 = sinphi * sinphi
        cosphi2 = oneT - sinphi2
        y = muladd(-m, sinphi2, oneT)
        drf, ierr1 = DRF(cosphi2, y, _one(T))
        drd, ierr2 = DRD(cosphi2, y, _one(T))

        if ierr1 == ierr2 == 0
            return sinphi * muladd(-m * sinphi2 * inv3, drd, drf)
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
    oneT = one(T)
    if isnan(m)
        return (T(NaN), T(NaN))
    elseif m < zero(T)
        # Transformation for m < 0
        m_abs = -m
        one_plus_m_abs = m_abs + oneT
        m_transformed = m_abs / one_plus_m_abs
        k_trans, e_trans = ellipke(m_transformed)
        sqrt_factor = _sqrt(one_plus_m_abs)
        k_result = k_trans / sqrt_factor
        e_result = e_trans * sqrt_factor
        return (k_result, e_result)
    elseif m > oneT
        # Reciprocal modulus transformation for m > 1
        k = _sqrt(m)
        k_inv = inv(k)
        m_inv = inv(m)
        k_complete, e_complete = ellipke(m_inv)
        k_transformed = k_inv * k_complete
        e_transformed = k * (e_complete - (1 - m_inv) * k_complete)
        return (k_transformed, e_transformed)
    elseif m == oneT
        return (T(Inf), oneT)
    else
        # 0 ≤ m < 1
        inv3 = T(1 / 3)
        y = oneT - m
        drf, ierr1 = DRF(_zero(T), y, _one(T))
        drd, ierr2 = DRD(_zero(T), y, _one(T))
        @assert ierr1 == 0 && ierr2 == 0
        return (drf, muladd(-m * inv3, drd, drf))
    end
end

E(m) = ellipke(m)[2]

@inline function _Pi(n::A, sinphi::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    (isnan(n) || isnan(sinphi) || isnan(m)) && return T(NaN)
    #!(0 ≤ m ≤ 1) && throw(DomainError(m, "argument m not in [0,1]"))
    #sinphi = sin(phi)
    oneT = _one(T)
    sinphi2 = sinphi * sinphi
    cosphi2 = oneT - sinphi2
    y = muladd(-m, sinphi2, oneT)
    p = muladd(-n, sinphi2, oneT)
    drf, ierr1, drj, ierr2 = DRFJ(cosphi2, y, oneT, p)

    if ierr1 == 0 && ierr2 == 0
        return sinphi * muladd(n * sinphi2 / 3, drj, drf)
    elseif ierr1 == 2 && ierr2 == 2
        # 2 - (1+m)*sinphi2 < tol
        return T(Inf)
    elseif ierr1 == 0 && ierr2 == 2
        # 1 - n*sinphi2 < tol
        return T(Inf)
    end
    return T(NaN)
end

function Pi(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    oneT = one(T)

    if m < zero(T) # Imaginary modulus transformation https://dlmf.nist.gov/19.7#iii
        mc = oneT - m
        imc = inv(mc)
        mN = -m * imc
        sinφ = sin(φ)
        sinφ2 = sinφ * sinφ
        φN = asin(_sqrt(mc / muladd(-m, sinφ2, oneT)) * sinφ)

        nN = (n - m) * imc

        return _sqrt(imc) / nN * (mN * F(φN, mN) + imc * n * Pi(nN, φN, mN))
    end # https://link.springer.com/book/10.1007/978-3-642-65138-0 117.01
    if n > one(T)
        nc = oneT - n
        sinφ, cosφ = sincos(φ)
        sinφ2 = sinφ * sinφ
        t1 = sinφ / (cosφ * _sqrt(muladd(-m, sinφ2, oneT)))
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
    oneT = one(T)
    n > oneT && return K(m) - Pi(m / n, m)
    n == zero(T) && return K(m)
    m == zero(T) || m == oneT && return T(Inf) #atanh(_sqrt(-1 + n)*tan(θ))/_sqrt(-1 + n)
    kc = _sqrt(oneT - m)
    nc = oneT - n
    return cel(kc, nc, oneT, oneT)
end

Π = Pi

function ellipj(u, m)
    phi = am(u, m)
    s, c = sincos(phi)
    d = _sqrt(muladd(-m, s * s, 1))
    return (s, c, d)
end

end # module
