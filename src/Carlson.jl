module CarlsonAlg
# elliptic integrals of 1st/2nd/3rd kind
export E, F, K, Pi

# matlab compatible
export ellipj, ellipke

@inline function _sqrt(x::Union{Float32,Float64})
    return Core.Intrinsics.sqrt_llvm(x)
end
_zero(T) = zero(T)
_one(T) = one(T)
function _isequals(A, B)
    return A == B
end

include("jacobi.jl")
include("slatec.jl")
include("complete_reduction.jl")

# Kernels for the shared complete-integral reduction, valid on `0 ≤ m < 1`.
@inline function _ellipke_base(m::T) where {T}
    inv3 = T(1 / 3)
    y = one(T) - m
    drf, ierr1 = DRF(_zero(T), y, _one(T))
    drd, ierr2 = DRD(_zero(T), y, _one(T))
    @assert ierr1 == 0 && ierr2 == 0
    return (drf, muladd(-m * inv3, drd, drf))
end

# `K` keeps only the first component, so this skips the `DRD` call and hands back a
# placeholder; the compiler drops it along with the unused half of the lift.
@inline function _K_base(m::T) where {T}
    drf, ierr = DRF(_zero(T), one(T) - m, _one(T))
    @assert ierr == 0
    return (drf, zero(T))
end

# assumes 0 ≤ m ≤ 1
@inline function _rawF(sinφ::A, m::B) where {A,B}
    T = promote_type(A, B)
    _isequals(abs(sinφ), 1) && _isequals(m, 1) && return sign(sinφ) * T(Inf)
    oneT = _one(T)
    sinφ2 = sinφ * sinφ
    drf, ierr = DRF(oneT - sinφ2, muladd(-m, sinφ2, oneT), oneT)
    @assert ierr == 0
    sinφ * drf
end

@inline function _F(φ::A, m::B) where {A,B}
    T = promote_type(A, B)

    (isnan(φ) || isnan(m)) && return T(NaN)
    # (m < 0 || m > 1) && throw(DomainError(m, "argument m not in [0,1]"))
    if abs(φ) > T(π / 2)
        # Abramowitz & Stegun (17.4.3)
        φ2 = φ + T(π / 2)
        return 2 * fld(φ2, T(π)) * K(m) - _rawF(cos(mod(φ2, T(π))), m)
    end
    return _rawF(sin(φ), m)
end

"""
``F(\\varphi |\\, m) = \\int_0^{\\varphi}\\frac{d\\theta}{\\sqrt{1-m\\sin(\\theta)^2}}.``

Returns the incomplete elliptic integral of the first kind.

# Arguments

- `φ` : Amplitude
- `m` : Elliptic parameter 
"""
function F(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)
    if m > 1
        ## Abramowitz & Stegum*(17.4.15)
        m12 = _sqrt(m)
        sinθ = m12 * sin(φ)
        # F(φ|m) is real only on the first arch |φ| ≤ asin(1/√m); past it 1 - m*sin(θ)^2
        # turns negative somewhere in [0, φ]. |sinθ| > 1 catches most of that, but once
        # |φ| > π/2 the sine folds back into range and the transformation below would
        # silently evaluate the wrong branch instead.
        (abs(φ) > T(π / 2) || abs(sinθ) > oneT) && return T(NaN)
        theta = asin(sinθ)
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

K(m::T) where {T} = _ellipke_reduced(_K_base, m)[1]

"""
``E(\\varphi|\\, m) = \\int_0^{\\varphi}d\\theta \\sqrt{1-m\\sin(\\theta)^2}.``

Returns the incomplete elliptic integral of the second kind.

# Arguments

- `φ` : Amplitude
- `m` : Elliptic parameter

"""
function E(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    (isnan(φ) || isnan(m)) && return T(NaN)
    # if !(0 ≤ m ≤ 1) throw(DomainError(m, "argument m not in [0,1]")) end
    if 2abs(φ) > T(π)
        φ2 = φ + T(π / 2)
        return 2fld(φ2, T(π)) * E(m) - _E(cos(mod(φ2, T(π))), m)
    end
    _E(sin(φ), m)
end

@inline function _E(sinφ::A, m::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)
    inv3 = T(1 / 3)
    if m > one(T)
        minv = inv(m)
        m12 = _sqrt(m)
        sinalpha = m12 * sinφ
        sinalpha2 = sinalpha * sinalpha
        cosalpha2 = oneT - sinalpha2
        yalpha = muladd(-minv, sinalpha2, oneT)

        drf, ierr1 = DRF(cosalpha2, yalpha, _one(T))
        drd, ierr2 = DRD(cosalpha2, yalpha, _one(T))

        if ierr1 == ierr2 == 0
            return sinalpha * muladd(-sinalpha2 * inv3, drd, drf) / m12
        elseif ierr1 == ierr2 == 2
            # 2 - (1+m)*sinφ2 < tol
            return sinφ
        else
            return T(NaN)
        end
    else
        sinφ2 = sinφ * sinφ
        cosφ2 = oneT - sinφ2
        y = muladd(-m, sinφ2, oneT)
        drf, ierr1 = DRF(cosφ2, y, _one(T))
        drd, ierr2 = DRD(cosφ2, y, _one(T))

        if ierr1 == ierr2 == 0
            return sinφ * muladd(-m * sinφ2 * inv3, drd, drf)
        elseif ierr1 == ierr2 == 2
            # 2 - (1+m)*sinφ2 < tol
            return sinφ
        else
            return T(NaN)
        end
    end
end

"""
`ellipke(m::Real)`
returns `(K(m), E(m))` for scalar `0 ≤ m ≤ 1`
"""
ellipke(m::T) where {T} = _ellipke_reduced(_ellipke_base, m)

E(m) = ellipke(m)[2]

@inline function _Pi(n::A, sinφ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    (isnan(n) || isnan(sinφ) || isnan(m)) && return T(NaN)
    #!(0 ≤ m ≤ 1) && throw(DomainError(m, "argument m not in [0,1]"))
    # sinφ = sin(φ)
    oneT = _one(T)
    sinφ2 = sinφ ^ 2
    cosφ2 = oneT - sinφ2
    y = muladd(-m, sinφ2, oneT)
    p = muladd(-n, sinφ2, oneT)
    drf, ierr1, drj, ierr2 = DRFJ(cosφ2, y, oneT, p)

    if ierr1 == 0 && ierr2 == 0
        return sinφ * muladd(n * sinφ2 / 3, drj, drf)
    elseif ierr1 == 2 && ierr2 == 2
        # 2 - (1+m)*sinφ2 < tol
        return T(Inf)
    elseif ierr1 == 0 && ierr2 == 2
        # 1 - n*sinφ2 < tol
        return T(Inf)
    end
    return T(NaN)
end

"""
``\\Pi (n;\\varphi \\,|\\,m)=\\int _{0}^{\\sin \\varphi }{\\frac {1}{1-nt^{2}}}{\\frac {dt}{\\sqrt {\\left(1-mt^{2}\\right)\\left(1-t^{2}\\right)}}}.``

Returns the incomplete elliptic integral of the third kind.

# Arguments

- `n` : Characteristic
- `φ` : Amplitude
- `m` : Elliptic parameter
"""
function Pi(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    oneT = one(T)

    if m > oneT
        # Reciprocal modulus transformation: Π(n;φ|m) = Π(n/m; θ|1/m)/√m with sinθ = √m sinφ.
        # Applied before the n > 1 reduction below, whose n₁ = m/n is an exact 2-cycle when
        # n and m both exceed 1 (n₂ = m/n₁ = n) and never terminates. Same real domain as F.
        k = _sqrt(m)
        sinθ = k * sin(φ)
        (abs(φ) > T(π / 2) || abs(sinθ) > oneT) && return T(NaN)
        return Pi(n / m, asin(sinθ), inv(m)) / k
    end

    if m < zero(T) # Imaginary modulus transformation https://dlmf.nist.gov/19.7# iii
        if n == m # the transformation below is 0/0 at n == m, so use the closed form
            sinφ, cosφ = sincos(φ)
            return (E(φ, m) - m * sinφ * cosφ / _sqrt(muladd(-m, sinφ * sinφ, oneT))) /
                   (oneT - m)
        end
        if 2abs(φ) > T(π)
            φ2 = φ + T(π / 2)
            return 2fld(φ2, T(π)) * Pi(n, m) + Pi(n, mod(φ2, T(π)) - T(π / 2), m)
        end
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
        h1 = nc * (n−m) / n
        n1 = m / n
        return (_FukushimaT(t1, h1) - Pi(n1, φ, m) + F(φ, m))
    end
    if 2abs(φ) > T(π)
        # j = floor(φ/T(π))
        # φ2= φ - j*T(π)
        # signφ = sign(φ2)
        # if abs(φ2) > T(π/2)
        #    j += signφ*one(T)
        #    φ2= φ2 - signφ*T(π)
        # end
        # signφ = sign(φ2)
        φ2 = φ + T(π / 2)
        return 2fld(φ2, T(π)) * Pi(n, m) - _Pi(n, cos(mod(φ2, T(π))), m)
        # return 2j * Pi(n, m) - signφ*Pi(n, abs(φ2), m)
    end

    return _Pi(n, sin(φ), m)
end

# https://doi.org/T(10).1016/j.cam.2011.1107
"""
``\\Pi(n|\\,m)=\\int_{0}^{1 }{\\frac{1}{1-nt^{2}}}{\\frac{dt}{\\sqrt{\\left(1-mt^{2}\\right)\\left(1-t^{2}\\right)}}}.``

Returns the complete elliptic integral of the third kind.

# Arguments

- `n` : Characteristic
- `m` : Elliptic parameter
"""
function Pi(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)
    zeroT = zero(T)
    m > oneT && return T(NaN) # Complex Branch
    n > oneT && return K(m) - Pi(m / n, m) # Cauchy principal value
    n == oneT && return T(Inf) # 1/(1 - n*sin(θ)^2) pole at θ = π/2
    n == zeroT && return K(m)
    m == oneT && return T(Inf)
    kc = _sqrt(oneT - m)
    nc = oneT - n
    return _cel(kc, nc, oneT, oneT)
end

Π = Pi

function ellipj(u, m)
    φ = am(u, m)
    s, c = sincos(φ)
    d = _sqrt(muladd(-m, s * s, 1))
    return (s, c, d)
end

function _J(n, sinφ, m)
    oneT = one(typeof(sinφ))
    sinφ2 = sinφ ^ 2
    cosφ2 = oneT - sinφ2
    y = muladd(-m, sinφ2, oneT)
    p = muladd(-n, sinφ2, oneT)
    return sinφ * (sinφ2 / 3) * DRJ(cosφ2, y, oneT, p)[1]
end

function J(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    n > one(T) && return m / n * J(m / n, m)
    kc = √(1 - m)
    nc = 1 - n
    return _cel(kc, nc, zero(T), one(T))
end

function J(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    oneT = one(T)
    # Reduction of Amplitude
    φ == zero(T) && return zero(T)
    φ == T(π / 2) && return J(n, m)
    φ < zero(T) && return -J(n, -φ, m)

    if abs(φ) > T(π / 2) && m < oneT
        j = floor(φ / T(π))
        newφ = φ - j * T(π)
        signφ = sign(newφ)
        if abs(newφ) > T(π / 2)
            j += signφ * oneT
            newφ = newφ - signφ * T(π)
        end
        signφ = sign(newφ)

        return 2 * j * J(n, m) + signφ * J(n, abs(newφ), m)
    end

    # Reduction of parameter
    if zero(T) < φ < T(π / 2)
        sinφ, cosφ = sincos(φ)
        sinφ2 = sinφ^2
        m_sinφ2 = m * sinφ2
        if m_sinφ2 ≤ oneT
            nc = oneT - n
            iszero(m) && !iszero(n) && return (_FukushimaT(sinφ / cosφ, nc) - φ) / n

            iszero(m) && iszero(n) && return φ / 2 - sin(2 * φ) / 4

            isone(m) && !isone(n) && return (atanh(sinφ) - _FukushimaT(sinφ, -n)) / nc

            isone(m) && isone(n) && return (sinφ / (cosφ * cosφ) - atanh(sinφ)) / 2
        end
    end

    if oneT < m < inv(sin(φ)^2)
        φR = asin(√m * sin(φ))
        nR = n / m
        mR = inv(m)
        # return NaN
        return mR * √mR * J(nR, φR, mR)
    elseif m < zero(T)
        mc = oneT - m
        sinφ = sin(φ)
        sinφ2 = sinφ * sinφ
        φN = asin(√(mc / muladd(-m, sinφ2, oneT)) * sinφ)
        nN = (n−m) / mc
        mN = -m / mc
        return mN * √mN * J(nN, φN, mN)
    end
    # Reduction of Characteristics
    if zero(T) < φ < oneT && zero(T) < m < oneT
        sinφ, cosφ = sincos(φ)
        sinφ2 = sinφ * sinφ
        sqrt_term = √(muladd(-m, sinφ2, oneT))
        if n > one(T)
            # t = inv(x^2)
            # h = y - x
            nc = oneT - n
            t1 = sinφ / (cosφ * sqrt_term)
            h1 = nc * (n−m) / n
            n1 = m / n
            return (-F(φ, m) + _FukushimaT(t1, h1) - n1 * J(n1, φ, m)) / n
        elseif n < zero(T)
            mc = oneT - m
            nc = oneT - n
            t2 = sinφ * cosφ / sqrt_term
            h2 = -n * (m - n) / nc
            n2 = (m−n) / nc
            # return NaN
            return (F(φ, m) - _FukushimaT(t2, h2) - (mc / nc) * J(n2, φ, m)) / nc
        end
    end

    zero(T) < φ < T(π / 2) &&
        zero(T) < m < oneT &&
        zero(T) < n < oneT &&
        return _J(n, sin(φ), m)
    return T(NaN)
end

end # module
