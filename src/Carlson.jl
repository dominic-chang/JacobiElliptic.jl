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
    drf, ierr1, drd, ierr2 = DRFD(_zero(T), y, _one(T))
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

        drf, ierr1, drd, ierr2 = DRFD(cosalpha2, yalpha, _one(T))

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
        drf, ierr1, drd, ierr2 = DRFD(cosφ2, y, _one(T))

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

# `DRFD` hands back the RF that F needs alongside the RD that E needs, so a caller wanting
# both at the same amplitude pays for one duplication rather than two. The amplitude and
# parameter reductions differ between F and E, so anything outside the principal arch (or
# any degenerate argument) falls back to the separate entry points.
@inline function _F_E(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    if isnan(φ) || isnan(m) || !(zero(T) <= m <= one(T)) || 2abs(φ) > T(π)
        return (F(φ, m), E(φ, m))
    end
    oneT = _one(T)
    sinφ = sin(φ)
    sinφ2 = sinφ * sinφ
    cosφ2 = oneT - sinφ2
    y = muladd(-m, sinφ2, oneT)
    drf, ierr1, drd, ierr2 = DRFD(cosφ2, y, oneT)
    if ierr1 == 0 && ierr2 == 0
        return (sinφ * drf, sinφ * muladd(-m * sinφ2 * T(1 / 3), drd, drf))
    end
    return (F(φ, m), E(φ, m))
end

# `DRFJ` shares one duplication between RF and RJ, and the RF it returns is exactly the
# one `_rawF` would compute for the same (sinφ, m). The reductions in `Pi` below need both
# integrals at the same argument, so this hands back the pair and they pay for one
# duplication instead of two.
@inline function _F_Pi(n::A, sinφ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    (isnan(n) || isnan(sinφ) || isnan(m)) && return (T(NaN), T(NaN))
    oneT = _one(T)
    sinφ2 = sinφ^2
    cosφ2 = oneT - sinφ2
    y = muladd(-m, sinφ2, oneT)
    p = muladd(-n, sinφ2, oneT)
    drf, ierr1, drj, ierr2 = DRFJ(cosφ2, y, oneT, p)

    f = sinφ * drf
    if ierr1 == 0 && ierr2 == 0
        return (f, sinφ * muladd(n * sinφ2 / 3, drj, drf))
    elseif ierr1 == 2 && ierr2 == 2
        # 2 - (1+m)*sinφ2 < tol; DRFJ zeroes its answers here, so F takes `_rawF`'s
        # signed pole rather than the zeroed drf.
        return (sinφ * T(Inf), T(Inf))
    elseif ierr1 == 0 && ierr2 == 2
        # 1 - n*sinφ2 < tol
        return (f, T(Inf))
    end
    return (T(NaN), T(NaN))
end

@inline _Pi(n::A, sinφ::B, m::C) where {A,B,C} = _F_Pi(n, sinφ, m)[2]

# `Pi` reduces any (n, φ, m) onto one of the three kernels below and then stops.

# Principle branch.
@inline function _Pi_principal(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    oneT = one(T)
    n <= oneT && return _Pi(n, sin(φ), m)
    # https://link.springer.com/book/10.1007/978-3-642-65138-0 117.01. n1 = m/n ∈ [0, 1),
    # so both integrals come from a single duplication.
    sinφ, cosφ = sincos(φ)
    sinφ2 = sinφ * sinφ
    t1 = sinφ / (cosφ * _sqrt(muladd(-m, sinφ2, oneT)))
    h1 = (oneT - n) * (n−m) / n
    f, pi_n = _F_Pi(m / n, sinφ, m)
    return _FukushimaT(t1, h1) - pi_n + f
end

# n ≤ 1, 0 ≤ m ≤ 1, 2|φ| > π: Π(n; φ + kπ|m) = Π(n; φ|m) + 2k Π(n|m).
@inline function _Pi_wrapped(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    φ2 = φ + T(π / 2)
    return 2fld(φ2, T(π)) * Pi(n, m) - _Pi(n, cos(mod(φ2, T(π))), m)
end

# m < 0, n ≠ m, 2|φ| ≤ π: imaginary modulus transformation,
# https://dlmf.nist.gov/19.7#iii
@inline function _Pi_negm(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    oneT = one(T)
    mc = oneT - m
    imc = inv(mc)
    mN = -m * imc
    sinφ = sin(φ)
    sinφ2 = sinφ * sinφ
    φN = asin(_sqrt(mc / muladd(-m, sinφ2, oneT)) * sinφ)
    nN = (n - m) * imc

    if nN <= oneT
        f, pi_n = _F_Pi(nN, sin(φN), mN)
        return _sqrt(imc) / nN * (mN * f + imc * n * pi_n)
    end
    # Otherwise mN and φN still land on the principal arch, with nN > 1.
    return _sqrt(imc) / nN * (mN * F(φN, mN) + imc * n * _Pi_principal(nN, φN, mN))
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
        return _Pi_principal(n / m, asin(sinθ), inv(m)) / k
    end

    if m < zero(T)
        if n == m # the transformation in `_Pi_negm` is 0/0 at n == m in this case. The closed form is used instead.
            sinφ, cosφ = sincos(φ)
            return (E(φ, m) - m * sinφ * cosφ / _sqrt(muladd(-m, sinφ * sinφ, oneT))) /
                   (oneT - m)
        end
        if 2abs(φ) > T(π)
            φ2 = φ + T(π / 2)
            return 2fld(φ2, T(π)) * Pi(n, m) + _Pi_negm(n, mod(φ2, T(π)) - T(π / 2), m)
        end
        return _Pi_negm(n, φ, m)
    end

    if n > oneT && 2abs(φ) > T(π)
        sinφ, cosφ = sincos(φ)
        sinφ2 = sinφ * sinφ
        t1 = sinφ / (cosφ * _sqrt(muladd(-m, sinφ2, oneT)))
        h1 = (oneT - n) * (n−m) / n
        return _FukushimaT(t1, h1) - _Pi_wrapped(m / n, φ, m) + F(φ, m)
    end

    2abs(φ) > T(π) && return _Pi_wrapped(n, φ, m)
    return _Pi_principal(n, φ, m)
end

# https://doi.org/T(10).1016/j.cam.2011.1107
"""
``\\Pi(n|\\,m)=\\int_{0}^{1 }{\\frac{1}{1-nt^{2}}}{\\frac{dt}{\\sqrt{\\left(1-mt^{2}\\right)\\left(1-t^{2}\\right)}}}.``

Returns the complete elliptic integral of the third kind.

# Arguments

- `n` : Characteristic
- `m` : Elliptic parameter
"""
@inline function _Pi_complete(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)
    n == oneT && return T(Inf) # 1/(1 - n*sin(θ)^2) pole at θ = π/2
    n == zero(T) && return K(m)
    m == oneT && return T(Inf)
    kc = _sqrt(oneT - m)
    nc = oneT - n
    return _cel(kc, nc, oneT, oneT)
end

function Pi(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)
    m > oneT && return T(NaN) # Complex Branch
    n > oneT && return K(m) - _Pi_complete(m / n, m)
    return _Pi_complete(n, m)
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

@inline function _J_complete(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    kc = √(1 - m)
    nc = 1 - n
    return _cel(kc, nc, zero(T), one(T))
end

function J(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    if n > one(T)
        n1 = m / n
        n1 > one(T) && return T(NaN) # The reduced characteristic does not terminate m > n,
        return n1 * _J_complete(n1, m)
    end
    return _J_complete(n, m)
end

function J(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n, φ, m = T(n), T(φ), T(m)
    oneT = one(T)
    zeroT = zero(T)
    halfπ = T(π / 2)

    acc = zeroT
    scale = oneT
    for _ = 1:6
        # Reduction of Amplitude
        φ == zeroT && return acc
        φ == halfπ && return acc + scale * J(n, m)
        if φ < zeroT
            scale = -scale
            φ = -φ
            continue
        end

        if φ > halfπ && m < oneT
            j = floor(φ / T(π))
            newφ = φ - j * T(π)
            signφ = sign(newφ)
            if abs(newφ) > halfπ
                j += signφ * oneT
                newφ = newφ - signφ * T(π)
            end
            signφ = sign(newφ)

            acc = acc + scale * (2 * j * J(n, m))
            scale = scale * signφ
            φ = abs(newφ)
            continue
        end

        # Reduction of parameter
        if zeroT < φ < halfπ
            sinφ, cosφ = sincos(φ)
            sinφ2 = sinφ^2
            m_sinφ2 = m * sinφ2
            if m_sinφ2 ≤ oneT
                nc = oneT - n
                iszero(m) &&
                    !iszero(n) &&
                    return acc + scale * ((_FukushimaT(sinφ / cosφ, nc) - φ) / n)

                iszero(m) && iszero(n) && return acc + scale * (φ / 2 - sin(2 * φ) / 4)

                isone(m) &&
                    !isone(n) &&
                    return acc + scale * ((atanh(sinφ) - _FukushimaT(sinφ, -n)) / nc)

                isone(m) &&
                    isone(n) &&
                    return acc + scale * ((sinφ / (cosφ * cosφ) - atanh(sinφ)) / 2)
            end
        end

        if oneT < m < inv(sin(φ)^2)
            mR = inv(m)
            scale = scale * (mR * √mR)
            n, φ, m = n / m, asin(√m * sin(φ)), mR
            continue
        elseif m < zeroT
            mc = oneT - m
            sinφ = sin(φ)
            sinφ2 = sinφ * sinφ
            mN = -m / mc
            scale = scale * (mN * √mN)
            n, φ, m = (n−m) / mc, asin(√(mc / muladd(-m, sinφ2, oneT)) * sinφ), mN
            continue
        end
        break
    end

    if zeroT < φ < oneT && zeroT < m < oneT
        sinφ, cosφ = sincos(φ)
        sinφ2 = sinφ * sinφ
        sqrt_term = √(muladd(-m, sinφ2, oneT))
        if n > oneT
            nc = oneT - n
            t1 = sinφ / (cosφ * sqrt_term)
            h1 = nc * (n−m) / n
            n1 = m / n
            return acc +
                   scale * ((-F(φ, m) + _FukushimaT(t1, h1) - n1 * _J(n1, sinφ, m)) / n)
        elseif n < zeroT
            mc = oneT - m
            nc = oneT - n
            t2 = sinφ * cosφ / sqrt_term
            h2 = -n * (m - n) / nc
            n2 = (m−n) / nc
            return acc +
                   scale *
                   ((F(φ, m) - _FukushimaT(t2, h2) - (mc / nc) * _J(n2, sinφ, m)) / nc)
        end
    end

    zeroT < φ < halfπ &&
        zeroT < m < oneT &&
        zeroT < n < oneT &&
        return acc + scale * _J(n, sin(φ), m)
    return T(NaN)
end

end # module
