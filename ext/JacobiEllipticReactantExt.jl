module JacobiEllipticReactantExt

using JacobiElliptic
using Reactant
import JacobiElliptic: CarlsonAlg, ArithmeticGeometricMeanAlg, StaticArrays

function _am_buffer(::T) where T
    zeroT = zero(T)
    #return StaticArrays.@MVector
    return [
        zeroT,
        zeroT,
        zeroT,
        zeroT,
        zeroT,
        zeroT,
        zeroT,
        zeroT,
        zeroT,
        zeroT,
        zeroT,
    ]
end

function _reactant_am_step(
    a::A,
    b::B,
    c::C,
    n::D,
    tol::E,
    _ambuf::G,
    idx::H,
) where {A,B,C,D,E,G,H}
    T = promote_type(A, B, C, E)
    zeroT = zero(T)
    next_a = (a + b) / 2
    next_b = sqrt(a * b)
    next_c = (a - b) / 2
    active = abs(c) > tol

    a = Base.ifelse(active, next_a, a)
    b = Base.ifelse(active, next_b, b)
    c = Base.ifelse(active, next_c, c)
    n = Base.ifelse(active, n + 1, n)
    _ambuf[idx] = Base.ifelse(active, next_c / next_a, zeroT)

    return a, b, c, n, _ambuf
end

function __reactant_am(u::A, m::B, tol::C) where {A,B,C}

    T = promote_type(A, B, C)
    zeroT = zero(T)

    # Use an immutable local buffer so GPU backends can keep this in registers.
    _ambuf = _am_buffer(zeroT)
    ans = zeroT

    sqrt_tol = sqrt(tol)
    m1 = one(T) - m
    t = tanh(u)
    flag = true 
    ans, flag = Base.ifelse(iszero(u), 
        (zeroT, false), 
        Base.ifelse(m < sqrt_tol, 
            (u - m * (u - sin(2 * u) / 2) / 4, false),
            Base.ifelse(m1 < sqrt_tol,
            (asin(t) + m1 * (t - u * (one(T) - t^2)) * cosh(u) / 4, false),
            (zeroT, true)
            )
        )
    )

    NaNT = T(NaN)
    @trace if flag
        a = one(T)
        b = sqrt(m1)
        c = sqrt(m)
        n = 0

        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 1)
        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 2)
        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 3)
        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 4)
        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 5)
        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 6)
        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 7)
        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 8)
        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 9)
        a, b, c, n, _ambuf = _reactant_am_step(a, b, c, n, tol, _ambuf, 10)
    
        ans = Base.ifelse(abs(c) > tol, NaNT, ans)

        #phi = ldexp(a*u, n) # Was slower on my benchmarks
        phi = a * u * (2^n)
        for i in 10:-1:1
            next_phi = (phi + asin(_ambuf[i] * sin(phi))) / 2
            phi = Base.ifelse(n >= i, next_phi, phi)
        end
        ans = phi

    end
    return ans
end

function __reactant_am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    return __reactant_am(u, m, eps(T))
end

"""
    am(u::Real, m::Real)

Returns amplitude, φ, such that u = F(φ | m)
"""
function _reactant_am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    ans = zero(T)
    piT = T(Base.π)
    halfpiT = T(Base.π / 2)
    @trace if m < 0
        mu1 = inv(1 - m)
        mu = -m * mu1
        sqrtmu1 = sqrt(mu1)
        v = u / sqrtmu1
        phi = __reactant_am(v, mu)
        s = sin(phi)
        t = floor((phi + halfpiT) / piT)

        ans = t * piT + cospi(t) * asin(sqrtmu1 * s / sqrt(1 - mu * s^2))
    elseif m <= 1 # 0 <= m <= 1
        ans = __reactant_am(u, m)
    else # m > 1
        k = sqrt(m)
        ans = asin(inv(k) * sin(__reactant_am(k * u, inv(m))))
    end
    return ans
end

CarlsonAlg.am(u::Reactant.TracedRNumber, m::Real) = _reactant_am(u, m)
CarlsonAlg.am(u::Real, m::Reactant.TracedRNumber) = _reactant_am(u, m)
CarlsonAlg.am(u::Reactant.TracedRNumber, m::Reactant.TracedRNumber) = _reactant_am(u, m)

for (f, a, b, c) in (
    (:sn, :(sin(phi)), :(sqrtmu1 * s), :(sqrt(mu) * sin(phi))),
    (:cn, :(cos(phi)), :(cos(phi)), :(sqrt(1 - mu * sin(phi)^2))),
    (:dn, :(sqrt(1 - m * sin(phi)^2)), :(one(T)), :(cos(phi))),
)
    reactant_f = Symbol(:_reactant_, f)
    @eval begin
        function $reactant_f(u::A, m::B) where {A,B}
            T = promote_type(A, B)
            ans = T(NaN)

            @trace if m < 0
                mu1 = inv(1 - m)
                mu = -m * mu1
                sqrtmu1 = sqrt(mu1)
                v = u / sqrtmu1
                phi = __reactant_am(v, mu)
                s = sin(phi)
                ans = ($b) / sqrt(1 - mu * s^2)
            elseif m > 1
                mu = inv(m)
                v = u * sqrt(m)
                phi = __reactant_am(v, mu)
                ans = $c
            else
                phi = __reactant_am(u, m)
                ans = $a
            end

            return ans
        end

        CarlsonAlg.$f(u::Reactant.TracedRNumber, m::Real) = $reactant_f(u, m)
        CarlsonAlg.$f(u::Real, m::Reactant.TracedRNumber) = $reactant_f(u, m)
        CarlsonAlg.$f(u::Reactant.TracedRNumber, m::Reactant.TracedRNumber) = $reactant_f(u, m)
    end
end

xn = ((:s, :(sn(u, m))), (:c, :(cn(u, m))), (:d, :(dn(u, m))), (:n, 1))
for (p, num) in xn, (q, den) in xn
    f = Symbol(p, q)

    if p == q
        @eval begin
            CarlsonAlg.$f(::Reactant.TracedRNumber{T}, ::S) where {T,S<:Real} = one(promote_type(T, S))
            CarlsonAlg.$f(::S, ::Reactant.TracedRNumber{T}) where {S<:Real,T} = one(promote_type(S, T))
            CarlsonAlg.$f(::Reactant.TracedRNumber{T1}, ::Reactant.TracedRNumber{T2}) where {T1,T2} =
                one(promote_type(T1, T2))
        end
    elseif q != :n
        @eval begin
            CarlsonAlg.$f(u::Reactant.TracedRNumber, m::Real) = ($num) / ($den)
            CarlsonAlg.$f(u::Real, m::Reactant.TracedRNumber) = ($num) / ($den)
            CarlsonAlg.$f(u::Reactant.TracedRNumber, m::Reactant.TracedRNumber) = ($num) / ($den)
        end
    end
end

function _reactant_DRD_ifbody(X::A, Y::B, Z::C) where {A,B,C}
    T = promote_type(A, B, C)
    oneT = one(T)
    threeT = T(3)
    sixT = T(6)
    inv4 = T(1 / 4)
    inv5 = T(1 / 5)
    ERRTOL = T((eps(Reactant.unwrapped_eltype(T)) / 6)^Reactant.unwrapped_eltype(T)(1 / 6))
    C1 = T(3 / 14)
    C2 = T(1 / 6)
    C3 = T(9 / 22)
    C4 = T(3 / 26)

    XN = X
    YN = Y
    ZN = Z
    SIGMA = zero(T)
    POWER4 = oneT
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)

    for _ in 1:12
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        YNROOTZNROOT = YNROOT * ZNROOT
        LAMDA = muladd(XNROOT, YNROOT + ZNROOT, YNROOTZNROOT)
        SIGMA = muladd(POWER4, inv(ZNROOT * (ZN + LAMDA)), SIGMA)
        POWER4 *= inv4
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4
    end

    MU = (XN + YN + threeT * ZN) * inv5
    invMU = inv(MU)
    XNDEV = (MU - XN) * invMU
    YNDEV = (MU - YN) * invMU
    ZNDEV = (MU - ZN) * invMU

    EA = XNDEV * YNDEV
    EB = ZNDEV * ZNDEV
    EC = EA - EB
    ED = EA - sixT * EB
    EF = ED + EC + EC
    S1 = ED * (-C1 + C3 * ED * inv4 - threeT * C4 * ZNDEV * EF / 2)
    S2 = ZNDEV * (C2 * EF + ZNDEV * (-C3 * EC + ZNDEV * C4 * EA))
    return threeT * SIGMA + POWER4 * (oneT + S1 + S2) / (MU * sqrt(MU))
end

function _reactant_DRD(X::A, Y::B, Z::C) where {A,B,C}
    T = promote_type(A, B, C)
    Tc = Reactant.unwrapped_eltype(T)
    X = T(X)
    Y = T(Y)
    Z = T(Z)
    zeroT = zero(T)

    ERRTOL = T((eps(Tc) / 6)^Tc(1 / 6))
    LOLIM = T(2 / (floatmax(Tc))^Tc(2 / 3))
    TUPLIM = T((ERRTOL / 10)^Tc(1 / 3) / floatmin(Tc)^Tc(1 / 3))
    UPLIM = TUPLIM^2

    ans = zeroT
    err = Base.ifelse(
        min(X, Y) < zeroT,
        1,
        Base.ifelse(max(X, Y, Z) > UPLIM, 3, Base.ifelse(min(X + Y, Z) < LOLIM, 2, 0)),
    )

    @trace if err == 0
        ans = _reactant_DRD_ifbody(X, Y, Z)
    end

    return (ans, err)
end

CarlsonAlg.DRD(X::Reactant.TracedRNumber, Y::Real, Z::Real) =
    _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(X::Real, Y::Reactant.TracedRNumber, Z::Real) =
    _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(X::Real, Y::Real, Z::Reactant.TracedRNumber) =
    _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(
    X::Reactant.TracedRNumber,
    Y::Reactant.TracedRNumber,
    Z::Real,
) = _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(
    X::Reactant.TracedRNumber,
    Y::Real,
    Z::Reactant.TracedRNumber,
) = _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(
    X::Real,
    Y::Reactant.TracedRNumber,
    Z::Reactant.TracedRNumber,
) = _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(
    X::Reactant.TracedRNumber,
    Y::Reactant.TracedRNumber,
    Z::Reactant.TracedRNumber,
) = _reactant_DRD(X, Y, Z)

function _reactant_DRF_ifbody(X::A, Y::B, Z::C, ERRTOL::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    C1 = T(1 / 24)
    C2 = T(3 / 44)
    C3 = T(1 / 14)
    C0 = T(1 / 10)
    inv3 = T(1 / 3)
    inv4 = T(1 / 4)

    XN = X
    YN = Y
    ZN = Z
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)

    for _ in 1:10
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        YNROOTZNROOT = YNROOT * ZNROOT
        LAMDA = muladd(XNROOT, YNROOT + ZNROOT, YNROOTZNROOT)
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4
    end

    MU = (XN + YN + ZN) * inv3
    ninvMU = -inv(MU)
    XNDEV = muladd(ninvMU, MU + XN, 2)
    YNDEV = muladd(ninvMU, MU + YN, 2)
    ZNDEV = muladd(ninvMU, MU + ZN, 2)

    XNDEVYNDEV = XNDEV * YNDEV
    E2 = muladd(-ZNDEV, ZNDEV, XNDEVYNDEV)
    E3 = XNDEVYNDEV * ZNDEV
    S = one(T) + muladd(E2, muladd(-C2, E3, muladd(C1, E2, -C0)), C3 * E3)
    return S / sqrt(MU)
end

function _reactant_DRF(X::A, Y::B, Z::C) where {A,B,C}
    T = promote_type(A, B, C)
    Tc = Reactant.unwrapped_eltype(T)
    X = T(X)
    Y = T(Y)
    Z = T(Z)
    zeroT = zero(T)

    ERRTOL = T((4 * eps(Tc) / 2)^Tc(1 / 6))
    LOLIM = T(5 * floatmin(Tc))
    UPLIM = T(floatmax(Tc) / 5)

    ans = zeroT
    err = Base.ifelse(
        min(X, Y, Z) < zeroT,
        1,
        Base.ifelse(max(X, Y, Z) > UPLIM, 3, Base.ifelse(min(X + Y, X + Z, Y + Z) < LOLIM, 2, 0)),
    )

    #ans = Base.ifelse(err == 0, _reactant_valid_DRF(X, Y, Z, ERRTOL), ans)
    @trace if err == 0
        ans = _reactant_DRF_ifbody(X, Y, Z, ERRTOL)
    end

    return (ans, err)
end

CarlsonAlg.DRF(X::Reactant.TracedRNumber, Y::Real, Z::Real) =
    _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(X::Real, Y::Reactant.TracedRNumber, Z::Real) =
    _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(X::Real, Y::Real, Z::Reactant.TracedRNumber) =
    _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(
    X::Reactant.TracedRNumber,
    Y::Reactant.TracedRNumber,
    Z::Real,
) = _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(
    X::Reactant.TracedRNumber,
    Y::Real,
    Z::Reactant.TracedRNumber,
) = _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(
    X::Real,
    Y::Reactant.TracedRNumber,
    Z::Reactant.TracedRNumber,
) = _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(
    X::Reactant.TracedRNumber,
    Y::Reactant.TracedRNumber,
    Z::Reactant.TracedRNumber,
) = _reactant_DRF(X, Y, Z)

function _reactant_ellipke_base(m::T) where T
    oneT = one(T)
    twoT = T(2)
    halfpi = T(π / 2)

    a = oneT
    b = sqrt(oneT - m)
    sum_c2 = m / twoT
    pow2 = oneT

    # Arithmetic Geometric mean implementation
    # https://www.math.emory.edu/~gliang7/AGM.pdf
    # https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean
    for _ in 1:8
        c = (a - b) / twoT
        next_a = (a + b) / twoT
        next_b = sqrt(a * b)
        sum_c2 += pow2 * c * c
        pow2 += pow2
        a = next_a
        b = next_b
    end

    k = halfpi / a
    e = k * (oneT - sum_c2)
    return (k, e)
end

#----------------------------------------------------------------------------------------
# Elliptic K(m)
#----------------------------------------------------------------------------------------

function _reactant_K(m::T) where T
    oneT = one(T)
    zeroT = zero(T)
    nanT = T(NaN)
    infT = T(Inf)
    ans = zeroT

    Reactant.@trace if isnan(m)
        ans = nanT
    elseif m < zeroT
        m_abs = -m
        one_plus_m_abs = m_abs + oneT
        m_transformed = m_abs / one_plus_m_abs
        k, _ = _reactant_ellipke_base(m_transformed)
        ans = k / sqrt(one_plus_m_abs)
    elseif m > oneT
        m_inv = inv(m)
        k, _ = _reactant_ellipke_base(m_inv)
        ans = inv(sqrt(m)) * k
    elseif m == oneT
        ans = infT
    else
        k, _ = _reactant_ellipke_base(m)
        ans = k
    end
    return ans
end

ArithmeticGeometricMeanAlg.K(m::Reactant.TracedRNumber) = _reactant_K(m)

#----------------------------------------------------------------------------------------
# Elliptic E(m)
#----------------------------------------------------------------------------------------
function _reactant_E(m::T) where T
    oneT = one(T)
    zeroT = zero(T)
    nanT = T(NaN)
    ans = zeroT

    Reactant.@trace if isnan(m)
        ans = nanT
    elseif m < zeroT
        m_abs = -m
        one_plus_m_abs = m_abs + oneT
        m_transformed = m_abs / one_plus_m_abs
        _, e_trans = _reactant_ellipke_base(m_transformed)
        ans = sqrt(one_plus_m_abs) * e_trans
    elseif m > oneT
        m_inv = inv(m)
        k_complete, e_complete = _reactant_ellipke_base(m_inv)
        ans = sqrt(m) * (e_complete - (oneT - m_inv) * k_complete)
    elseif m == oneT
        ans = oneT
    else
        _, ans = _reactant_ellipke_base(m)
    end

    return ans
end

ArithmeticGeometricMeanAlg.E(m::Reactant.TracedRNumber) = _reactant_E(m)
ArithmeticGeometricMeanAlg.ellipke(m::Reactant.TracedRNumber) = (_reactant_K(m), _reactant_E(m))

#----------------------------------------------------------------------------------------
# Elliptic F(ϕ, m)
#----------------------------------------------------------------------------------------

function _reactant_rawF(sinphi::A, m::B) where {A,B}
    T = promote_type(A, B)
    sinphi = T(sinphi)
    m = T(m)
    oneT = one(T)
    infT = T(Inf)
    sinphi2 = sinphi * sinphi
    drf, _ = _reactant_DRF(oneT - sinphi2, muladd(-m, sinphi2, oneT), oneT)
    return Base.ifelse((abs(sinphi) == oneT) & (m == oneT), sign(sinphi) * infT, sinphi * drf)
end

function _reactant_internal_F(phi::A, m::B) where {A,B}
    T = promote_type(A, B)
    phi = T(phi)
    m = T(m)
    halfpi = T(π / 2)
    piT = T(π)
    nanT = T(NaN)
    ans = zero(T)

    Reactant.@trace if isnan(phi) | isnan(m)
        ans = nanT
    elseif abs(phi) > halfpi
        phi2 = phi + halfpi
        ans = 2 * fld(phi2, piT) * K(m) - _reactant_rawF(cos(mod(phi2, piT)), m)
    else
        ans = _reactant_rawF(sin(phi), m)
    end

    return ans
end

function _reactant_F(phi::A, m::B) where {A,B}
    T = promote_type(A, B)
    phi = T(phi)
    m = T(m)
    oneT = one(T)
    halfpi = T(Base.π / 2)
    ans = zero(T)

    Reactant.@trace if m > oneT
        m12 = sqrt(m)
        theta = asin(m12 * sin(phi))
        ans = sign(theta) / m12 * _reactant_internal_F(abs(theta), inv(m))
    elseif m < zero(T)
        n = -m
        one_plus_n = oneT + n
        m12 = inv(sqrt(one_plus_n))
        m1m = n / one_plus_n
        newphi = halfpi - phi
        ans = m12 * K(m1m) - sign(newphi) * m12 * _reactant_internal_F(abs(newphi), m1m)
    else
        ans = sign(phi) * _reactant_internal_F(abs(phi), m)
    end

    return ans
end

CarlsonAlg.F(φ::Reactant.TracedRNumber, m::Real) = _reactant_F(φ, m)
CarlsonAlg.F(φ::Real, m::Reactant.TracedRNumber) = _reactant_F(φ, m)
CarlsonAlg.F(φ::Reactant.TracedRNumber, m::Reactant.TracedRNumber) = _reactant_F(φ, m)

#----------------------------------------------------------------------------------------
# Incomplete Elliptic E(ϕ, m)
#----------------------------------------------------------------------------------------

function _reactant_rawE(sinphi::A, m::B) where {A,B}
    T = promote_type(A, B)
    sinphi = T(sinphi)
    m = T(m)
    oneT = one(T)
    inv3 = T(1 / 3)
    nanT = T(NaN)
    ans = zero(T)

    Reactant.@trace if m > oneT
        minv = inv(m)
        m12 = sqrt(m)
        sinalpha = m12 * sinphi
        sinalpha2 = sinalpha * sinalpha
        cosalpha2 = oneT - sinalpha2
        yalpha = muladd(-minv, sinalpha2, oneT)
        drf, ierr1 = _reactant_DRF(cosalpha2, yalpha, oneT)
        drd, ierr2 = _reactant_DRD(cosalpha2, yalpha, oneT)

        ans = Base.ifelse(
            (ierr1 == 0) & (ierr2 == 0),
            sinalpha * muladd(-sinalpha2 * inv3, drd, drf) / m12,
            Base.ifelse((ierr1 == 2) & (ierr2 == 2), sinphi, nanT),
        )
    else
        sinphi2 = sinphi * sinphi
        cosphi2 = oneT - sinphi2
        y = muladd(-m, sinphi2, oneT)
        drf, ierr1 = _reactant_DRF(cosphi2, y, oneT)
        drd, ierr2 = _reactant_DRD(cosphi2, y, oneT)

        ans = Base.ifelse(
            (ierr1 == 0) & (ierr2 == 0),
            sinphi * muladd(-m * sinphi2 * inv3, drd, drf),
            Base.ifelse((ierr1 == 2) & (ierr2 == 2), sinphi, nanT),
        )
    end

    return ans
end

function _reactant_internal_E(phi::A, m::B) where {A,B}
    T = promote_type(A, B)
    phi = T(phi)
    m = T(m)
    halfpi = T(Base.π / 2)
    piT = T(Base.π)
    nanT = T(NaN)
    ans = zero(T)

    Reactant.@trace if isnan(phi) | isnan(m)
        ans = nanT
    elseif 2 * abs(phi) > piT
        phi2 = phi + halfpi
        ans = 2 * fld(phi2, piT) * E(m) - _reactant_rawE(cos(mod(phi2, piT)), m)
    else
        ans = _reactant_rawE(sin(phi), m)
    end

    return ans
end

function CarlsonAlg.E(φ::Reactant.TracedRNumber, m::Real)
    return _reactant_internal_E(φ, m)
end

function CarlsonAlg.E(φ::Real, m::Reactant.TracedRNumber)
    return _reactant_internal_E(φ, m)
end

function CarlsonAlg.E(φ::Reactant.TracedRNumber, m::Reactant.TracedRNumber)
    return _reactant_internal_E(φ, m)
end

function CarlsonAlg.E(
    phi::Reactant.AnyTracedRArray{Tp,N},
    m::Real,
) where {Tp,N}
    return _reactant_internal_E.(phi, m)
end

function CarlsonAlg.E(
    phi::Reactant.AnyTracedRArray{Tp,N},
    m::Reactant.TracedRNumber,
) where {Tp,N}
    return _reactant_internal_E.(phi, m)
end

end
