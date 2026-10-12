module JacobiEllipticReactantExt

using JacobiElliptic
using Reactant
import JacobiElliptic: CarlsonAlg, ArithmeticGeometricMeanAlg, StaticArrays

# Number of ARM iterations needed for conversion.
_reactant_agm_iters(::Type{<:Any}) = 8
_reactant_agm_iters(::Type{Float32}) = 5

@inline function _reactant_ellipke_base(m::T) where {T}
    a = one(T)
    b = sqrt(one(T) - m)
    sum_c2 = m / 2
    pow2 = one(T)
    for _ in 1:_reactant_agm_iters(Reactant.unwrapped_eltype(T))
        c = (a - b) / 2
        sum_c2 = muladd(pow2, c * c, sum_c2)
        pow2 += pow2
        a, b = (a + b) / 2, sqrt(a * b)
    end
    k = T(π / 2) / a
    e = k * (one(T) - sum_c2)
    return (k, e)
end

#----------------------------------------------------------------------------------------
# Elliptic K(m), E(m)
#----------------------------------------------------------------------------------------

# This code is purposefully branch-free. 
# Every case maps m into [0, 1) and runs the AGM algorithm.
# NaN inputs propagate through the arithmetic.
Reactant.@reactant_overlay function ArithmeticGeometricMeanAlg.ellipke(m::T) where {T}
    oneT = one(T)
    isneg = m < zero(T)
    isbig = m > oneT
    m_base = Base.ifelse(isneg, -m / (oneT - m), Base.ifelse(isbig, inv(m), m))
    # The reduction can land on exactly 1 for an input that is not 1: -m/(1 - m) rounds to
    # 1.0 once m <= -9.1e15. src/AGM.jl catches that through its recursion, so the reduced
    # parameter has to be tested too, or the AGM runs with b == 0 and returns garbage.
    isone = (m == oneT) | (m_base == oneT)
    # keep sqrt(1 - m_base) away from 0 at m == 1 so the unselected branch has finite gradients
    m_base = Base.ifelse(isone, zero(T), m_base)
    k, e = _reactant_ellipke_base(m_base)
    k = Base.ifelse(isone, T(Inf), k)
    e = Base.ifelse(isone, oneT, e)
    sq = sqrt(Base.ifelse(isneg, oneT - m, Base.ifelse(isbig, m, oneT)))
    K = Base.ifelse(isneg | isbig, k / sq, k)
    E = Base.ifelse(isneg, e * sq, Base.ifelse(isbig, sq * (e - (oneT - m_base) * k), e))
    return (K, E)
end
#----------------------------------------------------------------------------------------
# Elliptic K(m), E(m) 
#----------------------------------------------------------------------------------------

Reactant.@reactant_overlay ArithmeticGeometricMeanAlg.K(m::T) where {T} =
    ArithmeticGeometricMeanAlg.ellipke(m)[1]

Reactant.@reactant_overlay ArithmeticGeometricMeanAlg.E(m::T) where {T} =
    ArithmeticGeometricMeanAlg.ellipke(m)[2]

#----------------------------------------------------------------------------------------
# Elliptic F(ϕ, m)
#----------------------------------------------------------------------------------------

function _reactant_DRF_ifbody(X::A, Y::B, Z::C, ERRTOL::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    C1 = T(1 / 24)
    C2 = T(3 / 44)
    C3 = T(1 / 14)
    C0 = T(1 / 10)
    inv3 = T(1 / 3)
    inv4 = T(1 / 4)

    XN = copy(X)
    YN = copy(Y)
    ZN = copy(Z)
    MU = (XN + YN + ZN) * inv3
    ninvMU = -inv(MU)
    XNDEV = muladd(ninvMU, (MU + XN), 2)
    YNDEV = muladd(ninvMU, (MU + YN), 2)
    ZNDEV = muladd(ninvMU, (MU + ZN), 2)
    EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))

    Reactant.@trace while EPSLON >= ERRTOL
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        YNROOTZNROOT = YNROOT * ZNROOT
        LAMDA = muladd(XNROOT, YNROOT + ZNROOT, YNROOTZNROOT)
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4

        MU = (XN + YN + ZN) * inv3
        ninvMU = -inv(MU)
        XNDEV = muladd(ninvMU, (MU + XN), 2)
        YNDEV = muladd(ninvMU, (MU + YN), 2)
        ZNDEV = muladd(ninvMU, (MU + ZN), 2)
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
    end

    XNDEVYNDEV = XNDEV * YNDEV
    E2 = muladd(-ZNDEV, ZNDEV, XNDEVYNDEV)
    E3 = XNDEVYNDEV * ZNDEV
    S = one(T) + muladd(E2, muladd(-C2, E3, muladd(C1, E2, -C0)), C3 * E3)
    return S / sqrt(MU)
end

Reactant.@reactant_overlay function CarlsonAlg.DRF(X::A, Y::B, Z::C) where {A,B,C}
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
    minxyz = min(X, Y, Z) < zeroT
    maxxyz = max(X, Y, Z) > UPLIM
    mxyxzyz = min(X + Y, X + Z, Y + Z) < LOLIM
    err = Base.ifelse(minxyz, 1, Base.ifelse(maxxyz, 3, Base.ifelse(mxyxzyz, 2, 0)))

    @trace if err == 0
        ans = _reactant_DRF_ifbody(X, Y, Z, ERRTOL)
    end

    return (ans, err)
end

function _reactant_rawF(sinphi::A, m::B) where {A,B}
    T = promote_type(A, B)
    sinphi = T(sinphi)
    m = T(m)
    oneT = one(T)
    infT = T(Inf)
    sinphi2 = sinphi * sinphi
    drf, _ = CarlsonAlg.DRF(oneT - sinphi2, muladd(-m, sinphi2, oneT), oneT)
    return Base.ifelse(
        (abs(sinphi) == oneT) & (m == oneT),
        sign(sinphi) * infT,
        sinphi * drf,
    )
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
        ans =
            2 * fld(phi2, piT) * CarlsonAlg.K(m) -
            _reactant_rawF(cos(mod(phi2, piT)), m)
    else
        ans = _reactant_rawF(sin(phi), m)
    end

    return ans
end

Reactant.@reactant_overlay function CarlsonAlg.F(phi::A, m::B) where {A,B}
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
        ans =
            m12 * CarlsonAlg.K(m1m) -
            sign(newphi) * m12 * _reactant_internal_F(abs(newphi), m1m)
    else
        ans = sign(phi) * _reactant_internal_F(abs(phi), m)
    end

    return ans
end

function __reactant_am(u::A, m::B, tol::C) where {A,B,C}

    T = promote_type(A, B, C)

    # Use an immutable local buffer so GPU backends can keep this in registers.
    _ambuf = zeros(T, 10)

    sqrt_tol = sqrt(tol)
    zeroT = zero(T)
    t = tanh(u)
    m1 = one(T) - m
    flag = true
    ans, flag = Base.ifelse(
        iszero(u),
        (zeroT, false),
        Base.ifelse(
            m < sqrt_tol,
            (u - m * (u - sin(2 * u) / 2) / 4, false),
            Base.ifelse(
                m1 < sqrt_tol,
                (asin(t) + m1 * (t - u * (one(T) - t^2)) * cosh(u) / 4, false),
                (zeroT, true),
            ),
        ),
    )

    a = one(T)
    b = sqrt(m1)
    c = sqrt(m)
    n = zero(T)
    @trace track_numbers=true while (n < 10)
        a, b, c, n = ((a + b) / 2, sqrt(a * b), (a - b) / 2, n + 1)
        Reactant.@allowscalar _ambuf[unsafe_trunc(Int, n)] = c / a
    end

    phi = a * u * (2^n)
    @trace track_numbers=true while (n > 0) & flag
        temp = Reactant.@allowscalar _ambuf[unsafe_trunc(Int, n)]
        phi = (phi + asin(temp * sin(phi))) / 2
        n -= 1
    end

    # The AGM is only used on the `flag` path, so its non-convergence check must
    # not reach the closed forms above: the m1 < sqrt_tol branch sets m1 == 0 for
    # m == 1, which the AGM can never converge for.
    return Base.ifelse(flag, Base.ifelse(abs(c) > tol, T(NaN), phi), ans)
end

function __reactant_am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    return __reactant_am(u, m, eps(T))
end

"""
    am(u::Real, m::Real)

Returns amplitude, φ, such that u = F(φ | m)
"""
Reactant.@reactant_overlay function CarlsonAlg.am(u::A, m::B) where {A,B}
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

for (f, a, b, c) in (
    (:sn, :(sin(phi)), :(sqrtmu1 * s), :(sqrt(mu) * sin(phi))),
    (:cn, :(cos(phi)), :(cos(phi)), :(sqrt(1 - mu * sin(phi)^2))),
    (:dn, :(sqrt(1 - m * sin(phi)^2)), :(one(T)), :(cos(phi))),
)
    @eval begin
        Reactant.@reactant_overlay function CarlsonAlg.$f(u::A, m::B) where {A,B}
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
    end
end

#----------------------------------------------------------------------------------------
# Incomplete Elliptic E(ϕ, m)
#----------------------------------------------------------------------------------------

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

    XN = copy(X)
    YN = copy(Y)
    ZN = copy(Z)
    SIGMA = zero(T)
    POWER4 = one(T)

    MU = (XN + YN + threeT * ZN) * inv5
    invMU = inv(MU)
    XNDEV = (MU - XN) * invMU
    YNDEV = (MU - YN) * invMU
    ZNDEV = (MU - ZN) * invMU
    EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))

    Reactant.@trace while EPSLON >= ERRTOL
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        ynrootznroot = YNROOT * ZNROOT
        LAMDA = muladd(XNROOT, YNROOT + ZNROOT, ynrootznroot)
        SIGMA = muladd(POWER4, inv(ZNROOT * (ZN + LAMDA)), SIGMA)
        POWER4 *= inv4
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4

        MU = (XN + YN + threeT * ZN) * inv5
        invMU = inv(MU)
        XNDEV = (MU - XN) * invMU
        YNDEV = (MU - YN) * invMU
        ZNDEV = (MU - ZN) * invMU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
    end

    EA = XNDEV * YNDEV
    EB = ZNDEV * ZNDEV
    EC = EA - EB
    ED = EA - sixT * EB
    EF = ED + EC + EC
    S1 = ED * (-C1 + C3 * ED * inv4 - threeT * C4 * ZNDEV * EF / 2)
    S2 = ZNDEV * (C2 * EF + ZNDEV * (-C3 * EC + ZNDEV * C4 * EA))
    return threeT * SIGMA + POWER4 * (oneT + S1 + S2) / (MU * sqrt(MU))
end

Reactant.@reactant_overlay function CarlsonAlg.DRD(X::A, Y::B, Z::C) where {A,B,C}
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
    minxy = min(X, Y) < zeroT
    maxxyz = max(X, Y, Z) > UPLIM
    maxypz = min(X + Y, Z) < LOLIM
    err = Base.ifelse(minxy, 1, Base.ifelse(maxxyz, 3, Base.ifelse(maxypz, 2, 0)))

    ans = Base.ifelse(err == 0, _reactant_DRD_ifbody(X, Y, Z), ans)

    return (ans, err)
end

#----------------------------------------------------------------------------------------
# Complete Carlson K(m), E(m)
#----------------------------------------------------------------------------------------
function _reactant_complete_reduce_m(m::T) where {T}
# Maps m onto the [0, 1) for complete K and E
    oneT = one(T)
    zeroT = zero(T)
    reduced = Base.ifelse(
        m < zeroT,
        let m_abs = -m
            m_abs / (m_abs + oneT)
        end,
        Base.ifelse(m > oneT, inv(m), m),
    )
    return Base.ifelse(isnan(m) | (m == oneT), zeroT, reduced)
end

Reactant.@reactant_overlay function CarlsonAlg.K(m::T) where {T}
    oneT = one(T)
    zeroT = zero(T)
    m_eval = _reactant_complete_reduce_m(m)
    drf, _ = CarlsonAlg.DRF(zeroT, oneT - m_eval, oneT)

    return Base.ifelse(
        isnan(m),
        T(NaN),
        Base.ifelse(
            m == oneT,
            T(Inf),
            Base.ifelse(
                m < zeroT,
                drf / sqrt(oneT + abs(m)),
                Base.ifelse(m > oneT, drf / sqrt(m), drf),
            ),
        ),
    )
end

Reactant.@reactant_overlay function CarlsonAlg.ellipke(m::T) where {T}
    oneT = one(T)
    zeroT = zero(T)
    nanT = T(NaN)
    inv3 = T(1 / 3)
    m_eval = _reactant_complete_reduce_m(m)
    y = oneT - m_eval
    drf, _ = CarlsonAlg.DRF(zeroT, y, oneT)
    drd, _ = CarlsonAlg.DRD(zeroT, y, oneT)
    e_base = muladd(-m_eval * inv3, drd, drf)
    sqrt_factor = sqrt(oneT + abs(m))

    k = Base.ifelse(
        isnan(m),
        nanT,
        Base.ifelse(
            m == oneT,
            T(Inf),
            Base.ifelse(
                m < zeroT,
                drf / sqrt_factor,
                Base.ifelse(m > oneT, drf / sqrt(m), drf),
            ),
        ),
    )
    e = Base.ifelse(
        isnan(m),
        nanT,
        Base.ifelse(
            m == oneT,
            oneT,
            Base.ifelse(
                m < zeroT,
                e_base * sqrt_factor,
                Base.ifelse(
                    m > oneT,
                    sqrt(m) * muladd(-(oneT - inv(m)), drf, e_base),
                    e_base,
                ),
            ),
        ),
    )
    return (k, e)
end

function _reactant_rawE(sinphi::A, m::B) where {A,B}
    T = promote_type(A, B)
    sinphi = T(sinphi)
    m = T(m)
    inv3 = T(1 / 3)
    nanT = T(NaN)
    ans = zero(T)

    Reactant.@trace if m > one(T)
        minv = inv(m)
        m12 = sqrt(m)
        sinalpha = m12 * sinphi
        sinalpha2 = sinalpha * sinalpha
        cosalpha2 = one(T) - sinalpha2
        yalpha = muladd(-minv, sinalpha2, one(T))
        drf, ierr1 = CarlsonAlg.DRF(cosalpha2, yalpha, one(T))
        drd, ierr2 = CarlsonAlg.DRD(cosalpha2, yalpha, one(T))

        ans = Base.ifelse(
            (ierr1 == 0) & (ierr2 == 0),
            sinalpha * (drf - sinalpha2 * inv3 * drd) / m12,
            Base.ifelse((ierr1 == 2) & (ierr2 == 2), sinphi, nanT),
        )
    else
        sinphi2 = sinphi * sinphi
        cosphi2 = one(T) - sinphi2
        y = muladd(-m, sinphi2, one(T))
        drf, ierr1 = CarlsonAlg.DRF(cosphi2, y, one(T))
        drd, ierr2 = CarlsonAlg.DRD(cosphi2, y, one(T))

        ans = Base.ifelse(
            (ierr1 == 0) & (ierr2 == 0),
            sinphi * (drf - m * sinphi2 * inv3 * drd),
            Base.ifelse((ierr1 == 2) & (ierr2 == 2), sinphi, nanT),
        )
    end

    return ans
end

Reactant.@reactant_overlay function CarlsonAlg.E(phi::A, m::B) where {A,B}
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
        ans =
            2 * fld(phi2, piT) * CarlsonAlg.E(m) -
            _reactant_rawE(cos(mod(phi2, piT)), m)
    else
        ans = _reactant_rawE(sin(phi), m)
    end

    return ans
end


#----------------------------------------------------------------------------------------
# Incomplete Elliptic Pi(n, ϕ, m)
#----------------------------------------------------------------------------------------

function _reactant_DRC_ifbody(X::A, Y::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)
    twoT = T(2)
    inv3 = T(1 / 3)
    inv4 = T(1 / 4)
    C1 = T(1 / 7)
    C2 = T(9 / 22)
    C3 = T(3 / 8)
    C4 = T(3 / 10)
    ERRTOL = (eps(T) / 32)^T(1 / 6)

    XN = X
    YN = Y
    MU = zero(T)
    SN = zero(T)

    MU = (XN + YN + YN) * inv3
    invMU = inv(MU)
    SN = muladd(invMU, YN + MU, -twoT)

    @trace while abs(SN) >= ERRTOL
        LAMDA = muladd(twoT * sqrt(XN), sqrt(YN), YN)
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4

        MU = (XN + YN + YN) * inv3
        invMU = inv(MU)
        SN = muladd(invMU, YN + MU, -twoT)
    end

    S = SN^2 * muladd(SN, muladd(SN, muladd(SN, C2, C3), C1), C4)
    return (oneT + S) / sqrt(MU)
end

function _reactant_DRC(X::A, Y::B) where {A,B}
    T = promote_type(A, B)
    Tc = Reactant.unwrapped_eltype(T)
    X = T(X)
    Y = T(Y)
    zeroT = zero(T)

    LOLIM = T(5) * floatmin(Tc)
    UPLIM = T(floatmax(Tc) / 5)

    ans = zeroT
    err = Base.ifelse(
        (X < zeroT) | (Y <= zeroT),
        1,
        Base.ifelse(max(X, Y) > UPLIM, 3, Base.ifelse(X + Y < LOLIM, 2, 0)),
    )

    @trace if err == 0
        ans = _reactant_DRC_ifbody(X, Y)
    end

    return (ans, err)
end
function _reactant_DRJ_ifbody(X::A, Y::B, Z::C, P::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    oneT = one(T)
    twoT = T(2)
    threeT = T(3)
    inv4 = T(1 / 4)
    inv5 = T(1 / 5)
    ERRTOL = (eps((T)) / 6)^T(1 / 6)
    C1 = T(3 / 14)
    C2 = T(1 / 3)
    C3 = T(3 / 22)
    C4 = T(3 / 26)

    XN = copy(X)
    YN = copy(Y)
    ZN = copy(Z)
    PN = copy(P)
    SIGMA = zero(T)
    POWER4 = one(T)
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)
    PNDEV = zero(T)
    IER = zero(Int)

    XNYNZN = XN + YN + ZN
    MU = (XNYNZN + twoT * PN) * inv5
    invMU = inv(MU)
    XNDEV = (MU - XN) * invMU
    YNDEV = (MU - YN) * invMU
    ZNDEV = (MU - ZN) * invMU
    PNDEV = (MU - PN) * invMU
    EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV), abs(PNDEV))

    @trace while EPSLON >= ERRTOL
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        YNROOTZNROOT = YNROOT * ZNROOT
        rootsum = XNROOT + YNROOT + ZNROOT
        LAMDA = muladd(XNROOT, (YNROOT + ZNROOT), YNROOTZNROOT)
        pn_plus_lamda = PN + LAMDA
        alpha_base = muladd(PN, rootsum, XNROOT * YNROOTZNROOT)
        ALFA = alpha_base * alpha_base
        BETA = PN * pn_plus_lamda * pn_plus_lamda
        drc, IER = _reactant_DRC(ALFA, BETA)
        SIGMA = muladd(POWER4, drc, SIGMA)
        POWER4 *= inv4
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4
        ZN = (ZN + LAMDA) * inv4
        PN = pn_plus_lamda * inv4

        XNYNZN = XN + YN + ZN
        MU = (XNYNZN + twoT * PN) * inv5
        invMU = inv(MU)
        XNDEV = (MU - XN) * invMU
        YNDEV = (MU - YN) * invMU
        ZNDEV = (MU - ZN) * invMU
        PNDEV = (MU - PN) * invMU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV), abs(PNDEV))
    end

    YNDEVZNDEV = YNDEV * ZNDEV
    EA = muladd(XNDEV, YNDEV + ZNDEV, YNDEVZNDEV)
    EB = XNDEV * YNDEVZNDEV
    EC = PNDEV * PNDEV
    E2 = EA - threeT * EC
    E3 = muladd(twoT * PNDEV, EA - EC, EB)
    S1 = oneT + E2 * (-C1 + (threeT * C3 / 4) * E2 - (threeT * C4 / 2) * E3)
    S2 = EB * (C2 / twoT + PNDEV * (-twoT * C3 + PNDEV * C4))
    S3 = PNDEV * EA * (C2 - PNDEV * C3) - C2 * PNDEV * EC
    ans = threeT * SIGMA + POWER4 * (S1 + S2 + S3) / (MU * sqrt(MU))

    return (ans, IER)
end

function _reactant_DRJ(X::A, Y::B, Z::C, P::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    Tc = Reactant.unwrapped_eltype(T)
    X = T(X)
    Y = T(Y)
    Z = T(Z)
    P = T(P)
    zeroT = zero(T)

    LOLIM = T((5 * floatmin(Tc))^Tc(1 / 3))
    UPLIM = T(3 / 10) * T((floatmax(Tc) / 5)^Tc(1 / 3))

    ans = zeroT
    err = Base.ifelse(
        min(X, Y, Z) < zeroT,
        1,
        Base.ifelse(
            max(X, Y, Z, P) > UPLIM,
            3,
            Base.ifelse(min(X + Y, X + Z, Y + Z, P) < LOLIM, 2, 0),
        ),
    )

    @trace if err == 0
        ans, err = _reactant_DRJ_ifbody(X, Y, Z, P)
    end

    return (ans, err)
end
@inline function _reactant_custom_atanh(a::T) where {T}
    oneT = one(T)
    return log(abs(oneT + a) / abs(oneT - a)) / 2
end

@inline function _reactant_FukushimaT(t::A, h::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)

    return Base.ifelse(
        h > zero(T),
        atan(t * sqrt(h)) / sqrt(h),
        Base.ifelse(
            h == zero(T),
            t,
            begin
                sqrt_neg_h = sqrt(-h)
                arg = t * sqrt_neg_h
                Base.ifelse(abs(arg) < oneT, atanh(arg), _reactant_custom_atanh(arg)) /
                sqrt_neg_h
            end,
        ),
    )
end

function _reactant_cel(kc::A, p::B, a::C, b::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    ca = eps(Reactant.unwrapped_eltype(T))
    oneT = one(T)
    twoT = T(2)
    zeroT = zero(T)
    pi_over_2 = T(π / 2)

    kc = abs(kc)
    p = T(p)
    a = T(a)
    b = T(b)

    e = copy(kc)
    m = one(T)
    f = zero(T)
    g = zero(T)
    q = zero(T)

    Reactant.@trace if p > zeroT
        p = sqrt(p)
        b = b / p
    else
        f = kc^2
        q = oneT - f
        g = oneT - p
        f = f - p
        q = (b - a * p) * q
        p = sqrt(f / g)
        a = (a - b) / g
        b = -q * (g^2 * p) + a * p
    end

    count = 0
    f = copy(a)
    invp = inv(p)
    a = muladd(invp, b, a)
    g = e * invp
    b = twoT * muladd(f, g, b)
    p = g + p
    g = copy(m)
    m = kc + m

    Reactant.@trace while (abs(g - kc) >= g * ca) & (count < 1000)
        kc = twoT * √e
        e = kc * m

        f = copy(a)
        invp = inv(p)
        a = muladd(invp, b, a)
        g = e * invp
        b = twoT * muladd(f, g, b)
        p = g + p
        g = copy(m)
        m = kc + m
        count += 1
    end

    return pi_over_2 * muladd(a, m, b) / (m * (m + p))
end

function _reactant_rawPi(n::A, sinphi::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n = T(n)
    sinphi = T(sinphi)
    m = T(m)
    oneT = one(T)
    inv3 = T(1 / 3)
    infT = T(Inf)
    nanT = T(NaN)

    sinphi2 = sinphi * sinphi
    cosphi2 = oneT - sinphi2
    y = muladd(-m, sinphi2, oneT)
    p = muladd(-n, sinphi2, oneT)
    drf, ierr1 = CarlsonAlg.DRF(cosphi2, y, oneT)
    drj, ierr2 = _reactant_DRJ(cosphi2, y, oneT, p)

    return Base.ifelse(
        (ierr1 == 0) & (ierr2 == 0),
        sinphi * (drf + n * sinphi2 * inv3 * drj),
        Base.ifelse(
            ((ierr1 == 2) & (ierr2 == 2)) | ((ierr1 == 0) & (ierr2 == 2)),
            infT,
            nanT,
        ),
    )
end

function _reactant_incomplete_Pi_core(n::A, phi::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n = T(n)
    phi = T(phi)
    m = T(m)
    piT = T(π)
    halfpi = T(π / 2)
    nanT = T(NaN)
    ans = zero(T)

    Reactant.@trace if isnan(n) | isnan(phi) | isnan(m)
        ans = nanT
    elseif 2 * abs(phi) > piT
        phi2 = phi + halfpi
        ans = 2 * fld(phi2, piT) * Pi(n, m) - _reactant_rawPi(n, cos(mod(phi2, piT)), m)
    else
        ans = _reactant_rawPi(n, sin(phi), m)
    end

    return ans
end

function _reactant_incomplete_Pi_nonneg(n::A, phi::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n = T(n)
    phi = T(phi)
    m = T(m)
    oneT = one(T)
    ans = zero(T)

    Reactant.@trace if n > oneT
        nc = oneT - n
        sinphi, cosphi = sincos(phi)
        sinphi2 = sinphi * sinphi
        t1 = sinphi / (cosphi * sqrt(muladd(-m, sinphi2, oneT)))
        h1 = nc * (n - m) / n
        n1 = m / n
        ans =
            _reactant_FukushimaT(t1, h1) - _reactant_incomplete_Pi_core(n1, phi, m) +
            F(phi, m)
    else
        ans = _reactant_incomplete_Pi_core(n, phi, m)
    end

    return ans
end

Reactant.@reactant_overlay function CarlsonAlg.Pi(n::A, phi::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n = T(n)
    phi = T(phi)
    m = T(m)
    oneT = one(T)
    ans = zero(T)

    Reactant.@trace if (m < zero(T)) & (n == m)
        # The imaginary modulus transformation below is 0/0 at n == m, so use the closed
        # form Pi(m; phi | m) = (E(phi | m) - m*sin*cos/sqrt(1 - m*sin^2)) / (1 - m).
        sinphi, cosphi = sincos(phi)
        ans =
            (
                CarlsonAlg.E(phi, m) -
                m * sinphi * cosphi / sqrt(muladd(-m, sinphi * sinphi, oneT))
            ) / (oneT - m)
    elseif m < zero(T)
        # Pi(n; phi + k*pi | m) = Pi(n; phi | m) + 2k Pi(n | m).
        piT = T(Base.π)
        halfpi = T(Base.π / 2)
        phi2 = phi + halfpi
        wrapped = 2 * abs(phi) > piT
        phi_r = Base.ifelse(wrapped, mod(phi2, piT) - halfpi, phi)
        offset = Base.ifelse(wrapped, 2 * fld(phi2, piT) * CarlsonAlg.Pi(n, m), zero(T))

        mc = oneT - m
        imc = inv(mc)
        mN = -m * imc
        sinphi = sin(phi_r)
        sinphi2 = sinphi * sinphi
        phiN = asin(sqrt(mc / muladd(-m, sinphi2, oneT)) * sinphi)
        nN = (n - m) * imc

        ans =
            offset +
            sqrt(imc) / nN *
            (mN * F(phiN, mN) + imc * n * _reactant_incomplete_Pi_nonneg(nN, phiN, mN))
    else
        ans = _reactant_incomplete_Pi_nonneg(n, phi, m)
    end

    return ans
end

#----------------------------------------------------------------------------------------
# Complete Elliptic Pi(n, m)
#----------------------------------------------------------------------------------------

function _reactant_complete_Pi_core(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    n = T(n)
    m = T(m)
    oneT = one(T)
    zeroT = zero(T)
    infT = T(Inf)
    ans = zeroT

    Reactant.@trace if n == oneT
        # 1/(1 - n*sin(θ)^2) pole at θ = π/2.
        ans = infT
    elseif n == zeroT
        ans = CarlsonAlg.K(m)
    elseif m == oneT
        # Only m == 1 diverges. m == 0 falls through to `cel`, which gives the
        # finite pi/(2*sqrt(1 - n)).
        ans = infT
    else
        ans = _reactant_cel(sqrt(oneT - m), oneT - n, oneT, oneT)
    end

    return ans
end

Reactant.@reactant_overlay function CarlsonAlg.Pi(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    n = T(n)
    m = T(m)
    oneT = one(T)
    nanT = T(NaN)
    ans = zero(T)

    Reactant.@trace if m > oneT
        ans = nanT
    elseif n > oneT
        ans = CarlsonAlg.K(m) - _reactant_complete_Pi_core(m / n, m)
    else
        ans = _reactant_complete_Pi_core(n, m)
    end

    return ans
end




#----------------------------------------------------------------------------------------
# Jacobi Elliptic Functions
#----------------------------------------------------------------------------------------


for (p, num) in ((:s, :(sn(u, m))), (:c, :(cn(u, m))), (:d, :(dn(u, m))), (:n, 1)),
    (q, den) in ((:s, :(sn(u, m))), (:c, :(cn(u, m))), (:d, :(dn(u, m))), (:n, 1))

    f = Symbol(p, q)

    if p == q
        @eval begin
            Reactant.@reactant_overlay CarlsonAlg.$f(u::T, m::S) where {T,S} =
                one(promote_type(T, S))
        end
    elseif q != :n
        @eval begin
            Reactant.@reactant_overlay CarlsonAlg.$f(u::T, m::S) where {T,S} =
                ($num) / ($den)
        end
    end
end
end
