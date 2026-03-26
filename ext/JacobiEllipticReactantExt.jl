module JacobiEllipticReactantExt

using JacobiElliptic
using Reactant
import JacobiElliptic: CarlsonAlg, ArithmeticGeometricMeanAlg, StaticArrays

function _am_buffer(::T) where {T}
    zeroT = zero(T)
    #return StaticArrays.@MVector
    return [zeroT, zeroT, zeroT, zeroT, zeroT, zeroT, zeroT, zeroT, zeroT, zeroT, zeroT]
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
        for i = 10:-1:1
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
        CarlsonAlg.$f(u::Reactant.TracedRNumber, m::Reactant.TracedRNumber) =
            $reactant_f(u, m)
    end
end

xn = ((:s, :(sn(u, m))), (:c, :(cn(u, m))), (:d, :(dn(u, m))), (:n, 1))
for (p, num) in xn, (q, den) in xn
    f = Symbol(p, q)

    if p == q
        @eval begin
            CarlsonAlg.$f(::Reactant.TracedRNumber{T}, ::S) where {T,S<:Real} =
                one(promote_type(T, S))
            CarlsonAlg.$f(::S, ::Reactant.TracedRNumber{T}) where {S<:Real,T} =
                one(promote_type(S, T))
            CarlsonAlg.$f(
                ::Reactant.TracedRNumber{T1},
                ::Reactant.TracedRNumber{T2},
            ) where {T1,T2} = one(promote_type(T1, T2))
        end
    elseif q != :n
        @eval begin
            CarlsonAlg.$f(u::Reactant.TracedRNumber, m::Real) = ($num) / ($den)
            CarlsonAlg.$f(u::Real, m::Reactant.TracedRNumber) = ($num) / ($den)
            CarlsonAlg.$f(u::Reactant.TracedRNumber, m::Reactant.TracedRNumber) =
                ($num) / ($den)
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

    for _ = 1:12
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
    ERRTOL = T((eps(Reactant.unwrapped_eltype(T)) / 32)^Reactant.unwrapped_eltype(T)(1 / 6))

    XN = X
    YN = Y
    MU = zero(T)
    SN = zero(T)
    active = true

    for _ = 1:8
        next_MU = (XN + YN + YN) * inv3
        invMU = inv(next_MU)
        next_SN = muladd(invMU, YN + next_MU, -twoT)
        converged = abs(next_SN) < ERRTOL
        LAMDA = muladd(twoT * sqrt(XN), sqrt(YN), YN)
        next_XN = (XN + LAMDA) * inv4
        next_YN = (YN + LAMDA) * inv4
        continue_active = active & !converged

        MU = Base.ifelse(active, next_MU, MU)
        SN = Base.ifelse(active, next_SN, SN)
        XN = Base.ifelse(continue_active, next_XN, XN)
        YN = Base.ifelse(continue_active, next_YN, YN)
        active = continue_active
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

CarlsonAlg.DRD(X::Reactant.TracedRNumber, Y::Real, Z::Real) = _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(X::Real, Y::Reactant.TracedRNumber, Z::Real) = _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(X::Real, Y::Real, Z::Reactant.TracedRNumber) = _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(X::Reactant.TracedRNumber, Y::Reactant.TracedRNumber, Z::Real) =
    _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(X::Reactant.TracedRNumber, Y::Real, Z::Reactant.TracedRNumber) =
    _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(X::Real, Y::Reactant.TracedRNumber, Z::Reactant.TracedRNumber) =
    _reactant_DRD(X, Y, Z)
CarlsonAlg.DRD(
    X::Reactant.TracedRNumber,
    Y::Reactant.TracedRNumber,
    Z::Reactant.TracedRNumber,
) = _reactant_DRD(X, Y, Z)

function _reactant_DRJ_ifbody(X::A, Y::B, Z::C, P::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    oneT = one(T)
    twoT = T(2)
    threeT = T(3)
    inv4 = T(1 / 4)
    inv5 = T(1 / 5)
    ERRTOL = T((eps(Reactant.unwrapped_eltype(T)) / 6)^Reactant.unwrapped_eltype(T)(1 / 6))
    C1 = T(3 / 14)
    C2 = T(1 / 3)
    C3 = T(3 / 22)
    C4 = T(3 / 26)

    XN = X
    YN = Y
    ZN = Z
    PN = P
    SIGMA = zero(T)
    POWER4 = oneT
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)
    PNDEV = zero(T)
    IER = zero(Int)
    active = true

    for _ = 1:10
        XNYNZN = XN + YN + ZN
        next_MU = (XNYNZN + twoT * PN) * inv5
        invMU = inv(next_MU)
        next_XNDEV = (next_MU - XN) * invMU
        next_YNDEV = (next_MU - YN) * invMU
        next_ZNDEV = (next_MU - ZN) * invMU
        next_PNDEV = (next_MU - PN) * invMU
        EPSLON = max(abs(next_XNDEV), abs(next_YNDEV), abs(next_ZNDEV), abs(next_PNDEV))
        converged = EPSLON < ERRTOL

        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        YNROOTZNROOT = YNROOT * ZNROOT
        rootsum = XNROOT + YNROOT + ZNROOT
        LAMDA = muladd(XNROOT, YNROOT + ZNROOT, YNROOTZNROOT)
        pn_plus_lamda = PN + LAMDA
        alpha_base = muladd(PN, rootsum, XNROOT * YNROOTZNROOT)
        ALFA = alpha_base * alpha_base
        BETA = PN * pn_plus_lamda * pn_plus_lamda
        drc, next_IER = _reactant_DRC(ALFA, BETA)
        continue_active = active & !converged

        MU = Base.ifelse(active, next_MU, MU)
        XNDEV = Base.ifelse(active, next_XNDEV, XNDEV)
        YNDEV = Base.ifelse(active, next_YNDEV, YNDEV)
        ZNDEV = Base.ifelse(active, next_ZNDEV, ZNDEV)
        PNDEV = Base.ifelse(active, next_PNDEV, PNDEV)
        SIGMA = Base.ifelse(continue_active, muladd(POWER4, drc, SIGMA), SIGMA)
        POWER4 = Base.ifelse(continue_active, POWER4 * inv4, POWER4)
        XN = Base.ifelse(continue_active, (XN + LAMDA) * inv4, XN)
        YN = Base.ifelse(continue_active, (YN + LAMDA) * inv4, YN)
        ZN = Base.ifelse(continue_active, (ZN + LAMDA) * inv4, ZN)
        PN = Base.ifelse(continue_active, pn_plus_lamda * inv4, PN)
        IER = Base.ifelse(continue_active, next_IER, IER)
        active = continue_active
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

    for _ = 1:10
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
        Base.ifelse(
            max(X, Y, Z) > UPLIM,
            3,
            Base.ifelse(min(X + Y, X + Z, Y + Z) < LOLIM, 2, 0),
        ),
    )

    #ans = Base.ifelse(err == 0, _reactant_valid_DRF(X, Y, Z, ERRTOL), ans)
    @trace if err == 0
        ans = _reactant_DRF_ifbody(X, Y, Z, ERRTOL)
    end

    return (ans, err)
end

CarlsonAlg.DRF(X::Reactant.TracedRNumber, Y::Real, Z::Real) = _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(X::Real, Y::Reactant.TracedRNumber, Z::Real) = _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(X::Real, Y::Real, Z::Reactant.TracedRNumber) = _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(X::Reactant.TracedRNumber, Y::Reactant.TracedRNumber, Z::Real) =
    _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(X::Reactant.TracedRNumber, Y::Real, Z::Reactant.TracedRNumber) =
    _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(X::Real, Y::Reactant.TracedRNumber, Z::Reactant.TracedRNumber) =
    _reactant_DRF(X, Y, Z)
CarlsonAlg.DRF(
    X::Reactant.TracedRNumber,
    Y::Reactant.TracedRNumber,
    Z::Reactant.TracedRNumber,
) = _reactant_DRF(X, Y, Z)

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

    kc = abs(T(kc))
    p = T(p)
    a = T(a)
    b = T(b)

    e = kc
    m = oneT
    f = zeroT
    g = zeroT
    q = zeroT

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

    active = true
    for _ = 1:32
        current_f = a
        invp = inv(p)
        next_a = muladd(invp, b, a)
        next_g = e * invp
        next_b = twoT * muladd(current_f, next_g, b)
        next_p = next_g + p
        current_m = m
        next_m = kc + m
        converged = abs(current_m - kc) < current_m * ca
        continue_active = active & !converged
        next_kc = twoT * sqrt(e)
        next_e = next_kc * next_m

        a = Base.ifelse(active, next_a, a)
        b = Base.ifelse(active, next_b, b)
        p = Base.ifelse(active, next_p, p)
        m = Base.ifelse(active, next_m, m)
        kc = Base.ifelse(continue_active, next_kc, kc)
        e = Base.ifelse(continue_active, next_e, e)
        active = continue_active
    end

    return pi_over_2 * muladd(a, m, b) / (m * (m + p))
end

function _reactant_ellipke_base(m::T) where {T}
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
    for _ = 1:8
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

function _reactant_K(m::T) where {T}
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
function _reactant_E(m::T) where {T}
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
ArithmeticGeometricMeanAlg.ellipke(m::Reactant.TracedRNumber) =
    (_reactant_K(m), _reactant_E(m))

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

function CarlsonAlg.E(phi::Reactant.AnyTracedRArray{Tp,N}, m::Real) where {Tp,N}
    return _reactant_internal_E.(phi, m)
end

function CarlsonAlg.E(
    phi::Reactant.AnyTracedRArray{Tp,N},
    m::Reactant.TracedRNumber,
) where {Tp,N}
    return _reactant_internal_E.(phi, m)
end

#----------------------------------------------------------------------------------------
# Incomplete Elliptic Pi(n, ϕ, m)
#----------------------------------------------------------------------------------------

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
    drf, ierr1 = _reactant_DRF(cosphi2, y, oneT)
    drj, ierr2 = _reactant_DRJ(cosphi2, y, oneT, p)

    return Base.ifelse(
        (ierr1 == 0) & (ierr2 == 0),
        sinphi * muladd(n * sinphi2 * inv3, drj, drf),
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

function _reactant_incomplete_Pi(n::A, phi::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n = T(n)
    phi = T(phi)
    m = T(m)
    oneT = one(T)
    ans = zero(T)

    Reactant.@trace if m < zero(T)
        mc = oneT - m
        imc = inv(mc)
        mN = -m * imc
        sinphi = sin(phi)
        sinphi2 = sinphi * sinphi
        phiN = asin(sqrt(mc / muladd(-m, sinphi2, oneT)) * sinphi)
        nN = (n - m) * imc

        ans =
            sqrt(imc) / nN *
            (mN * F(phiN, mN) + imc * n * _reactant_incomplete_Pi_nonneg(nN, phiN, mN))
    else
        ans = _reactant_incomplete_Pi_nonneg(n, phi, m)
    end

    return ans
end

CarlsonAlg.Pi(n::Reactant.TracedRNumber, phi::Real, m::Real) =
    _reactant_incomplete_Pi(n, phi, m)
CarlsonAlg.Pi(n::Real, phi::Reactant.TracedRNumber, m::Real) =
    _reactant_incomplete_Pi(n, phi, m)
CarlsonAlg.Pi(n::Real, phi::Real, m::Reactant.TracedRNumber) =
    _reactant_incomplete_Pi(n, phi, m)
CarlsonAlg.Pi(n::Reactant.TracedRNumber, phi::Reactant.TracedRNumber, m::Real) =
    _reactant_incomplete_Pi(n, phi, m)
CarlsonAlg.Pi(n::Reactant.TracedRNumber, phi::Real, m::Reactant.TracedRNumber) =
    _reactant_incomplete_Pi(n, phi, m)
CarlsonAlg.Pi(n::Real, phi::Reactant.TracedRNumber, m::Reactant.TracedRNumber) =
    _reactant_incomplete_Pi(n, phi, m)
CarlsonAlg.Pi(
    n::Reactant.TracedRNumber,
    phi::Reactant.TracedRNumber,
    m::Reactant.TracedRNumber,
) = _reactant_incomplete_Pi(n, phi, m)

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

    Reactant.@trace if n == zeroT
        ans = K(m)
    elseif (m == zeroT) | (m == oneT)
        ans = infT
    else
        ans = _reactant_cel(sqrt(oneT - m), oneT - n, oneT, oneT)
    end

    return ans
end

function _reactant_complete_Pi(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    n = T(n)
    m = T(m)
    oneT = one(T)
    nanT = T(NaN)
    ans = zero(T)

    Reactant.@trace if m > oneT
        ans = nanT
    elseif n > oneT
        ans = K(m) - _reactant_complete_Pi_core(m / n, m)
    else
        ans = _reactant_complete_Pi_core(n, m)
    end

    return ans
end

CarlsonAlg.Pi(n::Reactant.TracedRNumber, m::Real) = _reactant_complete_Pi(n, m)
CarlsonAlg.Pi(n::Real, m::Reactant.TracedRNumber) = _reactant_complete_Pi(n, m)
CarlsonAlg.Pi(n::Reactant.TracedRNumber, m::Reactant.TracedRNumber) =
    _reactant_complete_Pi(n, m)

end
