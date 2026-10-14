module JacobiEllipticReactantExt

using JacobiElliptic
using Reactant
import JacobiElliptic: ArithmeticGeometricMeanAlg, CarlsonAlg

# Number of ARM iterations needed for conversion.
_reactant_agm_iters(::Type{<:Any}) = 8
_reactant_agm_iters(::Type{Float32}) = 5

# Number of Carlson iterations for DRF / DRD.
# The Carlson algorithm for these converge linearly instead of logarithmically as in AGM
# DRC / DRJ duplication passes. Scan over the (n, φ, m) domain Pi feeds them (30M
# triples): DRJ tops out at 12 passes for Float64 and 9 for Float32, the DRC nested inside
# its loop at 9 and 5. One spare each. DRJ x DRC unrolls to ~130 inner bodies, so the
# callers below are collapsed to a single DRJ evaluation to keep that affordable.
_reactant_drj_iters(::Type{<:Any}) = 13
_reactant_drj_iters(::Type{Float32}) = 10
_reactant_drc_iters(::Type{<:Any}) = 10
_reactant_drc_iters(::Type{Float32}) = 6

# Number of iterations needed for am. Also controls the buffer size
_reactant_am_iters(::Type{<:Any}) = 9
_reactant_am_iters(::Type{Float32}) = 7

_reactant_cel_iters(::Type{<:Any}) = 9
_reactant_cel_iters(::Type{Float32}) = 7

_reactant_duplication_iters(::Type{<:Any}) = 10
_reactant_duplication_iters(::Type{Float32}) = 6

@inline function _reactant_ellipke_base(m::T) where {T}
    a = one(T)
    b = √(one(T) - m)
    sum_c2 = m / 2
    pow2 = one(T)
    for _ = 1:_reactant_agm_iters(Reactant.unwrapped_eltype(T))
        c = (a - b) / 2
        sum_c2 = muladd(pow2, c * c, sum_c2)
        pow2 += pow2
        a, b = (a + b) / 2, √(a * b)
    end
    k = T(π / 2) / a
    e = k * (one(T) - sum_c2)
    return (k, e)
end

#----------------------------------------------------------------------------------------
# Complete AGM K(m), E(m)
#----------------------------------------------------------------------------------------

# Branchless counterpart of `src/complete_reduction.jl`, shared by the three complete
# integral overlays. The CPU ladder picks one of five cases; here every case is evaluated
# and then selected, so the kernel runs exactly once on a parameter always inside [0, 1).
#
# `_reactant_reduce_m` returns the reduced parameter together with the predicates that
# `_reactant_complete_lift` needs to map a kernel's `(K, E)` back to the original argument.
@inline function _reactant_reduce_m(m::T) where {T}
    oneT = one(T)
    isneg = m < zero(T)
    isbig = m > oneT
    reduced = Base.ifelse(isneg, -m / (oneT - m), Base.ifelse(isbig, inv(m), m))
    # The reduction does not always land strictly inside [0, 1): -m/(1 - m) rounds to
    # exactly 1 once m ≤ -9.1e15, and m == -Inf gives Inf/Inf == NaN. `src/AGM.jl` and
    # `src/Carlson.jl` screen both out before calling the kernel, so they are screened
    # here too, or the kernel runs with b == 0 and returns garbage.
    degenerate = (m == oneT) | (reduced == oneT)
    invalid = isnan(m) | isnan(reduced)
    # Clamping also keeps √(1 - reduced) away from 0, so the unselected branch of every
    # select below still has finite gradients.
    reduced = Base.ifelse(degenerate | invalid, zero(T), reduced)
    return (reduced, isneg, isbig, degenerate, invalid)
end

@inline function _reactant_complete_lift(m::T, reduction, k::T, e::T) where {T}
    oneT = one(T)
    nanT = T(NaN)
    reduced, isneg, isbig, degenerate, invalid = reduction
    k = Base.ifelse(degenerate, T(Inf), k)
    e = Base.ifelse(degenerate, oneT, e)
    sq = √(Base.ifelse(isneg, oneT - m, Base.ifelse(isbig, m, oneT)))
    K = Base.ifelse(isneg | isbig, k / sq, k)
    E = Base.ifelse(isneg, e * sq, Base.ifelse(isbig, sq * (e - (oneT - reduced) * k), e))
    return (Base.ifelse(invalid, nanT, K), Base.ifelse(invalid, nanT, E))
end

# The complete integrals are the one place where three overlays share an argument
# reduction, so it lives in `_reactant_reduce_m` / `_reactant_complete_lift` just above.
Reactant.@reactant_overlay function ArithmeticGeometricMeanAlg.ellipke(m::T) where {T}
    reduction = _reactant_reduce_m(m)
    k, e = _reactant_ellipke_base(reduction[1])
    return _reactant_complete_lift(m, reduction, k, e)
end

Reactant.@reactant_overlay ArithmeticGeometricMeanAlg.K(m::T) where {T} =
    ArithmeticGeometricMeanAlg.ellipke(m)[1]

Reactant.@reactant_overlay ArithmeticGeometricMeanAlg.E(m::T) where {T} =
    ArithmeticGeometricMeanAlg.ellipke(m)[2]

#----------------------------------------------------------------------------------------
# Elliptic F(ϕ, m)
#----------------------------------------------------------------------------------------

function _reactant_DRF_ifbody(X::A, Y::B, Z::C) where {A,B,C}
    T = promote_type(A, B, C)
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

    for _ = 1:_reactant_duplication_iters(Reactant.unwrapped_eltype(T))
        XNROOT = √(XN)
        YNROOT = √(YN)
        ZNROOT = √(ZN)
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
    return S / √(MU)
end

Reactant.@reactant_overlay function CarlsonAlg.DRF(X::A, Y::B, Z::C) where {A,B,C}
    T = promote_type(A, B, C)
    Tc = Reactant.unwrapped_eltype(T)
    X = T(X)
    Y = T(Y)
    Z = T(Z)
    zeroT = zero(T)

    LOLIM = T(5 * floatmin(Tc))
    UPLIM = T(floatmax(Tc) / 5)

    ans = zeroT
    minxyz = min(X, Y, Z) < zeroT
    maxxyz = max(X, Y, Z) > UPLIM
    mxyxzyz = min(X + Y, X + Z, Y + Z) < LOLIM
    err = Base.ifelse(minxyz, 1, Base.ifelse(maxxyz, 3, Base.ifelse(mxyxzyz, 2, 0)))

    # Safe to evaluate unconditionally now that the duplication loop has a fixed trip
    # count: bad arguments give NaN rather than spinning, and the ifelse discards them.
    ans = Base.ifelse(err == 0, _reactant_DRF_ifbody(X, Y, Z), zeroT)

    return (ans, err)
end

function _reactant_rawF(sinφ::A, m::B) where {A,B}
    T = promote_type(A, B)
    sinφ = T(sinφ)
    m = T(m)
    oneT = one(T)
    infT = T(Inf)
    sinφ2 = sinφ * sinφ
    drf, _ = CarlsonAlg.DRF(oneT - sinφ2, muladd(-m, sinφ2, oneT), oneT)
    return Base.ifelse((abs(sinφ) == oneT) & (m == oneT), sign(sinφ) * infT, sinφ * drf)
end

# Branch-free: both cases feed a single rawF. The Abramowitz & Stegun (17.4.3) reduction
# only picks that call's argument and adds the 2*fld(.)*K(m) offset, so selecting with
# ifelse costs one K but saves the second rawF -- and rawF is the expensive half.
function _reactant_internal_F(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    φ = T(φ)
    m = T(m)
    halfπ = T(π / 2)
    piT = T(π)

    wrapped = abs(φ) > halfπ
    φ2 = φ + halfπ
    sinarg = Base.ifelse(wrapped, cos(mod(φ2, piT)), sin(φ))
    raw = _reactant_rawF(sinarg, m)
    offset = 2 * fld(φ2, piT) * CarlsonAlg.K(m)
    ans = Base.ifelse(wrapped, offset - raw, raw)
    return Base.ifelse(isnan(φ) | isnan(m), T(NaN), ans)
end

Reactant.@reactant_overlay function CarlsonAlg.F(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    φ = T(φ)
    m = T(m)
    oneT = one(T)
    halfπ = T(Base.π / 2)
    zeroT = zero(T)

    # A&S 17.4.15 (m > 1) and 17.4.17 (m < 0) both reduce to one evaluation at |angle|
    # with a parameter in [0, 1], so select the pair and call _reactant_internal_F once
    # instead of once per branch. Unselected branches produce NaN (√ of a negative,
    # m1m = -Inf at m == 1); every one of them is discarded by an ifelse.
    isbig = m > oneT
    isneg = m < zeroT

    m12big = √(m)
    theta = asin(m12big * sin(φ))
    n = -m
    one_plus_n = oneT + n
    m12neg = inv(√(one_plus_n))
    m1m = n / one_plus_n
    newφ = halfπ - φ

    angle = Base.ifelse(isbig, abs(theta), Base.ifelse(isneg, abs(newφ), abs(φ)))
    param = Base.ifelse(isbig, inv(m), Base.ifelse(isneg, m1m, m))
    f = _reactant_internal_F(angle, param)

    scale = Base.ifelse(
        isbig,
        sign(theta) / m12big,
        Base.ifelse(isneg, -sign(newφ) * m12neg, sign(φ)),
    )
    sf = scale * f
    return Base.ifelse(isneg, m12neg * CarlsonAlg.K(m1m) + sf, sf)
end

function _reactant_am(u::A, m::B, tol::C) where {A,B,C}

    T = promote_type(A, B, C)

    sqrt_tol = √(tol)
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

    # size of buffer
    NP = _reactant_am_iters(Reactant.unwrapped_eltype(T))
    _ambuf = Vector{T}(undef, NP)
    a = one(T)
    b = √(m1)
    c = √(m)
    for i = 1:NP
        a, b, c = ((a + b) / 2, √(a * b), (a - b) / 2)
        _ambuf[i] = c / a
    end

    φ = a * u * (2^NP)
    for i = NP:-1:1
        φ = (φ + asin(_ambuf[i] * sin(φ))) / 2
    end

    return Base.ifelse(flag, Base.ifelse(abs(c) > tol, T(NaN), φ), ans)
end

function _reactant_am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    return _reactant_am(u, m, eps(T))
end

"""
    am(u::Real, m::Real)

Returns amplitude, φ, such that u = F(φ | m)
"""
Reactant.@reactant_overlay function CarlsonAlg.am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    piT = T(Base.π)
    halfpiT = T(Base.π / 2)

    oneT = one(T)
    isneg = m < zero(T)
    isbig = m > oneT

    mu1 = inv(oneT - m)
    mu = -m * mu1
    sqrtmu1 = √(mu1)
    k = √(m)
    minv = inv(m)

    φ = _reactant_am(
        Base.ifelse(isneg, u / sqrtmu1, Base.ifelse(isbig, k * u, u)),
        Base.ifelse(isneg, mu, Base.ifelse(isbig, minv, m)),
    )
    s = sin(φ)
    t = floor((φ + halfpiT) / piT)
    negans = t * piT + cospi(t) * asin(sqrtmu1 * s / √(oneT - mu * s^2))

    return Base.ifelse(isneg, negans, Base.ifelse(isbig, asin(s / k), φ))
end

for (f, a, b, c) in (
    (:sn, :(sin(φ)), :(sqrtmu1 * s), :(√(mu) * sin(φ))),
    (:cn, :(cos(φ)), :(cos(φ)), :(√(one(T) - mu * sin(φ)^2))),
    (:dn, :(√(one(T) - m * sin(φ)^2)), :(one(T)), :(cos(φ))),
)
    @eval begin
        Reactant.@reactant_overlay function CarlsonAlg.$f(u::A, m::B) where {A,B}
            T = promote_type(A, B)
            oneT = one(T)
            isneg = m < zero(T)
            isbig = m > oneT

            mu1 = inv(oneT - m)
            muneg = -m * mu1
            sqrtmu1 = √(mu1)
            mubig = inv(m)

            φ = _reactant_am(
                Base.ifelse(isneg, u / sqrtmu1, Base.ifelse(isbig, u * √(m), u)),
                Base.ifelse(isneg, muneg, Base.ifelse(isbig, mubig, m)),
            )
            s = sin(φ)

            negval = let mu = muneg
                ($b) / √(oneT - mu * s^2)
            end
            bigval = let mu = mubig
                $c
            end
            return Base.ifelse(isneg, negval, Base.ifelse(isbig, bigval, $a))
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

    for _ = 1:_reactant_duplication_iters(Reactant.unwrapped_eltype(T))
        XNROOT = √(XN)
        YNROOT = √(YN)
        ZNROOT = √(ZN)
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
    return threeT * SIGMA + POWER4 * (oneT + S1 + S2) / (MU * √(MU))
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
Reactant.@reactant_overlay function CarlsonAlg.K(m::T) where {T}
    reduction = _reactant_reduce_m(m)
    drf, _ = CarlsonAlg.DRF(zero(T), one(T) - reduction[1], one(T))
    # Only the first component is wanted; the placeholder second one and the `E` half of
    # the lift are dead and get folded away.
    return _reactant_complete_lift(m, reduction, drf, zero(T))[1]
end

Reactant.@reactant_overlay function CarlsonAlg.ellipke(m::T) where {T}
    reduction = _reactant_reduce_m(m)
    y = one(T) - reduction[1]
    drf, _ = CarlsonAlg.DRF(zero(T), y, one(T))
    drd, _ = CarlsonAlg.DRD(zero(T), y, one(T))
    e = muladd(-reduction[1] * T(1 / 3), drd, drf)
    return _reactant_complete_lift(m, reduction, drf, e)
end

function _reactant_rawE(sinφ::A, m::B) where {A,B}
    T = promote_type(A, B)
    sinφ = T(sinφ)
    m = T(m)
    oneT = one(T)
    inv3 = T(1 / 3)
    nanT = T(NaN)

    isbig = m > oneT

    m12 = √(m)
    sinalpha = m12 * sinφ
    sinalpha2 = sinalpha * sinalpha
    cosalpha2 = oneT - sinalpha2
    yalpha = muladd(-inv(m), sinalpha2, oneT)

    sinφ2 = sinφ * sinφ
    cosφ2 = oneT - sinφ2
    y = muladd(-m, sinφ2, oneT)

    X = Base.ifelse(isbig, cosalpha2, cosφ2)
    Y = Base.ifelse(isbig, yalpha, y)
    drf, ierr1 = CarlsonAlg.DRF(X, Y, oneT)
    drd, ierr2 = CarlsonAlg.DRD(X, Y, oneT)

    val = Base.ifelse(
        isbig,
        sinalpha * (drf - sinalpha2 * inv3 * drd) / m12,
        sinφ * (drf - m * sinφ2 * inv3 * drd),
    )
    return Base.ifelse(
        (ierr1 == 0) & (ierr2 == 0),
        val,
        Base.ifelse((ierr1 == 2) & (ierr2 == 2), sinφ, nanT),
    )
end

Reactant.@reactant_overlay function CarlsonAlg.E(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    φ = T(φ)
    m = T(m)
    halfπ = T(Base.π / 2)
    piT = T(Base.π)
    nanT = T(NaN)

    wrapped = 2 * abs(φ) > piT
    φ2 = φ + halfπ
    sinarg = Base.ifelse(wrapped, cos(mod(φ2, piT)), sin(φ))
    raw = _reactant_rawE(sinarg, m)
    offset = 2 * fld(φ2, piT) * CarlsonAlg.E(m)
    ans = Base.ifelse(wrapped, offset - raw, raw)
    return Base.ifelse(isnan(φ) | isnan(m), nanT, ans)
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

    XN = X
    YN = Y
    MU = zero(T)
    SN = zero(T)

    MU = (XN + YN + YN) * inv3
    invMU = inv(MU)
    SN = muladd(invMU, YN + MU, -twoT)

    for _ = 1:_reactant_drc_iters(Reactant.unwrapped_eltype(T))
        LAMDA = muladd(twoT * √(XN), √(YN), YN)
        XN = (XN + LAMDA) * inv4
        YN = (YN + LAMDA) * inv4

        MU = (XN + YN + YN) * inv3
        invMU = inv(MU)
        SN = muladd(invMU, YN + MU, -twoT)
    end

    S = SN^2 * muladd(SN, muladd(SN, muladd(SN, C2, C3), C1), C4)
    return (oneT + S) / √(MU)
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

    ans = Base.ifelse(err == 0, _reactant_DRC_ifbody(X, Y), zeroT)

    return (ans, err)
end
function _reactant_DRJ_ifbody(X::A, Y::B, Z::C, P::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    oneT = one(T)
    twoT = T(2)
    threeT = T(3)
    inv4 = T(1 / 4)
    inv5 = T(1 / 5)
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

    for _ = 1:_reactant_drj_iters(Reactant.unwrapped_eltype(T))
        XNROOT = √(XN)
        YNROOT = √(YN)
        ZNROOT = √(ZN)
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
    ans = threeT * SIGMA + POWER4 * (S1 + S2 + S3) / (MU * √(MU))

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

    ok = err == 0
    bodyans, bodyerr = _reactant_DRJ_ifbody(X, Y, Z, P)
    ans = Base.ifelse(ok, bodyans, zeroT)
    err = Base.ifelse(ok, bodyerr, err)

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
        atan(t * √(h)) / √(h),
        Base.ifelse(
            h == zero(T),
            t,
            begin
                sqrt_neg_h = √(-h)
                arg = t * sqrt_neg_h
                Base.ifelse(abs(arg) < oneT, atanh(arg), _reactant_custom_atanh(arg)) /
                sqrt_neg_h
            end,
        ),
    )
end

function _reactant_cel(kc::A, p::B, a::C, b::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    oneT = one(T)
    twoT = T(2)
    zeroT = zero(T)
    π_over_2 = T(π / 2)

    kc = abs(kc)
    p = T(p)
    a = T(a)
    b = T(b)

    e = copy(kc)
    m = one(T)
    f = zero(T)
    g = zero(T)
    q = zero(T)

    ispos = p > zeroT
    p_pos = √(p)
    b_pos = b / p_pos
    f0 = kc * kc
    g0 = oneT - p
    q0 = (b - a * p) * (oneT - f0)
    p_neg = √((f0 - p) / g0)
    a_neg = (a - b) / g0
    b_neg = -q0 * (g0^2 * p_neg) + a_neg * p_neg

    p = Base.ifelse(ispos, p_pos, p_neg)
    b = Base.ifelse(ispos, b_pos, b_neg)
    a = Base.ifelse(ispos, a, a_neg)

    f = copy(a)
    invp = inv(p)
    a = muladd(invp, b, a)
    g = e * invp
    b = twoT * muladd(f, g, b)
    p = g + p
    g = copy(m)
    m = kc + m

    for _ = 1:_reactant_cel_iters(Reactant.unwrapped_eltype(T))
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
    end

    return π_over_2 * muladd(a, m, b) / (m * (m + p))
end

function _reactant_rawPi(n::A, sinφ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n = T(n)
    sinφ = T(sinφ)
    m = T(m)
    oneT = one(T)
    inv3 = T(1 / 3)
    infT = T(Inf)
    nanT = T(NaN)

    sinφ2 = sinφ * sinφ
    cosφ2 = oneT - sinφ2
    y = muladd(-m, sinφ2, oneT)
    p = muladd(-n, sinφ2, oneT)
    drf, ierr1 = CarlsonAlg.DRF(cosφ2, y, oneT)
    drj, ierr2 = _reactant_DRJ(cosφ2, y, oneT, p)

    return Base.ifelse(
        (ierr1 == 0) & (ierr2 == 0),
        sinφ * (drf + n * sinφ2 * inv3 * drj),
        Base.ifelse(
            ((ierr1 == 2) & (ierr2 == 2)) | ((ierr1 == 0) & (ierr2 == 2)),
            infT,
            nanT,
        ),
    )
end

function _reactant_incomplete_Pi_core(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n = T(n)
    φ = T(φ)
    m = T(m)
    piT = T(π)
    halfπ = T(π / 2)
    nanT = T(NaN)
    ans = zero(T)

    # One rawPi -- and therefore one DRJ -- serves both cases.
    wrapped = 2 * abs(φ) > piT
    φ2 = φ + halfπ
    sinarg = Base.ifelse(wrapped, cos(mod(φ2, piT)), sin(φ))
    raw = _reactant_rawPi(n, sinarg, m)
    offset = 2 * fld(φ2, piT) * Pi(n, m)
    ans = Base.ifelse(wrapped, offset - raw, raw)
    return Base.ifelse(isnan(n) | isnan(φ) | isnan(m), nanT, ans)
end

function _reactant_incomplete_Pi_nonneg(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n = T(n)
    φ = T(φ)
    m = T(m)
    oneT = one(T)
    ans = zero(T)

    isbig = n > oneT
    sinφ, cosφ = sincos(φ)
    sinφ2 = sinφ * sinφ
    t1 = sinφ / (cosφ * √(muladd(-m, sinφ2, oneT)))
    h1 = (oneT - n) * (n - m) / n
    core = _reactant_incomplete_Pi_core(Base.ifelse(isbig, m / n, n), φ, m)
    return Base.ifelse(isbig, _reactant_FukushimaT(t1, h1) - core + F(φ, m), core)
end

Reactant.@reactant_overlay function CarlsonAlg.Pi(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n = T(n)
    φ = T(φ)
    m = T(m)
    oneT = one(T)

    isneg = m < zero(T)
    isnm = isneg & (n == m)
    piT = T(Base.π)
    halfπ = T(Base.π / 2)

    # n == m makes the imaginary modulus transformation 0/0, so use the closed form
    # Pi(m; φ | m) = (E(φ | m) - m*sin*cos/√(1 - m*sin^2)) / (1 - m).
    sc, cc = sincos(φ)
    closed = (CarlsonAlg.E(φ, m) - m * sc * cc / √(muladd(-m, sc * sc, oneT))) / (oneT - m)

    # Pi(n; φ + k*π | m) = Pi(n; φ | m) + 2k Pi(n | m); reduce before the transform,
    # whose asin cannot represent |φ| > π/2.
    φ2 = φ + halfπ
    wrapped = 2 * abs(φ) > piT
    φ_r = Base.ifelse(wrapped, mod(φ2, piT) - halfπ, φ)
    offset = Base.ifelse(wrapped, 2 * fld(φ2, piT) * CarlsonAlg.Pi(n, m), zero(T))

    mc = oneT - m
    imc = inv(mc)
    mN = -m * imc
    sinφ = sin(φ_r)
    sinφ2 = sinφ * sinφ
    φN = asin(√(mc / muladd(-m, sinφ2, oneT)) * sinφ)
    nN = (n - m) * imc

    nonneg = _reactant_incomplete_Pi_nonneg(
        Base.ifelse(isneg, nN, n),
        Base.ifelse(isneg, φN, φ),
        Base.ifelse(isneg, mN, m),
    )
    transformed = offset + √(imc) / nN * (mN * F(φN, mN) + imc * n * nonneg)

    return Base.ifelse(isnm, closed, Base.ifelse(isneg, transformed, nonneg))
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

    celv = _reactant_cel(√(oneT - m), oneT - n, oneT, oneT)
    return Base.ifelse(
        (n == oneT) | (m == oneT),
        infT,
        Base.ifelse(n == zeroT, CarlsonAlg.K(m), celv),
    )
end

Reactant.@reactant_overlay function CarlsonAlg.Pi(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    n = T(n)
    m = T(m)
    oneT = one(T)
    nanT = T(NaN)

    isbign = n > oneT
    core = _reactant_complete_Pi_core(Base.ifelse(isbign, m / n, n), m)
    ans = Base.ifelse(isbign, CarlsonAlg.K(m) - core, core)
    return Base.ifelse(m > oneT, nanT, ans)
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
