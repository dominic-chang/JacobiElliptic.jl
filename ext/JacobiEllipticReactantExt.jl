module JacobiEllipticReactantExt

using JacobiElliptic
using Reactant
import JacobiElliptic: CarlsonAlg

@inline _reactant_sqrt(x::Reactant.TracedRNumber) = sqrt(x)

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
    oneT = one(T)
    C1 = T(Tc(1 / 24))
    C2 = T(Tc(3 / 44))
    C3 = T(Tc(1 / 14))
    C0 = T(Tc(1 / 10))

    ans = zeroT
    err = Base.ifelse(
        min(X, Y, Z) < zeroT,
        1,
        Base.ifelse(max(X, Y, Z) > UPLIM, 3, Base.ifelse(min(X + Y, X + Z, Y + Z) < LOLIM, 2, 0)),
    )

    XN = X
    YN = Y
    ZN = Z
    MU = zeroT
    XNDEV = zeroT
    YNDEV = zeroT
    ZNDEV = zeroT
    XNROOT = X
    YNROOT = Y
    ZNROOT = Z
    YNROOTZNROOT = X
    LAMDA = X

    flag = err == 0
    for _ in 1:40
        MU_iter = (XN + YN + ZN) / 3
        ninvMU = -inv(MU_iter)
        XNDEV_iter = muladd(ninvMU, (MU_iter + XN), 2)
        YNDEV_iter = muladd(ninvMU, (MU_iter + YN), 2)
        ZNDEV_iter = muladd(ninvMU, (MU_iter + ZN), 2)
        EPSLON = max(abs(XNDEV_iter), abs(YNDEV_iter), abs(ZNDEV_iter))
        next_flag = flag & !(EPSLON < ERRTOL)

        next_XNROOT = _reactant_sqrt(XN)
        next_YNROOT = _reactant_sqrt(YN)
        next_ZNROOT = _reactant_sqrt(ZN)
        next_YNROOTZNROOT = next_YNROOT * next_ZNROOT
        next_LAMDA = muladd(next_XNROOT, next_YNROOT + next_ZNROOT, next_YNROOTZNROOT)

        MU = Base.ifelse(flag, MU_iter, MU)
        XNDEV = Base.ifelse(flag, XNDEV_iter, XNDEV)
        YNDEV = Base.ifelse(flag, YNDEV_iter, YNDEV)
        ZNDEV = Base.ifelse(flag, ZNDEV_iter, ZNDEV)
        XNROOT = Base.ifelse(next_flag, next_XNROOT, XNROOT)
        YNROOT = Base.ifelse(next_flag, next_YNROOT, YNROOT)
        ZNROOT = Base.ifelse(next_flag, next_ZNROOT, ZNROOT)
        YNROOTZNROOT = Base.ifelse(next_flag, next_YNROOTZNROOT, YNROOTZNROOT)
        LAMDA = Base.ifelse(next_flag, next_LAMDA, LAMDA)
        XN = Base.ifelse(next_flag, (XN + next_LAMDA) / 4, XN)
        YN = Base.ifelse(next_flag, (YN + next_LAMDA) / 4, YN)
        ZN = Base.ifelse(next_flag, (ZN + next_LAMDA) / 4, ZN)
        flag = next_flag
    end

    XNDEVYNDEV = XNDEV * YNDEV
    E2 = muladd(-ZNDEV, ZNDEV, XNDEVYNDEV)
    E3 = XNDEVYNDEV * ZNDEV
    S = oneT + muladd(E2, muladd(-C2, E3, muladd(C1, E2, -C0)), C3 * E3)
    ans = S / _reactant_sqrt(MU)

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

function CarlsonAlg.K(m::Reactant.TracedRNumber)
    T = typeof(m)
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
        drf, _ = _reactant_DRF(zeroT, oneT - m_transformed, oneT)
        ans = drf / _reactant_sqrt(one_plus_m_abs)
    elseif m > oneT
        m_inv = inv(m)
        drf, _ = _reactant_DRF(zeroT, oneT - m_inv, oneT)
        ans = inv(_reactant_sqrt(m)) * drf
    elseif m == oneT
        ans = infT
    else
        drf, _ = _reactant_DRF(zeroT, oneT - m, oneT)
        ans = drf
    end
    return ans
end

end
