module JacobiEllipticReactantExt

using JacobiElliptic
using Reactant
import JacobiElliptic: CarlsonAlg

function _reactant_DRF(y::Reactant.AnyTracedRArray{Tc,N}) where {Tc,N}
    zeroT = Tc(0)
    oneT = Tc(1)
    inv3 = Tc(1 / 3)
    inv4 = Tc(1 / 4)
    ERRTOL = Tc((4 * eps(Tc) / 2)^Tc(1 / 6))
    C1 = Tc(1 / 24)
    C2 = Tc(3 / 44)
    C3 = Tc(1 / 14)
    C0 = Tc(1 / 10)

    XN = zeroT .* y
    YN = y
    ZN = zeroT .* y .+ oneT

    MU = (XN .+ YN .+ ZN) .* inv3
    ninvMU = .-inv.(MU)
    XNDEV = muladd.(ninvMU, MU .+ XN, 2)
    YNDEV = muladd.(ninvMU, MU .+ YN, 2)
    ZNDEV = muladd.(ninvMU, MU .+ ZN, 2)
    EPSLON = max.(abs.(XNDEV), max.(abs.(YNDEV), abs.(ZNDEV)))
    active = EPSLON .>= ERRTOL

    Reactant.@trace while sum(active) > 0
        XNROOT = sqrt.(XN)
        YNROOT = sqrt.(YN)
        ZNROOT = sqrt.(ZN)
        YNROOTZNROOT = YNROOT .* ZNROOT
        LAMDA = muladd.(XNROOT, YNROOT .+ ZNROOT, YNROOTZNROOT)

        next_XN = (XN .+ LAMDA) .* inv4
        next_YN = (YN .+ LAMDA) .* inv4
        next_ZN = (ZN .+ LAMDA) .* inv4

        XN = ifelse.(active, next_XN, XN)
        YN = ifelse.(active, next_YN, YN)
        ZN = ifelse.(active, next_ZN, ZN)

        MU = (XN .+ YN .+ ZN) .* inv3
        ninvMU = .-inv.(MU)
        XNDEV = muladd.(ninvMU, MU .+ XN, 2)
        YNDEV = muladd.(ninvMU, MU .+ YN, 2)
        ZNDEV = muladd.(ninvMU, MU .+ ZN, 2)
        EPSLON = max.(abs.(XNDEV), max.(abs.(YNDEV), abs.(ZNDEV)))
        active = EPSLON .>= ERRTOL
    end

    XNDEVYNDEV = XNDEV .* YNDEV
    E2 = muladd.(-ZNDEV, ZNDEV, XNDEVYNDEV)
    E3 = XNDEVYNDEV .* ZNDEV
    S = oneT .+ muladd.(E2, muladd.(-C2, E3, muladd.(C1, E2, -C0)), C3 .* E3)
    return S ./ sqrt.(MU)
end

function _reactant_DRF_whilebody(XN::A, YN::B, ZN::C, ERRTOL::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    inv3 = T(1 / 3)
    inv4 = T(1 / 4)

    MU = (XN + YN + ZN) * inv3
    ninvMU = -inv(MU)
    XNDEV = muladd(ninvMU, MU + XN, 2)
    YNDEV = muladd(ninvMU, MU + YN, 2)
    ZNDEV = muladd(ninvMU, MU + ZN, 2)
    EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))

    XNROOT = sqrt(XN)
    YNROOT = sqrt(YN)
    ZNROOT = sqrt(ZN)
    YNROOTZNROOT = YNROOT * ZNROOT
    LAMDA = muladd(XNROOT, YNROOT + ZNROOT, YNROOTZNROOT)
    next_XN = (XN + LAMDA) * inv4
    next_YN = (YN + LAMDA) * inv4
    next_ZN = (ZN + LAMDA) * inv4

    return (next_XN, next_YN, next_ZN, XNDEV, YNDEV, ZNDEV, MU, EPSLON >= ERRTOL)
end

function _reactant_DRF_ifbody(X::A, Y::B, Z::C, ERRTOL::D) where {A,B,C,D}
    T = promote_type(A, B, C, D)
    C1 = T(1 / 24)
    C2 = T(3 / 44)
    C3 = T(1 / 14)
    C0 = T(1 / 10)

    XN = X
    YN = Y
    ZN = Z
    MU = zero(T)
    XNDEV = zero(T)
    YNDEV = zero(T)
    ZNDEV = zero(T)
    active = true

    Reactant.@trace while active
        XN, YN, ZN, XNDEV, YNDEV, ZNDEV, MU, active = _reactant_DRF_whilebody(XN, YN, ZN, ERRTOL)
    end

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


#----------------------------------------------------------------------------------------
# Elliptic K(m)
#----------------------------------------------------------------------------------------

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
        ans = drf / sqrt(one_plus_m_abs)
    elseif m > oneT
        m_inv = inv(m)
        drf, _ = _reactant_DRF(zeroT, oneT - m_inv, oneT)
        ans = inv(sqrt(m)) * drf
    elseif m == oneT
        ans = infT
    else
        drf, _ = _reactant_DRF(zeroT, oneT - m, oneT)
        ans = drf
    end
    return ans
end

function CarlsonAlg.K(m::Reactant.AnyTracedRArray{Tc,N}) where {Tc,N}
    zeroT = Tc(0)
    oneT = Tc(1)
    nanT = Tc(NaN)
    infT = Tc(Inf)
    safe_one = zeroT .* m .+ oneT
    safe_half = zeroT .* m .+ Tc(0.5)

    nanmask = isnan.(m)
    negmask = m .< zeroT
    gtmask = m .> oneT
    eqmask = m .== oneT

    m_abs = .-m
    one_plus_m_abs = m_abs .+ oneT
    safe_m_transformed = ifelse.(negmask, m_abs ./ one_plus_m_abs, zeroT .* m)
    safe_m_inv = ifelse.(gtmask, inv.(m), safe_one)

    drf_input = ifelse.(
        negmask,
        oneT .- safe_m_transformed,
        ifelse.(gtmask, oneT .- safe_m_inv, ifelse.(nanmask .| eqmask, safe_half, oneT .- m)),
    )
    drf = _reactant_DRF(drf_input)

    neg_res = drf ./ sqrt.(one_plus_m_abs)
    gt_res = inv.(sqrt.(ifelse.(gtmask, m, safe_one))) .* drf

    return ifelse.(
        nanmask,
        zeroT .* m .+ nanT,
        ifelse.(eqmask, zeroT .* m .+ infT, ifelse.(negmask, neg_res, ifelse.(gtmask, gt_res, drf))),
    )
end

end
