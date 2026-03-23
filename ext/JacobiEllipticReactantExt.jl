module JacobiEllipticReactantExt

using JacobiElliptic
using Reactant
import JacobiElliptic: CarlsonAlg, ArithmeticGeometricMeanAlg

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

function _reactant_DRF(
    X::Reactant.AnyTracedRArray{Tx,N},
    Y::Reactant.AnyTracedRArray{Ty,N},
    Z::Union{Real,Reactant.TracedRNumber},
) where {Tx,Ty,N}
    T = promote_type(Tx, Ty, typeof(Z))
    Tc = Reactant.unwrapped_eltype(T)
    zeroT = T(0)
    inv3 = T(1 / 3)
    inv4 = T(1 / 4)
    ERRTOL = T((4 * eps(Tc) / 2)^Tc(1 / 6))
    C1 = T(Tc(1 / 24))
    C2 = T(Tc(3 / 44))
    C3 = T(Tc(1 / 14))
    C0 = T(Tc(1 / 10))

    XN = T.(X)
    YN = T.(Y)
    ZN = zeroT .* XN .+ T(Z)

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
    S = one(T) .+ muladd.(E2, muladd.(-C2, E3, muladd.(C1, E2, -C0)), C3 .* E3)
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

function _reactant_rawF(
    sinphi::Reactant.AnyTracedRArray{Tp,N},
    m::Real,
) where {Tp,N}
    T = promote_type(Tp, typeof(m))
    sinphi = T.(sinphi)
    mT = T(m)
    oneT = one(T)
    infT = T(Inf)
    sinphi2 = sinphi .* sinphi
    drf = _reactant_DRF(oneT .- sinphi2, muladd.(-mT, sinphi2, oneT), oneT)
    infmask = (abs.(sinphi) .== oneT) .& (mT == oneT)
    return ifelse.(infmask, sign.(sinphi) .* infT, sinphi .* drf)
end

function _reactant_rawF(
    sinphi::Reactant.AnyTracedRArray{Tp,N},
    m::Reactant.TracedRNumber,
) where {Tp,N}
    T = promote_type(Tp, typeof(m))
    sinphi = T.(sinphi)
    mT = T(m)
    oneT = one(T)
    infT = T(Inf)
    sinphi2 = sinphi .* sinphi
    drf = _reactant_DRF(oneT .- sinphi2, muladd.(-mT, sinphi2, oneT), oneT)
    infmask = (abs.(sinphi) .== oneT) .& (mT == oneT)
    return ifelse.(infmask, sign.(sinphi) .* infT, sinphi .* drf)
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

function _reactant_internal_F(
    phi::Reactant.AnyTracedRArray{Tp,N},
    m::Real,
) where {Tp,N}
    T = promote_type(Tp, typeof(m))
    phi = T.(phi)
    mT = T(m)
    halfpi = T(π / 2)
    piT = T(π)
    nanT = T(NaN)
    zeroT = T(0)

    nanmask = isnan.(phi) .| isnan(mT)
    large_mask = abs.(phi) .> halfpi

    phi2 = phi .+ halfpi
    pi_arr = zeroT .* phi .+ piT
    reduced = 2 .* fld.(phi2, pi_arr) .* K(mT) .- _reactant_rawF(cos.(mod.(phi2, piT)), mT)
    standard = _reactant_rawF(sin.(phi), mT)
    ans = ifelse.(large_mask, reduced, standard)
    return ifelse.(nanmask, zeroT .* phi .+ nanT, ans)
end

function _reactant_F(phi::A, m::B) where {A,B}
    T = promote_type(A, B)
    phi = T(phi)
    m = T(m)
    oneT = one(T)
    halfpi = T(π / 2)
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

function CarlsonAlg.F(
    phi::Reactant.AnyTracedRArray{Tp,N},
    m::Real,
) where {Tp,N}
    T = promote_type(Tp, typeof(m))
    phi = T.(phi)
    mT = T(m)
    oneT = one(T)
    halfpi = T(π / 2)

    if mT > oneT
        m12 = sqrt(mT)
        theta = asin.(m12 .* sin.(phi))
        return sign.(theta) ./ m12 .* _reactant_internal_F(abs.(theta), inv(mT))
    elseif mT < zero(T)
        n = -mT
        one_plus_n = oneT + n
        m12 = inv(sqrt(one_plus_n))
        m1m = n / one_plus_n
        newphi = halfpi .- phi
        return m12 .* K(m1m) .- sign.(newphi) .* m12 .* _reactant_internal_F(abs.(newphi), m1m)
    else
        return sign.(phi) .* _reactant_internal_F(abs.(phi), mT)
    end
end

function CarlsonAlg.F(
    phi::Reactant.AnyTracedRArray{Tp,N},
    m::Reactant.TracedRNumber,
) where {Tp,N}
    T = promote_type(Tp, typeof(m))
    phi = T.(phi)
    mT = T(m)
    oneT = one(T)
    zeroT = zero(T)
    halfpi = T(π / 2)

    gtmask = mT > oneT
    negmask = mT < zeroT

    m12_gt = sqrt(ifelse(gtmask, mT, oneT))
    theta = asin.(m12_gt .* sin.(phi))
    gt_res = sign.(theta) ./ m12_gt .* _reactant_internal_F(abs.(theta), inv(ifelse(gtmask, mT, oneT)))

    n = -mT
    one_plus_n = oneT + n
    m12_neg = inv(sqrt(ifelse(negmask, one_plus_n, oneT)))
    m1m = ifelse(negmask, n / one_plus_n, zeroT)
    newphi = halfpi .- phi
    neg_res = m12_neg .* K(m1m) .- sign.(newphi) .* m12_neg .* _reactant_internal_F(abs.(newphi), m1m)

    std_res = sign.(phi) .* _reactant_internal_F(abs.(phi), mT)
    return ifelse.(gtmask, gt_res, ifelse.(negmask, neg_res, std_res))
end

CarlsonAlg.F(φ::Reactant.TracedRNumber, m::Real) = _reactant_F(φ, m)
CarlsonAlg.F(φ::Real, m::Reactant.TracedRNumber) = _reactant_F(φ, m)
CarlsonAlg.F(φ::Reactant.TracedRNumber, m::Reactant.TracedRNumber) = _reactant_F(φ, m)

end
