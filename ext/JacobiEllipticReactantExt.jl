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

function CarlsonAlg.F(
    phi::Reactant.AnyTracedRArray{Tp,N},
    m::Real,
) where {Tp,N}
    T = promote_type(Tp, typeof(m))
    phi = T.(phi)
    mT = T(m)
    oneT = one(T)
    halfpi = T(Base.π / 2)

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
    halfpi = T(Base.π / 2)

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
