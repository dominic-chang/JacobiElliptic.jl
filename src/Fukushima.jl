module FukushimaAlg
using StaticArrays

include(joinpath(@__DIR__, "kernels.jl"))

export K, E, F, Pi, J
export sn, cn, dn, sc, sd, asn, acn

const ONE_DIV_PI = inv(π)

############################################################################
# Complete Elliptic Integrals (https://doi.org/10.1007/s10569-009-9228-z)
############################################################################
@inline _K1(x::T) where {T} = evalpoly(
    x,
    (
        T(1.5910034537907922),
        T(0.41600074399178694),
        T(0.24579151426410342),
        T(0.17948148291490615),
        T(0.14455605708755515),
        T(0.12320099331242772),
        T(0.10893881157429353),
        T(0.09885340987159291),
        T(0.09143962920174975),
        T(0.0858425915954139),
        T(0.08154111871830322),
    ),
)
@inline _K2(x::T) where {T} = evalpoly(
    x,
    (
        T(1.63525673226458),
        T(0.4711906261487323),
        T(0.3097284108314996),
        T(0.2522083117731357),
        T(0.22672562321968465),
        T(0.21577444672958598),
        T(0.21310877187734892),
        T(0.21602912460518828),
        T(0.2232558316330579),
        T(0.23418050129420992),
        T(0.24855768297226408),
        T(0.26636380989261754),
    ),
)
@inline _K3(x::T) where {T} = evalpoly(
    x,
    (
        T(1.685750354812596),
        T(0.5417318486132803),
        T(0.40152443839069024),
        T(0.3696424734208891),
        T(0.37606071535458363),
        T(0.4052358870851259),
        T(0.45329438175399905),
        T(0.5205189476511842),
        T(0.609426039204995),
        T(0.7242635222829089),
        T(0.8710138477098124),
        T(1.057652872753547),
    ),
)
@inline _K4(x::T) where {T} = evalpoly(
    x,
    (
        T(1.7443505972256133),
        T(0.6348642753719353),
        T(0.5398425641644455),
        T(0.5718927051937874),
        T(0.6702951362654062),
        T(0.8325865900109772),
        T(1.0738574482479333),
        T(1.4220914606754977),
        T(1.9203871834023047),
        T(2.6325525483316543),
        T(3.6521097473190394),
        T(5.115867135558866),
        T(7.224080007363877),
    ),
)
@inline _K5(x::T) where {T} = evalpoly(
    x,
    (
        T(1.8138839368169826),
        T(0.7631632457005573),
        T(0.7619286053215958),
        T(0.9510746536684279),
        T(1.315180671703161),
        T(1.9285606934774109),
        T(2.9375093425313787),
        T(4.594894405442878),
        T(7.33007122188172),
        T(11.871512597425301),
        T(19.45851374822938),
        T(32.20638657246427),
        T(53.73749198700555),
        T(90.27388602941),
    ),
)
@inline _K6(x::T) where {T} = evalpoly(
    x,
    (
        T(1.8989249102715535),
        T(0.9505217946182445),
        T(1.1510775899590158),
        T(1.7502391069863006),
        T(2.952676812636875),
        T(5.285800396121451),
        T(9.83248571665998),
        T(18.787148683275596),
        T(36.61468615273698),
        T(72.45292395127771),
        T(145.1079577347069),
        T(293.4786396308497),
        T(598.385181505501),
        T(1228.4200130758634),
        T(2536.5297553827645),
    ),
)
@inline _K7(x::T) where {T} = evalpoly(
    x,
    (
        T(2.0075983984243764),
        T(1.2484572312123474),
        T(1.9262346570764797),
        T(3.7512896400875877),
        T(8.119944554932045),
        T(18.665721308735552),
        T(44.603924842914374),
        T(109.50920543094983),
        T(274.2779548232414),
        T(697.5598008606327),
        T(1795.7160145002472),
        T(4668.38171679039),
        T(12235.762468136643),
        T(32290.17809718321),
        T(85713.07608195965),
        T(228672.1890493117),
        T(612757.2711915852),
    ),
)
@inline _K8(x::T) where {T} = evalpoly(
    x,
    (
        T(2.1565156474996434),
        T(1.7918056418494632),
        T(3.8267512874657132),
        T(10.386724683637972),
        T(31.403314054680703),
        T(100.92370394986955),
        T(337.3268282632273),
        T(1158.7079305678278),
        T(4060.9907421936323),
        T(14454.001840343448),
        T(52076.661075994045),
        T(189493.65914621568),
        T(695184.5762413896),
        T(2.567994048255285e6),
        T(9.541921966748387e6),
        T(3.563492744218076e7),
        T(1.3366929846120408e8),
        T(5.0335218668662846e8),
        T(1.90197572953866e9),
        T(7.208915015330104e9),
    ),
)
@inline _K9(x::T) where {T} = evalpoly(
    x,
    (
        T(2.3181226217125106),
        T(2.6169201502912327),
        T(7.897935075731356),
        T(30.502397154466724),
        T(131.48693655235286),
        T(602.9847637356492),
        T(2877.024617809973),
        T(14110.519919151804),
        T(70621.4408815654),
        T(358977.266582531),
        T(1.8472382637239718e6),
        T(9.600515416049214e6),
        T(5.030767708502367e7),
        T(2.654441886527128e8),
        T(1.4088623250287027e9),
        T(7.515687935373775e9),
    ),
)
@inline _K10(x::T) where {T} = evalpoly(
    x,
    (
        T(2.473596173751344),
        T(3.727624244118099),
        T(15.607393035549306),
        T(84.12850842805888),
        T(506.98181970406137),
        T(3252.2770581451236),
        T(21713.242419574344),
        T(149037.04518909327),
        T(1.0439993310899908e6),
        T(7.427974817042039e6),
        T(5.350383967558661e7),
        T(3.8924988699487084e8),
        T(2.8552883511008105e9),
        T(2.109007703876684e10),
        T(1.566998339477902e11),
        T(1.17022224242244e12),
        T(8.7779483236689375e12),
        T(6.610124275248495e13),
        T(4.994880537133888e14),
        T(3.7859743397240296e16),
    ),
)
@inline _Kq(x::T) where {T} = evalpoly(
    x,
    (
        T(0.0),
        T(0.0625),
        T(0.03125),
        T(0.0205078125),
        T(0.01513671875),
        T(0.011934280395507812),
        T(0.009816169738769531),
        T(0.008315593004226685),
        T(0.007199153304100037),
        T(0.00633745662344154),
        T(0.00565311038371874),
        T(0.005097046040418718),
        T(0.004636680381850056),
        T(0.004249547423822886),
        T(0.003919665602267974),
    ),
)

@inline _E1(x::T) where {T} = evalpoly(
    x,
    (
        T(1.5509733517804722),
        T(-0.4003010201031985),
        T(-0.07849861944294194),
        T(-0.034318853117591995),
        T(-0.0197180433173655),
        T(-0.01305950773199331),
        T(-0.009442372874146548),
        T(-0.007246728512402157),
        T(-0.00580742401295609),
        T(-0.004809187786009338),
    ),
)
@inline _E2(x::T) where {T} = evalpoly(
    x,
    (
        T(1.5101218320928198),
        T(-0.41711633390586755),
        T(-0.09012382040477457),
        T(-0.04372994401908431),
        T(-0.027965493064761784),
        T(-0.020644781177568104),
        T(-0.016650786739707237),
        T(-0.01426196082884252),
        T(-0.012759847429264804),
        T(-0.011799303775587354),
        T(-0.011197445703074968),
    ),
)
@inline _E3(x::T) where {T} = evalpoly(
    x,
    (
        T(1.4674622093394272),
        T(-0.43657629094633776),
        T(-0.10515555766694255),
        T(-0.05737184359324173),
        T(-0.04139162772734022),
        T(-0.03452772850528084),
        T(-0.031495443512532785),
        T(-0.030527000890325277),
        T(-0.0309169840192389),
        T(-0.03237139531475812),
        T(-0.03478996038640416),
    ),
)
@inline _E4(x::T) where {T} = evalpoly(
    x,
    (
        T(1.4226911334908792),
        T(-0.4595135196210487),
        T(-0.12525053982206188),
        T(-0.07813854509440948),
        T(-0.06471427847205),
        T(-0.06208433913173031),
        T(-0.06519703281557247),
        T(-0.07279389536257878),
        T(-0.084959075171781),
        T(-0.102539850131046),
        T(-0.12705358515769605),
        T(-0.1607911206912746),
    ),
)
@inline _E5(x::T) where {T} = evalpoly(
    x,
    (
        T(1.3754019718711163),
        T(-0.4872021832731848),
        T(-0.15331170134854022),
        T(-0.11184944491702783),
        T(-0.10884095252313576),
        T(-0.12295422312026907),
        T(-0.15221716396203505),
        T(-0.20049532364269734),
        T(-0.27617433306775174),
        T(-0.39351311430437586),
        T(-0.5757544060278792),
        T(-0.8605232357272398),
        T(-1.3088332057585401),
    ),
)
@inline _E6(x::T) where {T} = evalpoly(
    x,
    (
        T(1.3250244979582302),
        T(-0.5217276475575667),
        T(-0.19490643048212622),
        T(-0.17162372682201127),
        T(-0.20275465292641914),
        T(-0.27879895311853475),
        T(-0.42069845728100574),
        T(-0.675948400853106),
        T(-1.1363431218392293),
        T(-1.9767211439543984),
        T(-3.5316967730957227),
        T(-6.446753640156048),
        T(-11.97703130208884),
    ),
)
@inline _E7(x::T) where {T} = evalpoly(
    x,
    (
        T(1.2707074796501499),
        T(-0.5668391682878666),
        T(-0.2621607934324926),
        T(-0.2922441735330774),
        T(-0.4403978408504232),
        T(-0.7749476413813975),
        T(-1.498870837987561),
        T(-3.089708310445187),
        T(-6.6675959033810015),
        T(-14.89436036517319),
        T(-34.18120574251449),
        T(-80.15895841905397),
        T(-191.34894807629848),
        T(-463.5938853480342),
        T(-1137.38082216936),
    ),
)
@inline _E8(x::T) where {T} = evalpoly(
    x,
    (
        T(1.2110560275684594),
        T(-0.6303064132874558),
        T(-0.38716640952066916),
        T(-0.5922782353119346),
        T(-1.23755558451305),
        T(-3.0320566617452474),
        T(-8.18168822157359),
        T(-23.55507217389693),
        T(-71.04099935893065),
        T(-221.879685319235),
        T(-712.1364793277636),
        T(-2336.1253314403966),
        T(-7801.945954775964),
        T(-26448.19586059192),
        T(-90799.48341621365),
        T(-315126.04064491636),
        T(-1.1040113443115912e6),
    ),
)
@inline _E9(x::T) where {T} = evalpoly(
    x,
    (
        T(1.1613071521962828),
        T(-0.7011002845552895),
        T(-0.5805514744654373),
        T(-1.2436930610777865),
        T(-3.679383613496635),
        T(-12.815909243378957),
        T(-49.25672530759985),
        T(-202.18187354340904),
        T(-869.8602699308701),
        T(-3877.0058473132895),
        T(-17761.7071017094),
        T(-83182.69029154233),
        T(-396650.4505013548),
        T(-1.9200334136826345e6),
    ),
)
@inline _E10(x::T) where {T} = evalpoly(
    x,
    (
        T(1.1246173251197522),
        T(-0.7708450563609095),
        T(-0.8447940536449113),
        T(-2.4900973094503946),
        T(-10.239717411543843),
        T(-49.7490054655148),
        T(-267.09866751957054),
        T(-1532.66588382523),
        T(-9222.313478526092),
        T(-57502.51612140314),
        T(-368596.11674161063),
        T(-2.4156110887010912e6),
        T(-1.6120097815816568e7),
        T(-1.0920993852030899e8),
        T(-7.493807581942496e8),
        T(-5.198725846725541e9),
        T(-3.64092568881214e10),
    ),
)

@inline function _Kseries(x::T) where {T}
    if x <= T(0.5)
        if x <= T(0.2)
            return x <= T(0.1) ? _K1(x - T(0.05)) : _K2(x - T(0.15))
        elseif x <= T(0.3)
            return _K3(x - T(0.25))
        elseif x <= T(0.4)
            return _K4(x - T(0.35))
        else
            return _K5(x - T(0.45))
        end
    elseif x <= T(0.8)
        if x <= T(0.6)
            return _K6(x - T(0.55))
        elseif x <= T(0.7)
            return _K7(x - T(0.65))
        else
            return _K8(x - T(0.75))
        end
    elseif x <= T(0.85)
        return _K9(x - T(0.825))
    elseif x <= T(0.9)
        return _K10(x - T(0.875))
    end
    # Fukushima (2009) Eq. 31: K(m) = -log(q(mc)) K(mc)/π, with K(mc) from the first fit.
    mc = one(T) - x
    return -log(_Kq(mc)) * (_K1(mc - T(0.05)) * T(ONE_DIV_PI))
end

@inline function _Eseries(x::T) where {T}
    if x <= T(0.5)
        if x <= T(0.2)
            return x <= T(0.1) ? _E1(x - T(0.05)) : _E2(x - T(0.15))
        elseif x <= T(0.3)
            return _E3(x - T(0.25))
        elseif x <= T(0.4)
            return _E4(x - T(0.35))
        else
            return _E5(x - T(0.45))
        end
    elseif x <= T(0.8)
        if x <= T(0.6)
            return _E6(x - T(0.55))
        elseif x <= T(0.7)
            return _E7(x - T(0.65))
        else
            return _E8(x - T(0.75))
        end
    elseif x <= T(0.85)
        return _E9(x - T(0.825))
    elseif x <= T(0.9)
        return _E10(x - T(0.875))
    end
    mc = one(T) - x
    kdm = _K1(mc - T(0.05))
    edm = _E1(mc - T(0.05))
    km = -log(_Kq(mc)) * (kdm * T(ONE_DIV_PI))
    return (T(π / 2) + km * (kdm - edm)) / kdm
end

"""
``K(m)=F(\\pi/2\\,|\\,m)``.

Returns the complete elliptic integral of the first kind.

# Arguments

- `m` : Elliptic parameter
"""
function K(m::T) where {T}
    isnan(m) && return T(NaN)
    m == one(T) && return T(Inf)

    # m < 0: K(m) = K(m/(m-1)) / √(1-m).
    negative = m < zero(T)
    x = negative ? m / (m - one(T)) : m

    x == zero(T) && return T(π / 2)
    x == one(T) && return T(Inf)
    x > one(T) && return T(NaN)

    t = _Kseries(x)
    return negative ? t / √(one(T) - m) : t
end

"""
``E(m)=E(\\pi/2\\,|\\,m)``.

Returns the complete elliptic integral of the second kind.

# Arguments

- `m` : Elliptic parameter
"""
function E(m::T) where {T}
    isnan(m) && return T(NaN)
    m == zero(T) && return T(π / 2)
    m == one(T) && return one(T)

    # A&S 17.4.16: E(m) = √(1-m) E(m/(m-1)) for m < 0.
    negative = m < zero(T)
    x = negative ? m / (m - one(T)) : m

    x == zero(T) && return T(π / 2)
    x == one(T) && return one(T)
    x > one(T) && return T(NaN)

    t = _Eseries(x)
    return negative ? t * √(one(T) - m) : t
end

############################################################################
# First Incomplete Elliptic integral and Inverse Jacobi Functions 
# (https://doi.org/T(10).1007/s00211-010-0321-8)
############################################################################
@inline function _serf(y::A, m::B) where {A,B}
    T = promote_type(A, B)
    y, m = T(y), T(m)
    return evalpoly(
        y,
        (
            T(1.0),
            evalpoly(m, (T(0.16666666666666666), T(0.16666666666666666))),
            evalpoly(m, (T(0.075), T(0.05), T(0.075))),
            evalpoly(
                m,
                (
                    T(0.044642857142857144),
                    T(0.026785714285714284),
                    T(0.026785714285714284),
                    T(0.044642857142857144),
                ),
            ),
            evalpoly(
                m,
                (
                    T(0.030381944444444444),
                    T(0.017361111111111112),
                    T(0.015625),
                    T(0.017361111111111112),
                    T(0.030381944444444444),
                ),
            ),
            evalpoly(
                m,
                (
                    T(0.022372159090909092),
                    T(0.012428977272727272),
                    T(0.01065340909090909),
                    T(0.01065340909090909),
                    T(0.012428977272727272),
                    T(0.022372159090909092),
                ),
            ),
            evalpoly(
                m,
                (
                    T(0.017352764423076924),
                    T(0.00946514423076923),
                    T(0.007887620192307692),
                    T(0.007512019230769231),
                    T(0.007887620192307692),
                    T(0.00946514423076923),
                    T(0.017352764423076924),
                ),
            ),
            evalpoly(
                m,
                (
                    T(0.01396484375),
                    T(0.00751953125),
                    T(0.00615234375),
                    T(0.005696614583333333),
                    T(0.005696614583333333),
                    T(0.00615234375),
                    T(0.00751953125),
                    T(0.01396484375),
                ),
            ),
            evalpoly(
                m,
                (
                    T(0.011551800896139705),
                    T(0.006160960477941176),
                    T(0.004976160386029412),
                    T(0.0045237821691176475),
                    T(0.004398121553308824),
                    T(0.0045237821691176475),
                    T(0.004976160386029412),
                    T(0.006160960477941176),
                    T(0.011551800896139705),
                ),
            ),
            evalpoly(
                m,
                (
                    T(0.009761609529194078),
                    T(0.005167910927220395),
                    T(0.004134328741776316),
                    T(0.0037102950246710527),
                    T(0.0035416452508223685),
                    T(0.0035416452508223685),
                    T(0.0037102950246710527),
                    T(0.004134328741776316),
                    T(0.005167910927220395),
                    T(0.009761609529194078),
                ),
            ),
        ),
    )
end

# The descending Landen ladder for asn, working in `y = s^2`. Returns |asn(√y|m)|; the
# caller restores the sign.
@inline function _asn_y(y::A, m::B) where {A,B}
    T = promote_type(A, B)
    y, m = T(y), T(m)
    yA = T(0.04094) - T(0.00652) * m
    p = one(T)
    for _ = 1:10
        y < yA && return p * √y * _serf(y, m)
        y = y * inv((one(T) + √(one(T) - y)) * (one(T) + √(one(T) - m * y)))
        p += p
    end
    return T(NaN)
end

"""
``\\text{asn}(u|\\,m)=\\text{sn}^{-1}(u,m)``.

Returns the inverse Jacobi Elliptic sn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic parameter
"""
function asn(s::A, m::B) where {A,B}
    T = promote_type(A, B)
    isnan(m) && return T(NaN)
    !(-1 ≤ m ≤ 1) && throw(DomainError("argument m not in [0,1]"))

    s, m = T(s), T(m)
    y = s * s
    # Only the direct series carries the sign of s; the descending ladder returns |asn|.
    y < T(0.04094) - T(0.00652) * m && return s * _serf(y, m)
    return _asn_y(y, m)
end

"""
``\\text{acn}(u|\\,m)=\\text{cn}^{-1}(u,m)``.

Returns the inverse Jacobi Elliptic cn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic parameter
"""
function acn(c::A, m::B) where {A,B}
    T = promote_type(A, B)
    isnan(m) && return T(NaN)
    !(-1 ≤ m ≤ 1) && throw(DomainError("argument m not in [0,1]"))

    c, m = T(c), T(m)
    mc = one(T) - m
    x = c * c
    p = one(T)
    for _ = 1:10
        2x > one(T) && return p * _asn_y(one(T) - x, m)

        d = √(mc + m * x)
        x = (√x + d) / (one(T) + d)
        p += p
    end
    return T(NaN)
end

@inline function _rawF(sinφ::A, m::B) where {A,B}
    T = promote_type(A, B)

    yS = T(0.9000308778823196)
    m == 0 && return asin(sinφ)
    m == 1 && return atanh(sinφ)

    sinφ2 = sinφ * sinφ
    sinφ2 ≤ yS && return asn(sinφ, m)

    oneT = one(T)
    mc = oneT - m

    c = √(oneT - sinφ2)
    x = c * c
    d2 = muladd(m, x, mc)
    # z = c/√(mc + m*c*c), so z^2 = x/d2.
    x < yS * d2 && return K(m) - _asn_y(x / d2, m)

    v = mc * (oneT - x)
    v < x * d2 && return acn(c, m)
    return K(m) - acn(√(v / d2), m)
end

#----------------------------------------------------------------------------------------
# Elliptic F
#----------------------------------------------------------------------------------------
function _F(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    abs(φ) ≤ T(π / 2) && return sign(φ) * _rawF(sin(φ), m)
    j = floor(φ / T(π))
    newφ = φ - j * T(π)
    signφ = sign(newφ)
    if abs(newφ) > T(π / 2)
        j += signφ * one(T)
        newφ = newφ - signφ * T(π)
    end
    signφ = sign(newφ)

    return 2j * K(m) + signφ * _rawF(sin(abs(newφ)), m)
end

function F(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)
    if m > 1
        ## Abramowitz & Stegum*(17.4.15)
        m12 = √(m)
        theta = asin(m12 * sin(φ))
        signθ = sign(theta)
        absθ = abs(theta)
        return signθ / m12 * _F(absθ, inv(m))
    elseif m < 0
        # Abramowitz & Stegum*(17.4.17)
        n = -m
        one_plus_n = oneT + n
        m12 = inv(√(one_plus_n))
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

#----------------------------------------------------------------------------------------
# Elliptic D
#----------------------------------------------------------------------------------------
function _D(m::T) where {T}
    return _cel2(√(one(T) - m), zero(T), one(T))
end

@inline function _Dsi(y::A, m::B) where {A,B}
    T = promote_type(A, B)
    y, m = T(y), T(m)
    return √(y) *
           y *
           evalpoly(
               y,
               (
                   T(0.3333333333333333),
                   evalpoly(m, (T(0.1), T(0.1))),
                   evalpoly(
                       m,
                       (
                           T(0.05357142857142857),
                           T(0.03571428571428571),
                           T(0.05357142857142857),
                       ),
                   ),
                   evalpoly(
                       m,
                       (
                           T(0.034722222222222224),
                           T(0.020833333333333332),
                           T(0.020833333333333332),
                           T(0.034722222222222224),
                       ),
                   ),
                   evalpoly(
                       m,
                       (
                           T(0.024857954545454544),
                           T(0.014204545454545454),
                           T(0.01278409090909091),
                           T(0.014204545454545454),
                           T(0.024857954545454544),
                       ),
                   ),
                   evalpoly(
                       m,
                       (
                           T(0.01893028846153846),
                           T(0.010516826923076924),
                           T(0.009014423076923076),
                           T(0.009014423076923076),
                           T(0.010516826923076924),
                           T(0.01893028846153846),
                       ),
                   ),
                   evalpoly(
                       m,
                       (
                           T(0.0150390625),
                           T(0.008203125),
                           T(0.0068359375),
                           T(0.006510416666666667),
                           T(0.0068359375),
                           T(0.008203125),
                           T(0.0150390625),
                       ),
                   ),
                   evalpoly(
                       m,
                       (
                           T(0.012321920955882353),
                           T(0.006634880514705882),
                           T(0.005428538602941176),
                           T(0.0050264246323529415),
                           T(0.0050264246323529415),
                           T(0.005428538602941176),
                           T(0.006634880514705882),
                           T(0.012321920955882353),
                       ),
                   ),
                   evalpoly(
                       m,
                       (
                           T(0.01033582185444079),
                           T(0.005512438322368421),
                           T(0.004452354029605263),
                           T(0.004047594572368421),
                           T(0.0039351613898026315),
                           T(0.004047594572368421),
                           T(0.004452354029605263),
                           T(0.005512438322368421),
                           T(0.01033582185444079),
                       ),
                   ),
               ),
           )
end

@inline function _Ds_cutoff(m::T) where {T}
    y0 = T == Float32 ? T(0.1888) : T(0.04094)
    m < zero(T) && return y0 / (one(T) - m)
    return y0 - (T == Float32 ? T(0.0378) : T(0.00652)) * m
end

function _Ds(s::A, m::B) where {A,B}
    T = promote_type(A, B)
    s, m = T(s), T(m)
    yA = _Ds_cutoff(m)

    y = s * s
    w = one(T)
    acc = zero(T)
    for _ = 1:10
        y < yA && break
        c = √(one(T) - y)
        d = √(one(T) - m * y)
        ynext = y / ((one(T) + c) * (one(T) + d))
        acc = muladd(w, √(y) * ynext, acc)
        w += w
        y = ynext
    end
    return muladd(w, _Dsi(y, m), acc)
end

function _rawD(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    ys = T == Float32 ? T(0.95) : T(0.9) # T(0.95) for single
    φs = T == Float32 ? T(1.345) : T(1.249) # T(1.345) for single
    s = sin(φ)

    if φ < φs
        return _Ds(s, m)
    end

    mc = one(T) - m
    c = cos(φ)
    z = c / √(mc + m * c^2)
    if z^2 < ys
        return _D(m) - _Ds(z, m) - s * z
    end

    w = √(1 - z^2)
    if c < w
        return _Ds(s, m)
    end

    return _D(m) - _Ds(√(one(T) - w^2), m) - √((one(T) - c^2) * (one(T) - w^2))
end

# https://www.sciencedirect.com/science/article/pii/S0377042711001270?ref=pdf_download&fr=RR-2&rr=80442ff79ecc4cef
function _D(φ::A, m::B) where {A,B}
    T = promote_type(A, B)

    # Reduction of Amplitude
    if abs(φ) > T(π / 2) && m < one(T)
        j = floor(φ / T(π))
        newφ = φ - j * T(π)
        signφ = sign(newφ)
        if abs(newφ) > T(π / 2)
            j += signφ * one(T)
            newφ = newφ - signφ * T(π)
        end
        signφ = sign(newφ)

        return 2 * j * _D(m) + signφ * _D(abs(newφ), m)
    end

    φ < zero(T) && return -_D(-φ, m)

    # Reduction of Parameter
    if m > one(T) && m * sin(φ)^2 ≤ one(T)
        mr = inv(m)
        φr = asin(√m * sin(φ))
        return mr * √mr * _D(φr, mr)
    end

    # if m < zero(T)
    #	mN = -m / (1 - m)
    #	φN = asin(√(mc / (1 - m * sin(φ)^2)) * sin(φ))
    #	return √(1 - mN) * (_D(φN, mN) - sin(φN) * cos(φN) / √(1 - mN * sin(φN)^2))
    # end


    return _rawD(φ, m)
end

#----------------------------------------------------------------------------------------
# Elliptic E
#----------------------------------------------------------------------------------------

# The three amplitudes at which E(φ|m) is in closed-form.
@inline function _E_closed_form(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    m == one(T) && abs(φ) <= T(π / 2) && return (true, sin(T(φ)))
    φ == T(π / 2) && return (true, E(T(m)))
    φ == -T(π / 2) && return (true, -E(T(m)))
    return (false, zero(T))
end

function E(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    (isnan(φ) || isnan(m)) && return T(NaN)
    closed, e = _E_closed_form(φ, m)
    closed && return e
    return F(φ, m) - m * _D(φ, m)
end

"""
``(F(\\varphi|\\,m),\\ E(\\varphi|\\,m))`` from a single evaluation of ``F``.

# Arguments

- `φ` : Amplitude
- `m` : Elliptic parameter
"""
@inline function _F_E(φ::A, m::B) where {A,B}
    T = promote_type(A, B)
    (isnan(φ) || isnan(m)) && return (T(NaN), T(NaN))
    closed, e = _E_closed_form(φ, m)
    closed && return (F(φ, m), e)
    f = F(φ, m)
    return (f, f - m * _D(φ, m))
end

#----------------------------------------------------------------------------------------
# Elliptic Π
#----------------------------------------------------------------------------------------
Π(n, m) = Pi(n, m)
Π(n, φ, m) = Pi(n, φ, m)
# https://doi.org/T(10).1016/j.cam.2011.1107

function Pi(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    oneT = one(T)
    zeroT = zero(T)
    m > oneT && return T(NaN) # Complex Branch
    n > oneT && return K(m) - Pi(m / n, m) # Cauchy principal value
    n == oneT && return T(Inf) # 1/(1 - n*sin(θ)^2) pole at θ = π/2
    n == zeroT && return K(m)
    m == oneT && return T(Inf)
    kc = √(oneT - m)
    nc = oneT - n
    return _cel(kc, nc, oneT, oneT)
end

function Pi(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    oneT = one(T)

    if m < 0 # Imaginary modulus transformation https://dlmf.nist.gov/19.7# iii
        if n == m # the transformation below is 0/0 at n == m, so use the closed form
            sinφ, cosφ = sincos(φ)
            return (E(φ, m) - m * sinφ * cosφ / √(muladd(-m, sinφ * sinφ, oneT))) /
                   (oneT - m)
        end
        # the transformation goes through asin, which folds |φ| > π/2 onto the principal
        # branch, so reduce first with Π(n; φ + kπ | m) = Π(n; φ | m) + 2k Π(n | m)
        if 2abs(φ) > T(π)
            φ2 = φ + T(π / 2)
            return 2fld(φ2, T(π)) * Pi(n, m) + Pi(n, mod(φ2, T(π)) - T(π / 2), m)
        end
        mc = oneT - m
        imc = inv(mc)
        mN = -m * imc
        sinφ = sin(φ)
        sinφ2 = sinφ * sinφ
        φN = asin(√(mc / muladd(-m, sinφ2, oneT)) * sinφ)

        nN = (n - m) * imc

        return √(imc) / nN * (mN * F(φN, mN) + imc * n * Pi(nN, φN, mN))
    end # https://link.springer.com/book/10.1007/978-3-642-65138-0 117.01
    if n > 1
        nc = oneT - n
        sinφ, cosφ = sincos(φ)
        sinφ2 = sinφ * sinφ
        t1 = sinφ / (cosφ * √(muladd(-m, sinφ2, oneT)))
        h1 = nc * (n−m) / n
        n1 = m / n
        return (_FukushimaT(t1, h1) - n1 * J(n1, φ, m))
    elseif n == 1
        φ == T(π/2) && return T(Inf)
        if m == 1
            sinφ, cosφ = sincos(φ)
            secφ = inv(cosφ)
            return 1 / 2 * (log(secφ * (1 + sinφ)) + secφ * sinφ / cosφ)
        else
            sinφ, cosφ = sincos(φ)
            sinφ2 = sinφ * sinφ
            sqrt_term = √(muladd(-m, sinφ2, oneT))
            return F(φ, m) + (-E(φ, m) + sqrt_term * sinφ / cosφ) / (oneT - m)
        end
    elseif n == 0
        return F(φ, m)
    end


    return n * J(n, φ, m) + F(φ, m)
end

"""
``J (n;\\varphi \\,|\\,m)=\\frac{\\Pi(n;\\varphi|\\, m) - F(\\varphi|\\,m)}{n}.``

Returns the associate incomplete elliptic integral of the third kind.

# Arguments

- `n` : Characteristic
- `φ` : Amplitude
- `m` : Elliptic parameter
"""
function J(n::A, φ::B, m::C) where {A,B,C} # Appendix A
    T = promote_type(A, B, C)
    oneT = one(T)
    # Reduction of Amplitude
    φ == zero(T) && return zero(T)
    φ == T(π / 2) && return J(n, m)

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
    φ < zero(T) && return -J(n, -φ, m)

    # Reduction of parameter
    if zero(T) < φ < T(π / 2)
        sinφ, cosφ = sincos(φ)
        sinφ2 = sinφ * sinφ
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
        return _rawJ(n, φ, m)
    return T(NaN)
end

function _rawJ(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    φs = T == Float32 ? T(1.345) : T(1.249) # T(1.345) for single
    s = sin(φ)
    φ < φs && return _Js(n, s, m)

    ys = T == Float32 ? T(0.95) : T(0.9) # T(0.95) for single
    mc = one(T) - m
    nc = one(T) - n
    h = n * nc * (n - m)
    c = cos(φ)
    z = c / √(mc + m * c * c)
    z * z < ys && return J(n, m) - _Js(n, z, m) - _FukushimaT(s * z / nc, h)

    w = √(mc / (one(T) - m * s * s)) * s
    c < w && return _Js(n, s, m)

    wc2 = one(T) - w * w
    return J(n, m) - _Js(n, √(wc2), m) - _FukushimaT(√((one(T) - c * c) * wc2) / nc, h)
end

"""
``J(n;\\varphi \\,|\\,m)=\\frac{\\Pi(n;\\pi/2|\\, m) - K(m)}{n}.``

Returns the associate complete elliptic integral of the third kind.

# Arguments

- `n` : Characteristic
- `m` : Elliptic parameter
"""
function J(n::A, m::B) where {A,B}
    T = promote_type(A, B)
    n > one(T) && return m / n * J(m / n, m)
    kc = √(1 - m)
    nc = 1 - n
    return _cel(kc, nc, zero(T), one(T))
end

# function Jc(n::T, c::T, m::T) where{T}
#    # _ybuf = @SVector [zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T)]#
#    _xbuf = @SVector [zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T)]#
#
#    nc = one(T) - n
#    mc = one(T) - m
#    h = n*nc*(n-m)
#
#    x0 = c^2
#   
#    xi = x0 
#
#    xB = T == Float32 ? T(0.94095) : T(0.98378)
#    I = 0
#    for _ in 1:10
#        if xi < xB
#            break
#        end
#        I += 1
#        # @set! _ybuf[I] = one(T) - xi
#        @set! _xbuf[I] = xi
#        ci = √xi
#        di = √(mc+m*xi)
#        xi = (ci +di)/(one(T)+di)
#    end
#    J0 = _JsI(n, one(T)-xi, m)
#    I == 0 && return J0
#    Ji = J0
#    yi = one(T) - xi 
#    ti = 0
#    for i in I:-1:1
#        # yim1 = _ybuf[i]
#        # sim1 = √yim1
#        # cim1 = √(one(T)-yim1)
#        # dim1 = √(one(T)-m*yim1)
#        
#        # ti = sim1*yi / (one(T)-n*(yim1 - cim1*dim1*yi))
#        # Ji = 2Ji + _FukushimaT(ti, h)
#
#        # yi = yim1
#
#        xim1 = _xbuf[i]
#        yim1 = one(T) - xim1
#        sim1 = √yim1
#        cim12 = (xim1)
#        dim12 = (one(T)-m*yim1)
#        
#        ti = sim1*yi / (one(T)-n*(yim1 - √(cim1*dim1)*yi))
#        Ji = 2Ji + _FukushimaT(ti, h)
#
#        xi = xim1
#
#    end
#    return Ji
#
# end

function _Js(n::A, s::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n, s, m = T(n), T(s), T(m)
    nc = one(T) - n
    h = n * nc * (n - m)
    yB = T == Float32 ? T(0.05905) : T(0.01622)

    y = s * s
    w = one(T)
    acc = zero(T)
    for _ = 1:10
        y < yB && break
        c = √(one(T) - y)
        d = √(one(T) - m * y)
        ynext = y / ((one(T) + c) * (one(T) + d))
        t = √(y) * ynext / (one(T) - n * (y - c * d * ynext))
        acc = muladd(w, _FukushimaT(t, h), acc)
        w += w
        y = ynext
    end
    return muladd(w, _JsI(n, y, m), acc)
end

@inline function _JsI(n::A, y::B, m::C) where {A,B,C}
    T = promote_type(A, B, C)
    n, y, m = T(n), T(y), T(m)
    return √(y) *
           y *
           evalpoly(
               y,
               (
                   T(0.3333333333333333),
                   evalpoly(n, (evalpoly(m, (T(0.1), T(0.1))), T(0.2))),
                   evalpoly(
                       n,
                       (
                           evalpoly(
                               m,
                               (
                                   T(0.05357142857142857),
                                   T(0.03571428571428571),
                                   T(0.05357142857142857),
                               ),
                           ),
                           evalpoly(m, (T(0.07142857142857142), T(0.07142857142857142))),
                           T(0.14285714285714285),
                       ),
                   ),
                   evalpoly(
                       n,
                       (
                           evalpoly(
                               m,
                               (
                                   T(0.034722222222222224),
                                   T(0.020833333333333332),
                                   T(0.020833333333333332),
                                   T(0.034722222222222224),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.041666666666666664),
                                   T(0.027777777777777776),
                                   T(0.041666666666666664),
                               ),
                           ),
                           evalpoly(m, (T(0.05555555555555555), T(0.05555555555555555))),
                           T(0.1111111111111111),
                       ),
                   ),
                   evalpoly(
                       n,
                       (
                           evalpoly(
                               m,
                               (
                                   T(0.024857954545454544),
                                   T(0.014204545454545454),
                                   T(0.01278409090909091),
                                   T(0.014204545454545454),
                                   T(0.024857954545454544),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.028409090909090908),
                                   T(0.017045454545454544),
                                   T(0.017045454545454544),
                                   T(0.028409090909090908),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.03409090909090909),
                                   T(0.022727272727272728),
                                   T(0.03409090909090909),
                               ),
                           ),
                           evalpoly(m, (T(0.045454545454545456), T(0.045454545454545456))),
                           T(0.09090909090909091),
                       ),
                   ),
                   evalpoly(
                       n,
                       (
                           evalpoly(
                               m,
                               (
                                   T(0.01893028846153846),
                                   T(0.010516826923076924),
                                   T(0.009014423076923076),
                                   T(0.009014423076923076),
                                   T(0.010516826923076924),
                                   T(0.01893028846153846),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.021033653846153848),
                                   T(0.01201923076923077),
                                   T(0.010817307692307692),
                                   T(0.01201923076923077),
                                   T(0.021033653846153848),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.02403846153846154),
                                   T(0.014423076923076924),
                                   T(0.014423076923076924),
                                   T(0.02403846153846154),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.028846153846153848),
                                   T(0.019230769230769232),
                                   T(0.028846153846153848),
                               ),
                           ),
                           evalpoly(m, (T(0.038461538461538464), T(0.038461538461538464))),
                           T(0.07692307692307693),
                       ),
                   ),
                   evalpoly(
                       n,
                       (
                           evalpoly(
                               m,
                               (
                                   T(0.0150390625),
                                   T(0.008203125),
                                   T(0.0068359375),
                                   T(0.006510416666666667),
                                   T(0.0068359375),
                                   T(0.008203125),
                                   T(0.0150390625),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.01640625),
                                   T(0.009114583333333334),
                                   T(0.0078125),
                                   T(0.0078125),
                                   T(0.009114583333333334),
                                   T(0.01640625),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.018229166666666668),
                                   T(0.010416666666666666),
                                   T(0.009375),
                                   T(0.010416666666666666),
                                   T(0.01818181818181818),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.020833333333333332),
                                   T(0.0125),
                                   T(0.0125),
                                   T(0.020833333333333332),
                               ),
                           ),
                           evalpoly(m, (T(0.025), T(0.016666666666666666), T(0.025))),
                           evalpoly(m, (T(0.03333333333333333), T(0.03333333333333333))),
                           T(0.06666666666666667),
                       ),
                   ),
                   evalpoly(
                       n,
                       (
                           evalpoly(
                               m,
                               (
                                   T(0.012321920955882353),
                                   T(0.006634880514705882),
                                   T(0.005428538602941176),
                                   T(0.0050264246323529415),
                                   T(0.0050264246323529415),
                                   T(0.005428538602941176),
                                   T(0.006634880514705882),
                                   T(0.012321920955882353),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.013269761029411764),
                                   T(0.007238051470588236),
                                   T(0.00603170955882353),
                                   T(0.0057444852941176475),
                                   T(0.00603170955882353),
                                   T(0.007238051470588236),
                                   T(0.013269761029411764),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.014476102941176471),
                                   T(0.008042279411764705),
                                   T(0.006893382352941176),
                                   T(0.006893382352941176),
                                   T(0.008042279411764705),
                                   T(0.014476102941176471),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.01608455882352941),
                                   T(0.009191176470588236),
                                   T(0.008272058823529412),
                                   T(0.009191176470588236),
                                   T(0.01608455882352941),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.01838235294117647),
                                   T(0.011029411764705883),
                                   T(0.011029411764705883),
                                   T(0.01838235294117647),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.022058823529411766),
                                   T(0.014705882352941176),
                                   T(0.022058823529411766),
                               ),
                           ),
                           evalpoly(m, (T(0.029411764705882353), T(0.029411764705882353))),
                           T(0.058823529411764705),
                       ),
                   ),
                   evalpoly(
                       n,
                       (
                           evalpoly(
                               m,
                               (
                                   T(0.01033582185444079),
                                   T(0.005512438322368421),
                                   T(0.004452354029605263),
                                   T(0.004047594572368421),
                                   T(0.0039351613898026315),
                                   T(0.004047594572368421),
                                   T(0.004452354029605263),
                                   T(0.005512438322368421),
                                   T(0.01033582185444079),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.011024876644736841),
                                   T(0.005936472039473684),
                                   T(0.004857113486842105),
                                   T(0.004497327302631579),
                                   T(0.004497327302631579),
                                   T(0.004857113486842105),
                                   T(0.005936472039473684),
                                   T(0.011024876644736841),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.011872944078947368),
                                   T(0.006476151315789474),
                                   T(0.005396792763157895),
                                   T(0.005139802631578948),
                                   T(0.005396792763157895),
                                   T(0.006476151315789474),
                                   T(0.011872944078947368),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.012952302631578948),
                                   T(0.007195723684210526),
                                   T(0.006167763157894737),
                                   T(0.006167763157894737),
                                   T(0.007195723684210526),
                                   T(0.012952302631578948),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.014391447368421052),
                                   T(0.008223684210526315),
                                   T(0.007401315789473684),
                                   T(0.008223684210526315),
                                   T(0.014391447368421052),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.01644736842105263),
                                   T(0.009868421052631578),
                                   T(0.009868421052631578),
                                   T(0.01644736842105263),
                               ),
                           ),
                           evalpoly(
                               m,
                               (
                                   T(0.019736842105263157),
                                   T(0.013157894736842105),
                                   T(0.019736842105263157),
                               ),
                           ),
                           evalpoly(m, (T(0.02631578947368421), T(0.02631578947368421))),
                           T(0.05263157894736842),
                       ),
                   ),
               ),
           )
end
# https://doi-org.ezp-prod1.hul.harvard.edu/T(10).031007/s10569-008-9177-y
# https://link.springer.com/article/T(10).1007/s10569-008-9177-y

function _cel2(kc::A, a::B, b::C) where {A,B,C}
    T = promote_type(A, B, C)
    return _cel(kc, one(T), a, b)
end

function _Kscreen(m::T) where {T}
    return T(π / 2) * (
        one(T) +
        m * (T(0.25) + m * (T(0.36) + m * (T(0.09765625) + m * T(0.07476806640625))))
    )
end

# Number of halvings that brings u into the radius where the four-term series in
# `_ΔXNloop` are accurate, i.e. max(6 + floor(log2(u)), 1). `frexp` splits a finite u > 0
@inline function _exp2_neg(::Type{T}, n::Integer) where {T}
    r = one(T)
    h = T(0.5)
    while n > 0
        isodd(n) && (r *= h)
        h *= h
        n >>= 1
    end
    return r
end

@inline function _dup_steps(u::T) where {T}
    u > zero(T) || return 0
    _, e = frexp(u)
    return max(5 + e, 1)
end

function _ΔXNloop(u::A, m::B, n::Integer) where {A,B}
    T = promote_type(A, B)
    m = T(m)
    oneT = one(T)
    twoT = T(2)

    up = T(u) * _exp2_neg(T, n)
    up2 = up * up

    Δsn =
        -up *
        up2 *
        evalpoly(
            up2,
            (
                evalpoly(m, (T(-1 / 6), T(-1 / 6))),
                evalpoly(m, (T(1 / 120), T(7 / 60), T(1 / 120))),
                evalpoly(m, (T(-1 / 5040), T(-3 / 112), T(-3 / 112), T(-1 / 5040))),
            ),
        )
    Δcn =
        -up2 * evalpoly(
            up2,
            (
                T(-1 / 2),
                evalpoly(m, (T(1 / 24), T(1 / 6))),
                evalpoly(m, (T(-1 / 720), T(-11 / 180), T(-1 / 45))),
                evalpoly(m, (T(-1 / 40320), T(-17 / 1680), T(-19 / 840), T(-1 / 630))),
            ),
        )
    Δdn =
        -up2 *
        m *
        evalpoly(
            up2,
            (
                T(-1 / 2),
                evalpoly(m, (T(1 / 6), T(1 / 24))),
                evalpoly(m, (T(-1 / 45), T(-11 / 180), T(-1 / 720))),
            ),
        )

    sn = up - Δsn
    cn = oneT - Δcn
    dn = oneT - Δdn

    for _ = 1:n
        sn2 = sn * sn
        sn4 = sn2 * sn2

        den = inv(oneT - m * sn4)
        Δsn = twoT * (Δsn * cn * dn + up * (Δcn * (oneT - Δdn) + Δdn - m * sn4)) * den
        Δcn = ((oneT + cn) * Δcn + sn2 - twoT * m * sn4) * den
        Δdn = ((oneT + dn) * Δdn + m * (sn2 - twoT * sn4)) * den

        up *= twoT
        sn = up - Δsn
        cn = oneT - Δcn
        dn = oneT - Δdn
    end
    return sn, cn, dn
end

function _fold_0_25(u1::A, m::B, kp::C) where {A,B,C}
    T = promote_type(A, B, C)
    oneT = one(T)
    u1 == 0 && return zero(T), oneT, oneT

    sn, cn, dn = _ΔXNloop(u1, m, _dup_steps(u1))
    one_plus_kp = oneT + kp
    den = inv(one_plus_kp - m * sn * sn)
    sqrt_one_plus_kp = √(one_plus_kp)
    sqrt_kp = √kp
    return den * sqrt_one_plus_kp * (cn * dn - kp * sn),
    den * sqrt_kp * sqrt_one_plus_kp * (cn + sn * dn),
    den * sqrt_kp * (one_plus_kp * dn + m * sn * cn)
end

function _fold_0_50(u1::A, m::B, Kscreen::C, Kactual::D, kp::E) where {A,B,C,D,E}
    T = promote_type(A, B, C, D, E)
    oneT = one(T)
    u1 == 0 && return zero(T), oneT, oneT

    if u1 > T(0.25) * Kscreen
        sn, cn, dn = _fold_0_25(Kactual / 2 - u1, m, kp)
    else
        sn, cn, dn = _ΔXNloop(u1, m, _dup_steps(u1))
    end
    invdn = inv(dn)
    return cn * invdn, kp * sn * invdn, kp * invdn
end

function _fold_1_00(u1::A, m::B, Kscreen::C, Kactual::D, kp::E) where {A,B,C,D,E}
    T = promote_type(A, B, C, D, E)
    oneT = one(T)
    u1 == 0 && return zero(T), oneT, oneT

    if u1 > Kscreen / 2
        sn, cn, dn = _fold_0_50(Kactual - u1, m, Kscreen, Kactual, kp)
    elseif u1 > Kscreen / 4
        sn, cn, dn = _fold_0_25(Kactual / 2 - u1, m, kp)
    else
        sn, cn, dn = _ΔXNloop(u1, m, _dup_steps(u1))
    end
    invdn = inv(dn)
    return cn * invdn, -kp * sn * invdn, kp * invdn
end

function _SN(u::T, m) where {T}
    tempK = K(m)
    kp = √(one(T) - m)
    fourK = 4 * tempK
    if u > fourK
        return _rawSN(u % fourK, m, tempK, tempK, kp)
    end
    return _rawSN(u, m, tempK, tempK, kp)
end

function _rawSN(u::A, m::B, Kscreen::C, Kactual::D, kp::E) where {A,B,C,D,E}
    T = promote_type(A, B, C, D, E)
    oneT = one(T)
    halfK = Kscreen / 2
    quarterK = Kscreen / 4
    check = u ≥ 2 * Kscreen
    sign = check ? -1 : 1
    u = check ? u - 2 * Kactual : u
    u > Kscreen && return sign * _fold_1_00(u - Kactual, m, Kscreen, Kactual, kp)[1]
    u == Kscreen && return sign * oneT
    u > halfK && return sign * _fold_0_50(Kactual - u, m, Kscreen, Kactual, kp)[1]
    u == halfK && return sign * inv(√(oneT + kp))
    u ≥ quarterK && return sign * _fold_0_25(Kactual / 2 - u, m, kp)[1]
    return sign * _ΔXNloop(u, m, _dup_steps(u))[1]
end

function _CN(u, m)
    tempK = K(m)
    kp = √(1 - m)
    fourK = 4 * tempK
    if u > fourK
        return _rawCN(u % fourK, m, tempK, tempK, kp)
    end
    return _rawCN(u, m, tempK, tempK, kp)
end
function _rawCN(u::A, m::B, Kscreen::C, Kactual::D, kp::E) where {A,B,C,D,E}
    T = promote_type(A, B, C, D, E)
    halfK = Kscreen / 2
    quarterK = Kscreen / 4
    check = u ≥ 2 * Kscreen
    sign = check ? -1 : 1
    u = check ? u - 2 * Kactual : u
    u > Kscreen && return sign * _fold_1_00(u - Kactual, m, Kscreen, Kactual, kp)[2]
    u == Kscreen && return sign * zero(T)
    u > halfK && return sign * _fold_0_50(Kactual - u, m, Kscreen, Kactual, kp)[2]
    u == halfK && return sign * √(kp / (one(T) + kp))
    u ≥ quarterK && return sign * _fold_0_25(Kactual / 2 - u, m, kp)[2]
    return sign * _ΔXNloop(u, m, _dup_steps(u))[2]
end

function _DN(u, m)
    tempK = K(m)
    kp = √(1 - m)
    fourK = 4 * tempK
    if u > fourK
        return _rawDN(u % fourK, m, tempK, tempK, kp)
    end
    return _rawDN(u, m, tempK, tempK, kp)
end
function _rawDN(u::A, m::B, Kscreen::C, Kactual::D, kp::E) where {A,B,C,D,E}
    halfK = Kscreen / 2
    quarterK = Kscreen / 4
    check = u ≥ 2 * Kscreen
    u = check ? u - 2 * Kactual : u
    u > Kscreen && return _fold_1_00(u - Kactual, m, Kscreen, Kactual, kp)[3]
    u == Kscreen && return kp
    u > halfK && return _fold_0_50(Kactual - u, m, Kscreen, Kactual, kp)[3]
    u == halfK && return √kp
    u ≥ quarterK && return _fold_0_25(Kactual / 2 - u, m, kp)[3]
    return _ΔXNloop(u, m, _dup_steps(u))[3]
end

# `_SN`, `_CN` and `_DN` each evaluate `K(m)` and then the full `(sn, cn, dn)` triple
# only to index one element out of it. Callers that want more than one of them pay for
# that three times over, so this hands back the whole triple from one pass.
function _SCD(u::T, m) where {T}
    tempK = K(m)
    kp = √(one(T) - m)
    fourK = 4 * tempK
    if u > fourK
        return _rawSCD(u % fourK, m, tempK, tempK, kp)
    end
    return _rawSCD(u, m, tempK, tempK, kp)
end

function _rawSCD(u::A, m::B, Kscreen::C, Kactual::D, kp::E) where {A,B,C,D,E}
    T = promote_type(A, B, C, D, E)
    oneT = one(T)
    halfK = Kscreen / 2
    quarterK = Kscreen / 4
    check = u ≥ 2 * Kscreen
    # `dn` is even about `2K` where `sn` and `cn` flip, so the fold sign reaches only
    # the first two components.
    sign = check ? -oneT : oneT
    u = check ? u - 2 * Kactual : u
    if u > Kscreen
        s, c, d = _fold_1_00(u - Kactual, m, Kscreen, Kactual, kp)
        return (sign * s, sign * c, d)
    end
    u == Kscreen && return (sign * oneT, sign * zero(T), T(kp))
    if u > halfK
        s, c, d = _fold_0_50(Kactual - u, m, Kscreen, Kactual, kp)
        return (sign * s, sign * c, d)
    end
    u == halfK && return (sign * inv(√(oneT + kp)), sign * √(kp / (oneT + kp)), √kp)
    if u ≥ quarterK
        s, c, d = _fold_0_25(Kactual / 2 - u, m, kp)
        return (sign * s, sign * c, d)
    end
    s, c, d = _ΔXNloop(u, m, _dup_steps(u))
    return (sign * s, sign * c, d)
end

_sc_helper(jacobituple) = jacobituple[1] / jacobituple[2]
function _rawSC(u::A, m::B, Kscreen::C, Kactual::D, kp::E) where {A,B,C,D,E}
    T = promote_type(A, B, C, D, E)
    fourK = 4 * Kscreen
    u = u ≥ fourK ? u % (4 * Kactual) : u
    u = u ≥ 2 * Kscreen ? u - 2 * Kactual : u
    u > Kscreen && return _sc_helper(_fold_1_00(u - Kactual, m, Kscreen, Kactual, kp))
    u == Kscreen && return T(Inf)
    u > Kscreen / 2 && return _sc_helper(_fold_0_50(Kactual - u, m, Kscreen, Kactual, kp))
    u == Kscreen / 2 && return inv(√(kp))
    u ≥ Kscreen / 4 && return _sc_helper(_fold_0_25(Kactual / 2 - u, m, kp))
    return _sc_helper(_ΔXNloop(u, m, _dup_steps(u)))
end

function _SC(u, m)
    tempK = K(m)
    return _rawSC(u, m, tempK, tempK, √(1 - m))
end

_sd_helper(jacobituple) = jacobituple[1] / jacobituple[3]
function _rawSD(u::A, m::B, Kscreen::C, Kactual::D, kp::E) where {A,B,C,D,E}
    T = promote_type(A, B, C, D, E)
    fourK = 4 * Kscreen
    u = u ≥ fourK ? u % (4 * Kactual) : u
    # `sn` flips sign across 2K and `dn` does not, so unlike `_rawSC` -- where the two
    # flips cancel in the ratio -- `sd` carries the fold sign through.
    check = u ≥ 2 * Kscreen
    sign = check ? -one(T) : one(T)
    u = check ? u - 2 * Kactual : u
    u > Kscreen &&
        return sign * _sd_helper(_fold_1_00(u - Kactual, m, Kscreen, Kactual, kp))
    u == Kscreen && return sign * inv(kp)
    u > Kscreen / 2 &&
        return sign * _sd_helper(_fold_0_50(Kactual - u, m, Kscreen, Kactual, kp))
    u == Kscreen / 2 && return sign * inv(√((1 + kp) * kp))
    u ≥ Kscreen / 4 && return sign * _sd_helper(_fold_0_25(Kactual / 2 - u, m, kp))
    return sign * _sd_helper(_ΔXNloop(u, m, _dup_steps(u)))
end

function _SD(u, m)
    tempK = K(m)
    return _rawSD(u, m, tempK, tempK, √(1 - m))
end

"""
``\\text{sn}(u|\\,m)=\\sin\\,\\text{am}(u,m)``, where ``\\text{am}(u|\\,m)=F^{-1}(u|\\,m)`` 
is the Jacobi amplitude.

Returns the Jacobi Elliptic sn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic parameter
"""
function sn(u, m)
    signu = sign(u)
    u = abs(u)
    m < 1 && return signu * _SN(u, m)
    sqrtm = √m
    return signu * _SN(u * sqrtm, inv(m)) / sqrtm
end

"""
``\\text{cn}(u|\\,m)=\\cos\\,\\text{am}(u,m)``, where ``\\text{am}(u|\\,m)=F^{-1}(u|\\,m)`` 
is the Jacobi amplitude.

Returns the Jacobi Elliptic cn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic parameter
"""
function cn(u, m)
    u = abs(u)
    m < 1 && return _CN(u, m)
    sqrtm = √m
    return _DN(u * sqrtm, inv(m))
end

"""
``\\text{dn}(u|\\,m)=\\frac{d}{d u}\\,\\text{am}(u,m)``, where ``\\text{am}(u|\\,m)=F^{-1}(u|\\,m)`` 
is the Jacobi amplitude.

Returns the Jacobi Elliptic dn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic parameter
"""
function dn(u, m)
    u = abs(u)
    m < 1 && return _DN(u, m)
    sqrtm = √m
    return _CN(u * sqrtm, inv(m))
end

# `sn`, `cn`, `dn` and `am` at one argument, sharing the single `_SCD` pass that the
# public entry points would otherwise repeat three times. The reciprocal-modulus branch
# maps all three onto the same reduced argument, so it fuses too. `am` comes from the
# shared AGM kernel and is independent of the duplication loop.
@inline function _scd_am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    signu = sign(u)
    absu = abs(u)
    if m < one(T)
        s, c, d = _SCD(absu, m)
        return (signu * s, c, d, am(u, m))
    end
    sqrtm = √m
    s, c, d = _SCD(absu * sqrtm, inv(m))
    return (signu * s / sqrtm, d, c, am(u, m))
end

function sc(u, m)
    signu = sign(u)
    u = abs(u)

    m < 1 && return signu * _SC(u, m)
    sqrtm = √m
    return signu * _SD(u * sqrtm, inv(m)) / sqrtm
end

function sd(u, m)
    signu = sign(u)
    u = abs(u)

    m < 1 && return signu * _SD(u, m)
    sqrtm = √m
    return signu * _SC(u * sqrtm, inv(m)) * inv(sqrtm)
end

function cd(u, m)
    return cn(u, m) / dn(u, m)
end

end
