const HALF_PI32 = 1.5707964f0
const PI32 = 3.1415927f0
const ONE_DIV_PI32 = 0.3183098861837907f0

############################################################################
# Complete Elliptic Integrals (https://doi.org/10.1007/s10569f0-009f0-9228f0-z)
############################################################################
function K( m::Float32 )::Float32
    if m < 1.0f0
        poly1(x::Float32)  = evalpoly(x, (1.59100345f0, 0.41600074f0, 0.24579151f0, 0.17948148f0, 0.14455606f0))
        poly2(x::Float32)  = evalpoly(x, (1.63525673f0, 0.47119063f0, 0.30972841f0, 0.25220831f0, 0.22672562f0))
        poly3(x::Float32)  = evalpoly(x, (1.68575035f0, 0.54173185f0, 0.40152444f0, 0.36964247f0, 0.37606072f0))
        poly4(x::Float32)  = evalpoly(x, (1.74435060f0, 0.63486428f0, 0.53984256f0, 0.57189271f0, 0.67029514f0, 0.83258659f0))
        poly5(x::Float32)  = evalpoly(x, (1.81388394f0, 0.76316325f0, 0.76192861f0, 0.95107465f0, 1.31518068f0, 1.92856069f0))
        poly6(x::Float32)  = evalpoly(x, (1.89892491f0, 0.95052179f0, 1.15107759f0, 1.75023911f0, 2.95267681f0, 5.28580040f0))
        poly7(x::Float32)  = evalpoly(x, (2.00759840f0, 1.24845723f0, 1.92623466f0, 3.75128964f0, 8.11994455f0, 18.6657213f0))
        poly8(x::Float32)  = evalpoly(x, (2.15651565f0, 1.79180564f0, 3.82675129f0, 10.3867247f0, 31.4033141f0, 100.923704f0, 337.326828f0, 1158.70793f0))
        poly9(x::Float32)  = evalpoly(x, (2.31812262f0, 2.61692015f0, 7.89793508f0, 30.5023972f0, 131.486937f0, 602.984764f0, 2877.02462f0))
        poly10(x::Float32) = evalpoly(x, (2.47359617f0, 3.72762424f0, 15.6073930f0, 84.1285084f0, 506.981820f0, 3252.27706f0, 21713.2424f0, 149037.045f0))
        poly11(x::Float32) = evalpoly(x, (0.0f0, 0.0625f0, 0.03125f0, 0.020507812f0, 0.015136719f0, 0.01193428f0, 0.00981617f0))
        poly12(x::Float32) = evalpoly(x, (1.59100345f0, 0.41600074f0, 0.24579151f0, 0.17948148f0, 0.14455606f0))


        flag::Bool = false
        kdm = 0.0f0
        td = 0.0f0
        qd = 0.0f0
        t = 0.0f0

        x = m
        if m < 0.0f0
            x = m / ( m - 1.0f0 );
            flag = true;
        end
        x == 0.0f0 && return HALF_PI32
        x == 1.0f0 && return Inf32
        x > 1.0f0 && return NaN32
        
        if x < 0.1f0 
            t = poly1( x - 0.05f0 );
        elseif ( x < 0.2f0 ) 
            t = poly2( x - 0.15f0 );
        elseif ( x < 0.3f0 ) 
            t = poly3( x - 0.25f0 );
        elseif ( x < 0.4f0 ) 
            t = poly4( x - 0.35f0 );
        elseif ( x < 0.5f0 ) 
            t = poly5( x - 0.45f0 );
        elseif ( x < 0.6f0 ) 
            t = poly6( x - 0.55f0 );
        elseif ( x < 0.7f0 ) 
            t = poly7( x - 0.65f0 );
        elseif ( x < 0.8f0 ) 
            t = poly8( x - 0.75f0 );
        elseif ( x < 0.85f0 ) 
            t = poly9( x - 0.825f0 );
        elseif ( x < 0.9f0 ) 
            t = poly10( x - 0.875f0 );
        else 
            td = 1.0f0 - x;
            qd = poly11( td );
            kdm = poly12( td - 0.05f0 );
            t = -log( qd ) * ( kdm * ONE_DIV_PI32 );
        end
        # Complete the transformation mentioned above for m < 0:
        flag && return t / sqrt( 1.0f0 - m );
        
        return t
    end
    m == 1.0f0 && return Inf32
    return NaN32
end

function E(m::Float32)::Float32
    m == 0.0f0 && return HALF_PI32
    if m < 1.0f0
        poly1(x) =  evalpoly(x, (1.5509733517804722f0, -0.4003010201031985f0, -0.07849861944294194f0, -0.034318853117591995f0, -0.0197180433173655f0, -0.01305950773199331f0, -0.009442372874146548f0, -0.007246728512402157f0, -0.00580742401295609f0, -0.004809187786009338f0))
        poly2(x) =  evalpoly(x, (1.5101218320928198f0, -0.41711633390586755f0, -0.09012382040477457f0, -0.04372994401908431f0, -0.027965493064761784f0, -0.020644781177568104f0, -0.016650786739707237f0, -0.01426196082884252f0, -0.012759847429264804f0, -0.011799303775587354f0, -0.011197445703074968f0))
        poly3(x) =  evalpoly(x, (1.4674622093394272f0, -0.43657629094633776f0, -0.10515555766694255f0, -0.05737184359324173f0, -0.04139162772734022f0, -0.03452772850528084f0, -0.031495443512532785f0, -0.030527000890325277f0, -0.0309169840192389f0, -0.03237139531475812f0, -0.03478996038640416f0))
        poly4(x) =  evalpoly(x, (1.4226911334908792f0, -0.4595135196210487f0, -0.12525053982206188f0, -0.07813854509440948f0, -0.06471427847205f0, -0.06208433913173031f0, -0.06519703281557247f0, -0.07279389536257878f0, -0.084959075171781f0, -0.102539850131046f0, -0.12705358515769605f0, -0.1607911206912746f0))
        poly5(x) =  evalpoly(x, (1.3754019718711163f0, -0.4872021832731848f0, -0.15331170134854022f0, -0.11184944491702783f0, -0.10884095252313576f0, -0.12295422312026907f0, -0.15221716396203505f0, -0.20049532364269734f0, -0.27617433306775174f0, -0.39351311430437586f0, -0.5757544060278792f0, -0.8605232357272398f0, -1.3088332057585401f0))
        poly6(x) =  evalpoly(x, (1.3250244979582302f0, -0.5217276475575667f0, -0.19490643048212622f0, -0.17162372682201127f0, -0.20275465292641914f0, -0.27879895311853475f0, -0.42069845728100574f0, -0.675948400853106f0, -1.1363431218392293f0, -1.9767211439543984f0, -3.5316967730957227f0, -6.446753640156048f0, -11.97703130208884f0))
        poly7(x) =  evalpoly(x, (1.2707074796501499f0, -0.5668391682878666f0, -0.2621607934324926f0, -0.2922441735330774f0, -0.4403978408504232f0, -0.7749476413813975f0, -1.498870837987561f0, -3.089708310445187f0, -6.6675959033810015f0, -14.89436036517319f0, -34.18120574251449f0, -80.15895841905397f0, -191.34894807629848f0, -463.5938853480342f0, -1137.38082216936f0))
        poly8(x) =  evalpoly(x, (1.2110560275684594f0, -0.6303064132874558f0, -0.38716640952066916f0, -0.5922782353119346f0, -1.23755558451305f0, -3.0320566617452474f0, -8.18168822157359f0, -23.55507217389693f0, -71.04099935893065f0, -221.879685319235f0, -712.1364793277636f0, -2336.1253314403966f0, -7801.945954775964f0, -26448.19586059192f0, -90799.48341621365f0, -315126.04064491636f0, -1104011.3443115912f0))
        poly9(x) =  evalpoly(x, (1.1613071521962828f0, -0.7011002845552895f0, -0.5805514744654373f0, -1.2436930610777865f0, -3.679383613496635f0, -12.815909243378957f0, -49.25672530759985f0, -202.18187354340904f0, -869.8602699308701f0, -3877.0058473132895f0, -17761.7071017094f0, -83182.69029154233f0, -396650.4505013548f0, -1920033.4136826345f0))
        poly10(x) = evalpoly(x, (1.1246173251197522f0, -0.7708450563609095f0, -0.8447940536449113f0, -2.4900973094503946f0, -10.239717411543843f0, -49.7490054655148f0, -267.09866751957054f0, -1532.66588382523f0, -9222.313478526092f0, -57502.51612140314f0, -368596.11674161063f0, -2415611.0887010912f0, -16120097.815816568f0, -109209938.52030899f0, -749380758.1942496f0, -5198725846.725541f0, -36409256888.1214f0))
        poly11(x) = evalpoly(x, (1.5910034537907922f0, 0.41600074399178694f0, 0.24579151426410342f0, 0.17948148291490615f0, 0.14455605708755515f0, 0.12320099331242772f0, 0.10893881157429353f0, 0.09885340987159291f0, 0.09143962920174975f0, 0.0858425915954139f0, 0.08154111871830322f0))
        poly12(x) = evalpoly(x, (1.5509733517804722f0, -0.4003010201031985f0, -0.07849861944294194f0, -0.034318853117591995f0, -0.0197180433173655f0, -0.01305950773199331f0, -0.009442372874146548f0, -0.007246728512402157f0, -0.00580742401295609f0, -0.004809187786009338f0))

        flag = false;
        kdm = 0f0;
        edm = 0f0;
        td = 0f0;
        km = 0f0;
        t = 0f0;
        x = 0f0;

        x = m;
        if  m < 0.0f0 
            x = m / ( m - 1.0f0 );
            flag::Bool = true;
        end
        x === 0.0f0 && return HALF_PI32
        x === 1.0f0 && return 1.0f0
        x > 1.0f0 && return NaN32

        if ( x < 0.1f0 ) 
            t = poly1( x - 0.05f0 );
        elseif ( x < 0.2f0 ) 
            t = poly2( x - 0.15f0 );
        elseif ( x < 0.3f0 ) 
            t = poly3( x - 0.25f0 );
        elseif ( x < 0.4f0 ) 
            t = poly4( x - 0.35f0 );
        elseif ( x < 0.5f0 ) 
            t = poly5( x - 0.45f0 );
        elseif ( x < 0.6f0 ) 
            t = poly6( x - 0.55f0 );
        elseif ( x < 0.7f0 ) 
            t = poly7( x - 0.65f0 );
        elseif ( x < 0.8f0 ) 
            t = poly8( x - 0.75f0 );
        elseif ( x < 0.85f0 ) 
            t = poly9( x - 0.825f0 );
        elseif ( x < 0.9f0 ) 
            t = poly10( x - 0.875f0 );
        else 
            td = 0.95f0 - x;
            kdm = poly11(td);
            edm = poly12(td);
            km = K( x );

            #// To avoid precision loss near 1f0, we use Eq. 33f0 f0rom f0ukushima (2009f0):
            t = ( HALF_PI32 + ( km * (kdm - edm) ) ) / kdm;
        end

        #// Complete the transformation mentioned above for m < 0:
        flag && return t * sqrt( 1.0f0 - m );

        return t;
    end
    m == 1.0f0 && return 1f0
    return NaN32
end

############################################################################
# f0irst Incomplete Elliptic Integral and Inverse Jacobi functions 
# (https://doi.org/10.1007/s00211f0-010f0-0321f0-8f0)
############################################################################
function serf(y::Float32, m::Float32)::Float32
    return 1.0f0 + y*(0.166667f0 + 0.166667f0*m + 
    y*(0.075f0 + (0.05f0 + 0.075f0*m)*m + 
       y*(0.0446429f0 + m*(0.0267857f0 + (0.0267857f0 + 0.0446429f0*m)*m) + 
          y*(0.0303819f0 + 
             m*(0.0173611f0 + 
                m*(0.015625f0 + (0.0173611f0 + 0.0303819f0*m)*m)) + 
             y*(0.0223722f0 + 
                m*(0.012429f0 + 
                   m*(0.0106534f0 + 
                    m*(0.0106534f0 + (0.012429f0 + 0.0223722f0*m)*m))) + 
                y*(0.0173528f0 + 
                   m*(0.00946514f0 + 
                    m*(0.00788762f0 + 
                    m*(0.00751202f0 + 
                    m*(0.00788762f0 + (0.00946514f0 + 0.0173528f0*m)*m)))) +
                    y*(0.0139648f0 + 
                    m*(0.00751953f0 + 
                    m*(0.00615234f0 + 
                    m*(0.00569661f0 + 
                    m*(0.00569661f0 + 
                    m*(0.00615234f0 + (0.00751953f0 + 
                    0.0139648f0*m)*m))))) + 
                    y*(0.0115518f0 + 
                    m*(0.00616096f0 + 
                    m*(0.00497616f0 + 
                    m*(0.00452378f0 + 
                    m*(0.00439812f0 + m*(0.00452378f0 + 
                    m*(0.00497616f0 + (0.00616096f0 + 
                    0.0115518f0*m)*m)))))) + (0.00976161f0 + 
                    m*(0.00516791f0 + 
                    m*(0.00413433f0 + 
                    m*(0.0037103f0 + m*(0.00354165f0 + 
                    m*(0.00354165f0 + m*(0.0037103f0 + 
                    m*(0.00413433f0 + (0.00516791f0 + 
                    0.00976161f0*m)*m))))))))*y))))))))
 
end

function asn(s::Float32, m::Float32)::Float32
    yA = 0.04094f0 - 0.00652f0*m
    y = s * s
    if y < yA
        return s*serf(y, m)
    end
    p = 1f0
    for _ in 1:10
        y = y / ((1f0+√(1f0-y))*(1f0+√(1f0-m*y)))
        p *= 2f0
        y < yA && return p*√y*serf(y, m)
    end
    return NaN32 
end

function acn(c::Float32, mc::Float32)::Float32
    m = 1f0 - mc
    x = c^2f0
    p = 1f0
    for _ in 1:10
        if (x > 0.5f0) 
            return p*asn(√(1f0-x), m)
        end
        d = √(mc + m * x)
        x = (√x + d)/(1f0+d)
        p *= 2f0
    end
    return NaN32
end

function _rawF(sinφ::Float32, m::Float32)
    yS = 0.90003085f0
    m == 0.0f0 && return asin(sinφ)
    m == 1.0f0 && return atanh(sinφ)
    
    sinφ^2f0 ≤ yS && return asn(sinφ, m)

    mc = 1f0 - m

    c = √(1f0-sinφ^2f0)
    x = c * c
    d2 = mc + m*x
    x <  yS*d2 && return (K(m) - asn(c/√(d2), m))

    v = mc*(1f0-x)
    v < x*d2 && return acn(c, m)
    return (K(m) - acn(√(v/d2), m))

end

function _F(φ::Float32, m::Float32)
    abs(φ) < HALF_PI32 && sign(φ)*return _rawF(sin(abs(φ)), m)
    j = round(φ/PI32)

    newφ = φ - j*PI32
    return 2f0*j*K(m) + sign(newφ)*_rawF(sin(mod(φ, PI32)), m)
end

function F(φ::Float32, m::Float32)
    if m > 1.f0
        ## Abramowitz & Stegum (17.4.15)
        m12 = sqrt(m)
        theta = asin(m12*sin(φ))
        signθ = sign(θ)
        absθ = abs(theta)
        return signθ/m12*_F(absθ, 1f0/m)
        #return NaN
    elseif m < 0f0
        # Abramowitz & Stegum (17.4.17)
        n = -m
        m12 = 1f0/sqrt(1f0+n)
        m1m = n/(1f0+n)
        newφ = HALF_PI32-φ
        signφ = sign(newφ)
        absφ = abs(newφ)
        return (m12*K(m1m) - signφ*m12*_F(absφ, m1m)) 
    end
    return _F(φ, m)
end

#Elliptic.jl Implementation
export am,
    sn, cn, dn, nn,
    sd, cd, dd, nd,
    sc, cc, dc, nc,
    ss, cs, ds, ns

# Abramowitz & Stegun, section 16.4, p571
#const _ambuf = Array{Float32}(undef, 10)
function _am(u::Float32, m::Float32, tol::Float32)::Float32
    _ambuf = @SVector Float32[0,0,0,0,0,0,0,0,0,0]#
    if u == 0.f0 return 0.f0 end

    sqrt_tol = sqrt(tol)
    if m < sqrt_tol
        # A&S 16.13.4
        return u - 0.25f0*m*(u - 0.5f0*sin(2.0f0*u))
    end
    m1 = 1.f0 - m
    if m1 < sqrt_tol
        # A&S 16.15.4
        t = tanh(u)
        return asin(t) + 0.25f0*m1*(t - u*(1.f0 - t^2f0))*cosh(u)
    end

    a,b,c,n = 1.f0, sqrt(m1), sqrt(m), 0
    while abs(c) > tol
        @assert n < 10
        a,b,c,n = 0.5f0*(a+b), sqrt(a*b), 0.5f0*(a-b), n+1
        @set! _ambuf[n] = c/a
        #_ambuf[n] = c/a

    end

    phi = ldexp(a*u, n)
    for i = n:-1:1
        phi = 0.5f0*(phi + asin(_ambuf[i]*sin(phi)))
    end
    return phi
end
_am(u::Float32, m::Float32)::Float32 = _am(u, m, eps(Float32))

"""
    am(u::Real, m::Real, [tol::Real=eps(Float32)])
Returns amplitude, φ, such that u = F(φ | m)
Landen sequence with convergence to `tol` used if `√(tol) ≤ m ≤ 1 - √(tol)`
"""
function am(u::Float32, m::Float32, tol::Float32)::Float32
    _am(u, m, tol)
end
am(u::Float32, m::Float32)::Float32 = am(u, m, eps(Float32))

function sn(u::Float32, m::Float32)::Float32
    # Abramowitz & Stegun, section 16.10, p573
    lt0 = m < 0.f0
    gt1 = m > 1.f0
    if !(lt0 || gt1)
        phi = _am(u,m)
        return sin(phi)
    elseif lt0
        mu1 = 1.0f0/(1.f0 - m)
        mu = -m*mu1
        sqrtmu1 = sqrt(mu1)
        v = u/sqrtmu1
        phi = _am(v,mu)
        s = sin(phi)
        return (sqrtmu1*s)/sqrt(1.f0 - mu*s^2f0)
    elseif gt1
        mu = 1f0/m
        v = u*sqrt(m)
        phi = _am(v,mu)
        return (sqrt(mu)*sin(phi))
    end
    return Inf32
end

function cn(u::Float32, m::Float32)::Float32
    # Abramowitz & Stegun, section 16.10, p573
    lt0 = m < 0.f0
    gt1 = m > 1.f0
    if !(lt0 || gt1)
        phi = _am(u,m)
        return cos(phi)
    elseif lt0
        mu1 = 1.0f0/(1.f0 - m)
        mu = -m*mu1
        sqrtmu1 = sqrt(mu1)
        v = u/sqrtmu1
        phi = _am(v,mu)
        s = sin(phi)
        return (cos(phi))/sqrt(1.f0 - mu*s^2f0)
    elseif gt1
        mu = 1f0/m
        v = 2.0f0#u*sqrt(m)
        phi = _am(v,mu)
        return (sqrt(1.f0 - mu*sin(phi)^2f0))
    end
    return Inf32
end

function dn(u::Float32, m::Float32)::Float32
    # Abramowitz & Stegun, section 16.10, p573
    lt0 = m < 0.f0
    gt1 = m > 1.f0
    if !(lt0 || gt1)
        phi = _am(u,m)
        return (sqrt(1.f0 - m*sin(phi)^2f0))
    elseif lt0
        mu1 = 1.0f0/(1.f0 - m)
        mu = -m*mu1
        sqrtmu1 = sqrt(mu1)
        v = u/sqrtmu1
        phi = _am(v,mu)
        s = sin(phi)
        return (1.0f0)/sqrt(1.f0 - mu*s^2f0)
    elseif gt1
        mu = 1f0/m
        v = 2.0f0#u*sqrt(m)
        phi = _am(v,mu)
        return cos(phi)
    end
    return Inf32
end

xn = ((:s,:(sn(u,m))), (:c,:(cn(u,m))), (:d,:(dn(u,m))), (:n,:(1.f0)))
for (p,num) in xn, (q,den) in xn
    f = Symbol(p, q)
    if (p == q)
        @eval ($f)(::Float32, ::Float32)::Float32 = 1.0f0
    elseif (q != :n)
        @eval ($f)(u::Float32, m::Float32)::Float32 = ($num)/($den)
    end
end

##https://doi-org.ezp-prod1.hul.harvard.edu/10.031007/s10569f0-008f0-9177f0-y
##https://link.springer.com/article/10.1007/s10569f0-008f0-9177f0-y
##TODO: Implement Δsn algorithm for half0 angle transformation
#
#_Kscreen(m::Float32)::Float32 = HALF_PI32*(1.0f0 + m*(0.25f0 + m*(0.36f0 + m*(0.09765625f0 + m*0.07476806640625f0))))
#
#function _rawΔSN(up::Float32, m::Float32, Kscreen::Float32, Kactual::Float32, kp::Float32)
#    if up > 0.031f0
#        u = up/2f0
#        _sn = _rawSN(up, m, Kscreen, Kactual, kp)
#        _sn2 = _sn*_sn
#        _sn4 = _sn2*_sn2
#        return 2f0/(1f0-m*_sn4)*(_rawCN(u, m, Kscreen, Kactual, kp)*_rawDN(u, m, Kscreen, Kactual, kp)*_rawΔSN(u, m, Kscreen, Kactual, kp)+u*(_rawΔCN(u, m, Kscreen, Kactual, kp)+_rawΔDN(u, m, Kscreen, Kactual, kp)-_rawΔCN(u, m, Kscreen, Kactual, kp)*_rawΔDN(u, m, Kscreen, Kactual, kp)-m*_sn4))
#    else
#        return up - _rawSN(up, m, Kscreen, Kactual, kp)
#    end
#end
#
#function _ΔSN(u::Float32, m)
#    return _rawΔSN(u, m, K(m), K(m), √(1f0-m))
#end
#
#function _rawΔCN(up::Float32, m::Float32, Kscreen::Float32, Kactual::Float32, kp::Float32)
#    if up > 0.031f0
#       u = up/2f0
#       _sn = _rawSN(u,m,Kscreen,Kactual, kp)
#       _sn2 = _sn*_sn
#       return ((1f0+_rawCN(u,m,Kscreen,Kactual, kp))*_rawΔCN(u,m,Kscreen,Kactual, kp)+(1f0-2f0*m*_sn2)*_sn2)/(1f0-m*_sn2*_sn2)
#    else 
#        return 1f0 - _rawCN(up, m, Kscreen, Kactual, kp)
#    end
#end
#
#function _ΔCN(u::Float32, m)
#    return _rawΔCN(u, m, K(m), K(m), √(1f0-m))
#end
#
#function _rawΔDN(up::Float32, m::Float32, Kscreen::Float32, Kactual::Float32, kp::Float32) 
#    if up > 0.031f0
#       u = up/2f0 
#       _sn = _rawSN(u,m,Kscreen,Kactual, kp)
#       _sn2 = _sn*_sn
#
#       return ((1f0+_rawDN(u,m,Kscreen,Kactual, kp))*_rawΔDN(u,m,Kscreen,Kactual, kp)+m*(1f0-2*_sn2)*_sn2)/(1f0-m*_sn2*_sn2)
#    else 
#        return 1f0 - _rawCN(up, m, Kscreen, Kactual, kp)
#    end
#end
#
#function _ΔDN(u::Float32, m)
#    return _rawΔDN(u, m, K(m), K(m), √(1f0-m))
#end
#
#function _XNloop(u::Float32, m::Float32, n::Float32)
#    up = u / (2.0f0^Float32(n))
#    up2 = up^2f0
#    up4 = up2 * up2
#    sn = up*(up2*(up2*((m*((-(m/5040f0)-3f0/112f0)*m-3f0/112f0)-1f0/5040f0)*up2+(m/120f0+7f0/60f0)*m+1f0/120f0)-m/6f0-1f0/6f0)+1f0)
#    cn = 1f0 + up2*(-(1f0/2f0)+up2*(1f0/24f0+m/6f0+up2*(-(1f0/720f0)+(-(11f0/180f0)-m/45f0)*m+(-(1f0/40320f0)+m*(-(17f0/1680f0)+(-(19f0/840f0)-m/630f0)*m))*up2)))
#    dn = m*(m*(up4*(1f0/24f0-(11f0*up2)/180f0)-(m*up4*up2)/720f0) + (up2*(1f0/6f0-up2/45f0)-1f0/2f0)*up2) + 1f0 
#
#    
#    i = 0
#    while i < n
#        sn2 = sn * sn
#        cn2 = cn * cn
#        dn2 = dn * dn
#        sn4 = sn2 * sn2
#
#        den = 1f0/(1.0f0-m*sn4)
#        sn = 2.0f0*(sn*cn*dn)*den
#        cn = (cn2-sn2*dn2)*den
#        dn = (dn2-m*sn2*cn2)*den
#        i += 1
#    end
#    return sn, cn, dn
#end
#
#function _ΔXNloop(u::Float32, m::Float32, n::Float32) #https://doi-org.ezp-prod1.hul.harvard.edu/10.1007/s10569-009-9228-z
#    up = 0
#
#    up = u / (2.0f0^n)
#    up2 = up^2f0
#    up4 = up2 * up2
#    sn = up*(up2*(up2*((m*((-(m/5040f0)-3f0/112f0)*m-3f0/112f0)-1f0/5040f0)*up2+(m/120f0+7f0/60f0)*m+1f0/120f0)-m/6f0-1f0/6f0)+1f0)
#    cn = 1f0 + up2*(-(1f0/2f0)+up2*(1f0/24f0+m/6f0+up2*(-(1f0/720f0)+(-(11f0/180f0)-m/45f0)*m+(-(1f0/40320f0)+m*(-(17f0/1680f0)+(-(19f0/840f0)-m/630f0)*m))*up2)))
#    dn = m*(m*(up4*(1f0/24f0-(11f0*up2)/180f0)-(m*up4*up2)/720f0) + (up2*(1f0/6f0-up2/45f0)-1f0/2f0)*up2) + 1f0 
#    sn = up*(up2*(up2*((m/120f0+7f0/60f0)*m+1f0/120f0)-m/6f0-1f0/6f0)+1f0)
#    cn = 1f0 + up2*(-(1f0/2f0)+up2*(1f0/24f0+m/6f0+up2*(-(1f0/720f0)+(-(11f0/180f0)-m/45f0)*m)))
#    dn = m*(m*(up4*(1f0/24f0-(11f0*up2)/180f0)-(m*up4*up2)/720f0) - (1f0/2f0)*up2) + 1f0 
#
#    Δsn = up-sn
#    Δcn = 1f0-cn
#    Δdn = 1f0-dn
#
#
#    i = 0
#    while i < n
#        sn2 = sn * sn
#        sn4 = sn2 * sn2
#
#        den = 1f0/(1.0f0-m*sn4)
#        Δsn = 2.0f0*(Δsn*cn*dn + up*(Δcn*(1f0 - Δdn) + Δdn  - m*sn4))*den
#        Δcn = ((1f0+cn)*Δcn + sn2-2f0*m*sn4)*den
#        Δdn = ((1f0+dn)*Δdn + m*(sn2-2f0*sn4))*den
#
#        up *= 2f0
#        sn = up - Δsn
#        cn = 1f0 - Δcn
#        dn = 1f0 - Δdn
#        
#        i += 1f0
#    end
#    return sn, cn, dn
#end
#
#function fold_0_25(u1::Float32, m::Float32, kp::Float32) 
#    u1 == 0f0 && return 0f0, 1f0, 1f0
#    #return (u1 > 0f0 ? max(6f0 + floor(log2(u1)), 1f0) : 0f0), (u1 > 0f0 ? max(6f0 + floor(log2(u1)), 1f0) : 0f0), (u1 > 0f0 ? max(6f0 + floor(log2(u1)), 1f0) : 0f0)
#    sn, cn, dn = _ΔXNloop(u1, m, (u1 > 0f0 ? max(6f0 + floor(log2(u1)), 1f0) : 0f0))
#    den = 1.0f0/(1.0f0+kp-m*sn*sn)
#    #return sn, cn, dn
#    return den*√(1f0+kp)*(cn*dn-kp*sn), den*√(kp*(1f0+kp))*(cn+sn*dn), den*√kp*((1f0+kp)*dn+m*sn*cn)
#end
#
#function fold_0_50(u1::Float32, m::Float32, Kscreen::Float32, Kactual::Float32, kp::Float32)
#    u1 == 0f0 && return 0f0, 1f0, 1f0
#
#    if u1 > 0.25f0*Kscreen 
#        sn, cn, dn = fold_0_25(Kactual/2f0 - u1, m, kp)
#    else
#        sn, cn, dn =  _ΔXNloop(u1, m, u1 > 0f0 ?  max(6f0 + floor(log2(u1)), 1f0) : 0f0)
#    end
#    return cn/dn, kp*sn/dn, kp/dn
#end
#
#function fold_1_00(u1::Float32, m::Float32, Kscreen::Float32, Kactual::Float32, kp::Float32)
#    u1 == 0f0 && return 0f0, 1f0, 1f0
#
#    if u1 > 0.5f0Kscreen
#        sn, cn, dn = fold_0_50(Kactual - u1, m, Kscreen, Kactual, kp)
#    elseif u1 > 0.25f0Kscreen 
#        sn, cn, dn = fold_0_25(Kactual/2f0 - u1, m, kp)
#    else
#        sn, cn, dn = _ΔXNloop(u1, m, u1 > 0f0 ?  max(6f0+floor(log2(u1)), 1f0) : 0f0)
#    end
#    return cn/dn, -kp*sn/dn, kp/dn
#end
#
#funcs = ((:_SN, :_rawSN), (:_CN, :_rawCN), (:_DN, :_rawDN))
#for (enum, funcpair) in enumerate(funcs)
#    (func, helper) = funcpair
#    @eval begin
#        function $(helper)(u::Float32, m::Float32, Kscreen::Float32, Kactual::Float32, kp::Float32) 
#            u = u > 4f0*Kscreen ?  u % 4f0*Kactual : u
#            u = u > 2f0*Kscreen ? u - 2f0*Kactual : u
#            u > Kscreen && return fold_1_00(u - Kactual, m, Kscreen, Kactual, kp)[$enum]
#            u > 0.5f0*Kscreen && return fold_0_50(Kactual - u, m, Kscreen, Kactual, kp)[$enum]
#            u > 0.25f0*Kscreen && return fold_0_25(Kactual/2f0 - u, m, kp)[$enum]
#            return _ΔXNloop(u, m, u > 0f0 ? max(6f0+log2(u), 1f0) : 0f0)[$enum] 
#        end
#
#        function $(func)(u::Float32, m::Float32)
#             return $(helper)(u, m, K(m), K(m), √(1f0-m)) 
#        end
#    end
#end
#
#_sc_helper(jacobituple::Tuple{Float32, Float32}) = jacobituple[1]/jacobituple[2]
#function _rawSC(u::Float32, m::Float32, Kscreen::Float32, Kactual::Float32, kp::Float32) 
#    u = u > 4f0*Kscreen ?  u % 4f0*Kactual : u
#    u = u > 2f0*Kscreen ? u - 2f0*Kactual : u
#    u > Kscreen && return _sc_helper(fold_1_00(u - Kactual, m, Kscreen, Kactual, kp))
#    u > 0.5f0Kscreen && return _sc_helper(fold_0_50(Kactual - u, m, Kscreen, Kactual, kp))
#    u > 0.25f0Kscreen && return _sc_helper(fold_0_25(Kactual/2f0 - u, m, kp))
#    return _sc_helper(_ΔXNloop(u, m, u > 0f0 ? max(6f0+floor(log2(u)), 1f0) : 0f0))
#end
#_SC(u::Float32, m::Float32) = _rawSC(u, m, K(m), K(m), √(1f0-m))
#
#_sd_helper(jacobituple::Tuple{Float32, Float32}) = jacobituple[1]/jacobituple[3]
#function _rawSD(u::Float32, m::Float32, Kscreen::Float32, Kactual::Float32, kp::Float32) 
#    u = u > 4Kscreen ?  u % 4Kactual : u
#    u = u > 2Kscreen ? u - 2Kactual : u
#    u > Kscreen && return _sd_helper(fold_1_00(u - Kactual, m, Kscreen, Kactual, kp))
#    u > 0.5f0Kscreen && return _sd_helper(fold_0_50(Kactual - u, m, Kscreen, Kactual, kp))
#    u > 0.25f0Kscreen && return _sd_helper(fold_0_25(Kactual/2f0 - u, m, kp))
#    return _sd_helper(_ΔXNloop(u, m, u > 0f0 ? max(6f0+floor(log2(u)), 1f0) : 0f0))
#end
#
#_SD(u::Float32, m::Float32) = _rawSD(u, m, K(m), K(m), √(1f0-m))
#
#
#function sn(u::Float32, m::Float32)  
#    signu = sign(u)
#    u = abs(u)
#    #m < 1f0 && return signu*(_SN(u, m))
#    m < 1f0 && return signu*(u - _ΔSN(u, m))
#    sqrtm = √m
#    #return signu*_SN(u*sqrtm, 1f0/m)/sqrtm
#    return signu*(u*sqrtm -_ΔSN(u*sqrtm, 1f0/m))/sqrtm
#end
#
#function cn(u::Float32, m::Float32)  
#    u = abs(u)
#    m < 1f0 && return _CN(u, m)
#    sqrtm = √m
#    return _DN(u*sqrtm, 1/m)
#end
#
#function dn(u::Float32, m::Float32)  
#    u = abs(u)
#    m < 1f0 && return _DN(u, m)
#    sqrtm = √m
#    return _CN(u*sqrtm, 1/m)
#end
#
#function sc(u::Float32, m::Float32)
#    signu = sign(u)
#    u = abs(u)
#
#    m < 1f0 && return signu*_SC(u, m)
#    sqrtm = √m
#    return signu*_SD(u*sqrtm, 1/m)/sqrtm
#
#end
#
#function sd(u::Float32, m::Float32)
#    signu = sign(u)
#    u = abs(u)
#
#    m < 1f0 && return signu*_SD(u, m)
#    sqrtm = √m
#    return signu*_SC(u*sqrtm, 1/m)/sqrtm
#
#end
#
#