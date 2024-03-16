module FukushimaAlg
using StaticArrays, Setfield

export K, E, F, Pi, J
export sn, cn, dn, sc, sd, asn, acn

const ONE_DIV_PI = inv(π)

############################################################################
# Complete Elliptic Integrals (https://doi.org/10.1007/s10569-009-9228-z)
############################################################################
"""
``K(m) = \\int_0^{\\pi/2}\\frac{d\\theta}{\\sqrt{1-k^2\\sin(\\theta)^2}}.``

Returns the complete elliptic integral of the first kind.
    
# Arguments

- `m` : Elliptic modulus
"""
function K(m::T) where T
    if m < one(T)
        # I didn't really see any speedup from evalpoly, so I left the evaluation in this form
        poly1(x::T)  =  T(1.59100345379079220) + (x * (T(0.41600074399178694) + (x * (T(0.24579151426410342) + (x * (T(0.17948148291490615) + (x * (T(0.14455605708755515) + (x * (T(0.123200993312427720) + (x * (T(0.108938811574293530) + (x * (T(0.098853409871592910) + (x * (T(0.091439629201749750) + (x * (T(0.08584259159541390) + (x * (T(0.08154111871830322)))))))))))))))))))))
        poly2(x::T)  =  T(1.63525673226458000) + (x * (T(0.47119062614873230) + (x * (T(0.30972841083149960) + (x * (T(0.25220831177313570) + (x * (T(0.22672562321968465) + (x * (T(0.215774446729585980) + (x * (T(0.213108771877348920) + (x * (T(0.216029124605188280) + (x * (T(0.223255831633057900) + (x * (T(0.23418050129420992) + (x * (T(0.24855768297226408) + (x * (T(0.266363809892617540)))))))))))))))))))))))
        poly3(x::T)  =  T(1.68575035481259600) + (x * (T(0.54173184861328030) + (x * (T(0.40152443839069024) + (x * (T(0.36964247342088910) + (x * (T(0.37606071535458363) + (x * (T(0.405235887085125900) + (x * (T(0.453294381753999050) + (x * (T(0.520518947651184200) + (x * (T(0.609426039204995000) + (x * (T(0.72426352228290890) + (x * (T(0.87101384770981240) + (x * (T(1.057652872753547000)))))))))))))))))))))))
        poly4(x::T)  =  T(1.74435059722561330) + (x * (T(0.63486427537193530) + (x * (T(0.53984256416444550) + (x * (T(0.57189270519378740) + (x * (T(0.67029513626540620) + (x * (T(0.832586590010977200) + (x * (T(1.073857448247933300) + (x * (T(1.422091460675497700) + (x * (T(1.920387183402304700) + (x * (T(2.63255254833165430) + (x * (T(3.65210974731903940) + (x * (T(5.115867135558866000) + (x * (T(7.224080007363877000)))))))))))))))))))))))))
        poly5(x::T)  =  T(1.81388393681698260) + (x * (T(0.76316324570055730) + (x * (T(0.76192860532159580) + (x * (T(0.95107465366842790) + (x * (T(1.31518067170316100) + (x * (T(1.928560693477410900) + (x * (T(2.937509342531378700) + (x * (T(4.594894405442878000) + (x * (T(7.330071221881720000) + (x * (T(11.8715125974253010) + (x * (T(19.4585137482293800) + (x * (T(32.20638657246427000) + (x * (T(53.73749198700555000) + (x * (T(90.27388602941000000)))))))))))))))))))))))))))
        poly6(x::T)  =  T(1.89892491027155350) + (x * (T(0.95052179461824450) + (x * (T(1.15107758995901580) + (x * (T(1.75023910698630060) + (x * (T(2.95267681263687500) + (x * (T(5.285800396121451000) + (x * (T(9.832485716659980000) + (x * (T(18.78714868327559600) + (x * (T(36.61468615273698000) + (x * (T(72.4529239512777100) + (x * (T(145.107957734706900) + (x * (T(293.4786396308497000) + (x * (T(598.3851815055010000) + (x * (T(1228.420013075863400) + (x * (T(2536.529755382764500)))))))))))))))))))))))))))))
        poly7(x::T)  =  T(2.00759839842437640) + (x * (T(1.24845723121234740) + (x * (T(1.92623465707647970) + (x * (T(3.75128964008758770) + (x * (T(8.11994455493204500) + (x * (T(18.66572130873555200) + (x * (T(44.60392484291437400) + (x * (T(109.5092054309498300) + (x * (T(274.2779548232414000) + (x * (T(697.559800860632700) + (x * (T(1795.71601450024720) + (x * (T(4668.381716790390000) + (x * (T(12235.76246813664300) + (x * (T(32290.17809718321000) + (x * (T(85713.07608195965000) + (x * (T(228672.1890493117) + (x * (T(612757.27119158520)))))))))))))))))))))))))))))))))
        poly8(x::T)  =  T(2.15651564749964340) + (x * (T(1.79180564184946320) + (x * (T(3.82675128746571320) + (x * (T(10.3867246836379720) + (x * (T(31.4033140546807030) + (x * (T(100.9237039498695500) + (x * (T(337.3268282632273000) + (x * (T(1158.707930567827800) + (x * (T(4060.990742193632300) + (x * (T(14454.0018403434480) + (x * (T(52076.6610759940450) + (x * (T(189493.6591462156800) + (x * (T(695184.5762413896000) + (x * (T(2567994.048255285000) + (x * (T(9541921.966748387000) + (x * (T(35634927.44218076) + (x * (T(133669298.46120408) + (x * (T(503352186.68662846) + (x * (T(1901975729.538660) + (x * T(7208915015.33010400))))))))))))))))))))))))))))))))))))))
        poly9(x::T)  =  T(2.31812262171251060) + (x * (T(2.61692015029123270) + (x * (T(7.89793507573135600) + (x * (T(30.5023971544667240) + (x * (T(131.486936552352860) + (x * (T(602.9847637356492000) + (x * (T(2877.024617809973000) + (x * (T(14110.51991915180400) + (x * (T(70621.44088156540000) + (x * (T(358977.266582531000) + (x * (T(1847238.26372397180) + (x * (T(9600515.416049214000) + (x * (T(50307677.08502367000) + (x * (T(265444188.6527128000) + (x * (T(1408862325.028702700) + (x * (T(7515687935.373775)))))))))))))))))))))))))))))))
        poly10(x::T) =  T(2.47359617375134400) + (x * (T(3.72762424411809900) + (x * (T(15.6073930355493060) + (x * (T(84.1285084280588800) + (x * (T(506.981819704061370) + (x * (T(3252.277058145123600) + (x * (T(21713.24241957434400) + (x * (T(149037.0451890932700) + (x * (T(1043999.331089990800) + (x * (T(7427974.81704203900) + (x * (T(53503839.6755866100) + (x * (T(389249886.9948708400) + (x * (T(2855288351.100810500) + (x * (T(21090077038.76684000) + (x * (T(156699833947.7902000) + (x * (T(1170222242422.440) + (x * (T(8777948323668.9375) + (x * (T(66101242752484.950) + (x * (T(499488053713388.8) + (x * T(37859743397240296.0))))))))))))))))))))))))))))))))))))))
        poly11(x::T) =                           (x * (T(0.06250000000000000) + (x * (T(0.03125000000000000) + (x * (T(0.02050781250000000) + (x * (T(0.01513671875000000) + (x * (T(0.011934280395507812) + (x * (T(0.009816169738769531) + (x * (T(0.008315593004226685) + (x * (T(0.007199153304100037) + (x * (T(0.00633745662344154) + (x * (T(0.00565311038371874) + (x * (T(0.005097046040418718) + (x * (T(0.004636680381850056) + (x * (T(0.004249547423822886) + (x * (T(0.003919665602267974)))))))))))))))))))))))))))))
        poly12(x::T) =  T(1.59100345379079220) + (x * (T(0.41600074399178694) + (x * (T(0.24579151426410342) + (x * (T(0.17948148291490615) + (x * (T(0.14455605708755515) + (x * (T(0.123200993312427720) + (x * (T(0.108938811574293530) + (x * (T(0.098853409871592910) + (x * (T(0.091439629201749750) + (x * (T(0.08584259159541390) + (x * (T(0.08154111871830322)))))))))))))))))))))

        flag = false
        kdm = zero(T)
        td = zero(T)
        qd = zero(T)
        t = zero(T)

        x = m
        if m < zero(T)
            x = m / ( m - one(T) );
            flag = true;
        end
        x == zero(T) && return T(π/2)
        x == one(T) && return T(Inf)
        x > one(T) && return T(NaN)
        
        if x < T(0.1)
            t = poly1( x - T(0.05));
        elseif ( x < T(0.2)) 
            t = poly2( x - T(0.15));
        elseif ( x < T(0.3)) 
            t = poly3( x - T(0.25));
        elseif ( x < T(0.4)) 
            t = poly4( x - T(0.35));
        elseif ( x < T(0.5)) 
            t = poly5( x - T(0.45));
        elseif ( x < T(0.6)) 
            t = poly6( x - T(0.55));
        elseif ( x < T(0.7)) 
            t = poly7( x - T(0.65));
        elseif ( x < T(0.8)) 
            t = poly8( x - T(0.75));
        elseif ( x < T(0.85))
            t = poly9( x - T(0.825));
        elseif ( x < T(0.9)) 
            t = poly10( x - T(0.875));
        else 
            td = one(T) - x;
            qd = poly11( td );
            kdm = poly12(td - T(0.05));
            t = -log( qd ) * ( kdm * T(ONE_DIV_PI));
        end
        # Complete the transformation mentioned above for m < 0:
        flag && return t / sqrt( one(T) - m );
        
        return t
    end
    m == one(T) && return T(Inf)
    return T(NaN)
end

"""
``E(m) = \\int_0^{\\pi/2}\\sqrt{1-k^2\\sin(\\theta)^2}d\\theta.``

Returns the complete elliptic integral of the second kind.

# Arguments

- `m` : Elliptic modulus
"""
function E(m::T) where T
    m == zero(T) && return T(π/2)
    if m < one(T)
        # I didn't really see any speedup from evalpoly, so I left the evaluation in this form
        poly1(x::T) =  T(1.5509733517804722) + (x * (T(-0.40030102010319850) + (x * (T(-0.07849861944294194) + (x * (T(-0.034318853117591995) + (x * (T(-0.019718043317365500) + (x * (T(-0.013059507731993310) + (x * (T(-0.009442372874146548) + (x * (T(-0.007246728512402157) + (x * (T(-0.005807424012956090) + (x * (T(-0.004809187786009338)))))))))))))))))))
        poly2(x::T) =  T(1.5101218320928198) + (x * (T(-0.41711633390586755) + (x * (T(-0.09012382040477457) + (x * (T(-0.043729944019084310) + (x * (T(-0.027965493064761784) + (x * (T(-0.020644781177568104) + (x * (T(-0.016650786739707237) + (x * (T(-0.014261960828842520) + (x * (T(-0.012759847429264804) + (x * (T(-0.011799303775587354) + (x * (T(-0.011197445703074968)))))))))))))))))))))
        poly3(x::T) =  T(1.4674622093394272) + (x * (T(-0.43657629094633776) + (x * (T(-0.10515555766694255) + (x * (T(-0.057371843593241730) + (x * (T(-0.041391627727340220) + (x * (T(-0.034527728505280840) + (x * (T(-0.031495443512532785) + (x * (T(-0.030527000890325277) + (x * (T(-0.030916984019238900) + (x * (T(-0.032371395314758120) + (x * (T(-0.034789960386404160)))))))))))))))))))))
        poly4(x::T) =  T(1.4226911334908792) + (x * (T(-0.45951351962104870) + (x * (T(-0.12525053982206188) + (x * (T(-0.078138545094409480) + (x * (T(-0.064714278472050000) + (x * (T(-0.062084339131730310) + (x * (T(-0.065197032815572470) + (x * (T(-0.072793895362578780) + (x * (T(-0.084959075171781000) + (x * (T(-0.102539850131046000) + (x * (T(-0.127053585157696050) + (x * (T(-0.1607911206912746)))))))))))))))))))))))
        poly5(x::T) =  T(1.3754019718711163) + (x * (T(-0.48720218327318480) + (x * (T(-0.15331170134854022) + (x * (T(-0.111849444917027830) + (x * (T(-0.108840952523135760) + (x * (T(-0.122954223120269070) + (x * (T(-0.152217163962035050) + (x * (T(-0.200495323642697340) + (x * (T(-0.276174333067751740) + (x * (T(-0.393513114304375860) + (x * (T(-0.575754406027879200) + (x * (T(-0.8605232357272398) + (x * (T(-1.3088332057585401)))))))))))))))))))))))))
        poly6(x::T) =  T(1.3250244979582302) + (x * (T(-0.52172764755756670) + (x * (T(-0.19490643048212622) + (x * (T(-0.171623726822011270) + (x * (T(-0.202754652926419140) + (x * (T(-0.278798953118534750) + (x * (T(-0.420698457281005740) + (x * (T(-0.675948400853106000) + (x * (T(-1.136343121839229300) + (x * (T(-1.976721143954398400) + (x * (T(-3.531696773095722700) + (x * (T(-6.4467536401560480) + (x * (T(-11.977031302088840)))))))))))))))))))))))))
        poly7(x::T) =  T(1.2707074796501499) + (x * (T(-0.56683916828786660) + (x * (T(-0.26216079343249260) + (x * (T(-0.292244173533077400) + (x * (T(-0.440397840850423200) + (x * (T(-0.774947641381397500) + (x * (T(-1.498870837987561000) + (x * (T(-3.089708310445187000) + (x * (T(-6.667595903381001500) + (x * (T(-14.89436036517319000) + (x * (T(-34.18120574251449000) + (x * (T(-80.158958419053970) + (x * (T(-191.34894807629848) + (x * (T(-463.59388534803420) + (x * (T(-1137.380822169360)))))))))))))))))))))))))))))
        poly8(x::T) =  T(1.2110560275684594) + (x * (T(-0.63030641328745580) + (x * (T(-0.38716640952066916) + (x * (T(-0.592278235311934600) + (x * (T(-1.237555584513050000) + (x * (T(-3.032056661745247400) + (x * (T(-8.181688221573590000) + (x * (T(-23.55507217389693000) + (x * (T(-71.04099935893065000) + (x * (T(-221.8796853192350000) + (x * (T(-712.1364793277636000) + (x * (T(-2336.1253314403966) + (x * (T(-7801.9459547759640) + (x * (T(-26448.195860591920) + (x * (T(-90799.48341621365) + (x * (T(-315126.04064491636) + (x * T(-1104011.3443115912))))))))))))))))))))))))))))))))
        poly9(x::T) =  T(1.1613071521962828) + (x * (T(-0.70110028455528950) + (x * (T(-0.58055147446543730) + (x * (T(-1.243693061077786500) + (x * (T(-3.679383613496635000) + (x * (T(-12.81590924337895700) + (x * (T(-49.25672530759985000) + (x * (T(-202.1818735434090400) + (x * (T(-869.8602699308701000) + (x * (T(-3877.005847313289500) + (x * (T(-17761.70710170940000) + (x * (T(-83182.690291542330) + (x * (T(-396650.45050135480) + (x * T(-1920033.41368263450))))))))))))))))))))))))))
        poly10(x::T) = T(1.1246173251197522) + (x * (T(-0.77084505636090950) + (x * (T(-0.84479405364491130) + (x * (T(-2.490097309450394600) + (x * (T(-10.23971741154384300) + (x * (T(-49.74900546551480000) + (x * (T(-267.0986675195705400) + (x * (T(-1532.665883825230000) + (x * (T(-9222.313478526092000) + (x * (T(-57502.51612140314000) + (x * (T(-368596.1167416106300) + (x * (T(-2415611.0887010912) + (x * (T(-16120097.815816568) + (x * (T(-109209938.52030899) + (x * (T(-749380758.1942496) + (x * (T(-5198725846.7255410) + (x * T(-36409256888.121400))))))))))))))))))))))))))))))))
        poly11(x::T) = T(1.5910034537907922) + (x * (T(0.416000743991786940) + (x * (T(0.245791514264103420) + (x * (T(0.1794814829149061500) + (x * (T(0.1445560570875551500) + (x * (T(0.1232009933124277200) + (x * (T(0.1089388115742935300) + (x * (T(0.0988534098715929100) + (x * (T(0.0914396292017497500) + (x * (T(0.0858425915954139000) + (x * (T(0.0815411187183032200)))))))))))))))))))))
        poly12(x::T) = T(1.5509733517804722) + (x * (T(-0.40030102010319850) + (x * (T(-0.07849861944294194) + (x * (T(-0.034318853117591995) + (x * (T(-0.019718043317365500) + (x * (T(-0.013059507731993310) + (x * (T(-0.009442372874146548) + (x * (T(-0.007246728512402157) + (x * (T(-0.005807424012956090) + (x * (T(-0.004809187786009338)))))))))))))))))))

        flag = false;
        kdm = zero(T);
        edm = zero(T);
        td  = zero(T);
        km  = zero(T);
        t   = zero(T);
        x   = zero(T);

        x = m;
        if  m < zero(T)
            x = m / ( m - one(T) );
            flag = true;
        end
        x === zero(T) && return T(π/2)
        x === one(T) && return one(T)
        x > one(T) && return T(NaN)

        if ( x < T(0.1)) 
            t = poly1( x - T(0.05));
        elseif ( x < T(0.2)) 
            t = poly2( x - T(0.15));
        elseif ( x < T(0.3)) 
            t = poly3( x - T(0.25));
        elseif ( x < T(0.4)) 
            t = poly4( x - T(0.35));
        elseif ( x < T(0.5)) 
            t = poly5( x - T(0.45));
        elseif ( x < T(0.6)) 
            t = poly6( x - T(0.55));
        elseif ( x < T(0.7)) 
            t = poly7( x - T(0.65));
        elseif ( x < T(0.8)) 
            t = poly8( x - T(0.75));
        elseif ( x < T(0.85)) 
            t = poly9( x - T(0.825));
        elseif ( x < T(0.9)) 
            t = poly10( x - T(0.875));
        else 
            td = T(0.95) - x;
            kdm = poly11(td);
            edm = poly12(td);
            km = K( x );

            #// To avoid precision loss near 1, we use Eq. 33 from Fukushima (2009):
            t = ( T(π/2) + ( km * (kdm - edm) ) ) / kdm;
        end

        #// Complete the transformation mentioned above for m < 0:
        flag && return t * sqrt( one(T) - m );

        return t;
    end
    m == one(T) && return one(T)
    return T(NaN)
end

############################################################################
# First Incomplete Elliptic integral and Inverse Jacobi Functions 
# (https://doi.org/T(10).1007/s00211-010-0321-8)
############################################################################
function serf(y::A, m::B) where {A,B}
    T = promote_type(A,B)
    return 1 + y*(T(0.166667) + T(0.166667)*m + 
    y*(T(0.075) + (T(0.05) + T(0.075)*m)*m + 
       y*(T(0.0446429) + m*(T(0.0267857) + (T(0.0267857) + T(0.0446429)*m)*m) + 
          y*(T(0.0303819) + 
             m*(T(0.0173611) + 
                m*(T(0.015625) + (T(0.0173611) + T(0.0303819)*m)*m)) + 
             y*(T(0.0223722) + 
                m*(T(0.012429) + 
                   m*(T(0.0106534) + 
                    m*(T(0.0106534) + (T(0.012429) + T(0.0223722)*m)*m))) + 
                y*(T(0.0173528) + 
                    m*(T(0.00946514) + 
                    m*(T(0.00788762) + 
                    m*(T(0.00751202) + 
                    m*(T(0.00788762) + (T(0.00946514) + T(0.0173528)*m)*m)))) +
                    y*(T(0.01396480) + 
                    m*(T(0.00751953) + 
                    m*(T(0.00615234) + 
                    m*(T(0.00569661) + 
                    m*(T(0.00569661) + 
                    m*(T(0.00615234) + (T(0.00751953) + 
                    T(0.0139648)*m)*m))))) + 
                    y*(T(0.01155180) + 
                    m*(T(0.00616096) + 
                    m*(T(0.00497616) + 
                    m*(T(0.00452378) + 
                    m*(T(0.00439812) + m*(T(0.00452378) + 
                    m*(T(0.00497616) +   (T(0.00616096) + 
                    T(0.0115518)*m)*m)))))) + (T(0.00976161) + 
                    m*(T(0.00516791) + 
                    m*(T(0.00413433) + 
                    m*(T(0.00371030) + m*(T(0.00354165) + 
                    m*(T(0.00354165) + m*(T(0.00371030) + 
                    m*(T(0.00413433) +   (T(0.00516791) + 
                    T(0.00976161)*m)*m))))))))*y))))))))
end

"""
``\\text{asn}(u|\\,m)=\\text{sn}^{-1}(u,m)``.

Returns the inverse Jacobi Elliptic sn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic modulus
"""
function asn(s::A, m::B) where {A,B}
    T = promote_type(A,B)
    yA = T(0.04094) - T(0.00652)*m
    y = s * s
    y < yA && return s*serf(y, m)

    p = one(T)
    for _ in 1:10
        y = y * inv((1+√(1-y))*(1+√(1-m*y)))
        p += p
        y < yA && return p*√y*serf(y, m)
    end
    return T(NaN)
end

"""
``\\text{acn}(u|\\,m)=\\text{cn}^{-1}(u,m)``.

Returns the inverse Jacobi Elliptic cn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic modulus
"""
function acn(c::A, m::B) where {A,B}
    T = promote_type(A,B)
    mc = 1 - m
    x = c*c
    p = one(T)
    for _ in 1:10
        (2x > 1) && return p*asn(√(1-x), m)
        
        d = √(mc + m * x)
        x = (√x + d)/(1+d)
        p += p
    end
    return T(NaN)
end

@inline function rawF(sinφ::A, m::B) where {A,B}
    T = promote_type(A,B)

    yS = T(0.9000308778823196)
    m == 0 && return asin(sinφ)
    m == 1 && return atanh(sinφ)
    
    sinφ2 = sinφ*sinφ
    sinφ2 ≤ yS && return asn(sinφ, m)

    mc = 1 - m

    c = √(1-sinφ2)
    x = c * c
    d2 = mc + m*x
    #z = c/√(mc+m*c*c)
    x <  yS*d2 && return*(K(m) - asn(c/√(d2), m))

    v = mc*(1-x)
    v < x*d2 && return acn(c, m)
    return*(K(m) - acn(√(v/d2), m))

end

#----------------------------------------------------------------------------------------
# Elliptic F
#----------------------------------------------------------------------------------------
function _F(φ::A, m::B) where {A,B}
    T = promote_type(A,B)
    abs(φ) ≤ T(π/2) && return sign(φ)*rawF(sin(φ), m)
    j = floor(φ/T(π))
    newφ = φ - j*T(π)
    signφ = sign(newφ)
    if abs(newφ) > T(π/2)
        j += signφ*one(T)
        newφ = newφ - signφ*T(π)
    end
    signφ = sign(newφ)

    return 2j*K(m) + signφ*rawF(sin(abs(newφ)), m)
end

"""
``F(\\varphi |\\, m) = \\int_0^{\\varphi}\\frac{d\\theta}{\\sqrt{1-m\\sin(\\theta)^2}}.``

Returns the incomplete elliptic integral of the first kind.

# Arguments

- `φ` : Amplitude
- `m` : Elliptic modulus
"""
function F(φ::A, m::B) where {A,B}
    T = promote_type(A,B)
    if m > 1
        ## Abramowitz & Stegum*(17.4.15)
        m12 = sqrt(m)
        theta = asin(m12*sin(φ))
        signθ = sign(theta)
        absθ = abs(theta)
        return signθ/m12*_F(absθ, inv(m))
    elseif m < 0
        # Abramowitz & Stegum*(17.4.17)
        n = -m
        m12 = inv(sqrt(1+n))
        m1m = n/(1+n)
        newφ = T(π/2)-φ
        signφ = sign(newφ)
        absφ = abs(newφ)
        return (m12*K(m1m) - signφ*m12*_F(absφ, m1m)) 
    end
    absφ = abs(φ)
    signφ = sign(φ)
    return signφ*_F(absφ, m)
end

#----------------------------------------------------------------------------------------
# Elliptic D
#----------------------------------------------------------------------------------------
function D(m::T) where T
	return cel2(√(one(T) - m), zero(T), one(T))
end

function Dsi(s::A, m::B) where {A,B}
    T = promote_type(A,B)
	s2 = s^2
	return s * s2 *
		   (
			   T(1 / 3) +
			   s2 * (
				   T(1 / 10) + m / 10 +
				   s2 * (
					   T(3 / 56) + (T(1 / 28) + (3 * m) / 56) * m +
					   s2 * (
						   T(5 / 144) + m * (T(1 / 48) + (T(1 / 48) + (5 * m) / 144) * m) +
						   s2 * (
							   T(35 / 1408) + m * (T(5 / 352) + m * (T(9 / 704) + (T(5 / 352) + (35 * m) / 1408) * m)) +
							   s2 * (
								   T(63 / 3328) + m * (T(35 / 3328) + m * (T(15 / 1664) + m * (T(15 / 1664) +
																							   (T(35 / 3328) + (63 * m) / 3328) * m))) +
								   s2 * (
									   T(77 / 5120) + m * (T(21 / 2560) + m * (T(7 / 1024) + m * (T(5 / 768) +
																								  m * (T(7 / 1024) + (T(21 / 2560) + (77 * m) / 5120) * m)))) +
									   s2 * (
										   T(429 / 34816) + m * (T(231 / 34816) + m * (T(189 / 34816) +
																					   m * (T(175 / 34816) + m * (T(175 / 34816) + m * (T(189 / 34816) +
																																		(T(231 / 34816) + (429 * m) / 34816) * m))))) +
										   (T(6435 / 622592) + m * (T(429 / 77824) + m * (T(693 / 155648) +
																						  m * (T(315 / 77824) + m * (T(1225 / 311296) + m * (T(315 / 77824) +
																																			 m * (T(693 / 155648) + (T(429 / 77824) + (6435 * m) / 622592) *
																																									m))))))) * s2
									   )
								   )
							   )
						   )
					   )
				   )
			   )
		   )
end

function Ds(s::A, m::B) where {A,B}
    T = promote_type(A,B)
	_ybuf = @MArray [zero(T) for _ in 1:10]#

	y0 = s^2
	yi = y0

	yA(m::T) where T = T == Float32 ? (T(0.1888) - T(0.0378) * m) : (T(0.04094) - T(0.00652) * m)
	I = 0
	for _ in 1:10
		yi < yA(m) && break
		I += 1
		_ybuf[I] = yi

		ci = √(one(T) - yi)
		di = √(one(T) - m * yi)
		yi = yi / ((1 + ci) * (1 + di))
	end

	D0 = Dsi(sqrt(yi), m)
	I == zero(T) && return D0
	Di = D0

	for i in I:-1:1
		yim1 = _ybuf[i]
		sim1 = √yim1

		Di = 2Di + sim1 * yi

		yi = yim1
	end
	return Di
end

function rawD(φ::A, m::B) where {A,B}
    T = promote_type(A,B)
	ys = T == Float32 ? T(0.95) : T(0.9) # T(0.95) for single
	φs = T == Float32 ? T(1.345) : T(1.249) # T(1.345) for single
	s = sin(φ)

	if φ < φs
		return Ds(s, m)
	end

	mc = one(T) - m
	c = cos(φ)
	z = c / √(mc + m * c^2)
	if z^2 < ys
		return D(m) - Ds(z, m) - s * z
	end

	w = √(1 - z^2)
	if c < w
		return Ds(s, m)
	end

	return D(m) - Ds(sqrt(one(T) - w^2), m) - √((one(T) - c^2) * (one(T) - w^2))
end

# https://www.sciencedirect.com/science/article/pii/S0377042711001270?ref=pdf_download&fr=RR-2&rr=80442ff79ecc4cef
function D(φ::A, m::B) where {A,B}
    T = promote_type(A,B)

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

		return 2 * j * D(m) + signφ * D(abs(newφ), m)
	end

	φ < zero(T) && return -D(-φ, m)

	#Reduction of Parameter
	if m > one(T) && m * sin(φ)^2 ≤ one(T)
		mr = inv(m)
		φr = asin(√m * sin(φ))
		return mr * √mr * D(φr, mr)
    end

	#if m < zero(T)
	#	mN = -m / (1 - m)
	#	φN = asin(√(mc / (1 - m * sin(φ)^2)) * sin(φ))
	#	return √(1 - mN) * (D(φN, mN) - sin(φN) * cos(φN) / √(1 - mN * sin(φN)^2))
	#end


	return rawD(φ, m)
end

#----------------------------------------------------------------------------------------
# Elliptic E
#----------------------------------------------------------------------------------------

"""
``E(\\varphi|\\, m) = \\int_0^{\\varphi}d\\theta \\sqrt{1-m\\sin(\\theta)^2}.``

Returns the incomplete elliptic integral of the second kind.

# Arguments

- `φ` : Amplitude
- `m` : Elliptic modulus

"""
function E(φ::A, m::B) where {A,B}
    T = promote_type(A,B)
	return F(φ, m) - m * D(φ, m)
end

#----------------------------------------------------------------------------------------
# Elliptic Π
#----------------------------------------------------------------------------------------
Π(n, m) = Pi(n, m)
Π(n, φ, m) = Pi(n, φ, m)
#https://doi.org/T(10).1016/j.cam.2011.1107
"""
``\\Pi(n|\\,m)=\\int_{0}^{1 }{\\frac{1}{1-nt^{2}}}{\\frac{dt}{\\sqrt{\\left(1-mt^{2}\\right)\\left(1-t^{2}\\right)}}}.``

Returns the complete elliptic integral of the third kind.

# Arguments

- `n` : Characteristic
- `m` : Elliptic modulus
"""
function Pi(n::A, m::B) where {A,B}
    T = promote_type(A,B)
    n > one(T) && return K(m) - Pi(m/n, m)
    n == zero(T) && return K(m)
    m == zero(T) || m == one(T) && return T(Inf) #atanh(√(-1 + n)*tan(θ))/√(-1 + n)
    kc = √(one(T)-m)
    nc = one(T)-n
    return cel(kc, nc, one(T), one(T))
end

"""
``\\Pi (n;\\varphi \\,|\\,m)=\\int _{0}^{\\sin \\varphi }{\\frac {1}{1-nt^{2}}}{\\frac {dt}{\\sqrt {\\left(1-mt^{2}\\right)\\left(1-t^{2}\\right)}}}.``

Returns the incomplete elliptic integral of the third kind.

# Arguments

- `n` : Characteristic
- `φ` : Amplitude
- `m` : Elliptic modulus
"""
function Pi(n::A, φ::B, m::C) where {A,B,C}
    T = promote_type(A,B,C)
    if m < zero(T) # Imaginary modulus transformation https://dlmf.nist.gov/19.7#iii
        mc = one(T) - m
        imc = inv(mc) 
        mN = -m*imc
        φN = asin(sqrt(mc/(one(T)−m*sin(φ)^2)) * sin(φ) ) 

        nN = (n-m)*imc

        return sqrt(imc)/nN*(mN*F(φN, mN) + imc*n*Pi(nN, φN, mN))
    end # https://link.springer.com/book/10.1007/978-3-642-65138-0 117.01
    if n > one(T)
        nc = one(T) - n
        t1 = tan(φ)/sqrt(one(T) − m*sin(φ)^2) 
        h1 = nc*(n−m)/n 
        n1 = m/n
        return (FukushimaT(t1, h1) - n1*J(n1, φ, m))
    end
    return n*J(n, φ, m) + F(φ, m)
    
end

"""
``J (n;\\varphi \\,|\\,m)=\\frac{\\Pi(n;\\varphi|\\, m) - F(\\varphi|\\,m)}{n}.``

Returns the associate incomplete elliptic integral of the third kind.

# Arguments

- `n` : Characteristic
- `φ` : Amplitude
- `m` : Elliptic modulus
"""
function J(n::A, φ::B, m::C) where {A, B, C} #Appendix A
    T = promote_type(A,B,C)
    # Reduction of Amplitude
    φ == zero(T) && return zero(T)
    φ == T(π/2) && return J(n,m)

    if abs(φ) > T(π/2) && m < one(T)
        j = floor(φ/T(π))
        newφ = φ - j*T(π)
        signφ = sign(newφ)
        if abs(newφ) > T(π/2)
            j += signφ*one(T)
            newφ = newφ - signφ*T(π)
        end
        signφ = sign(newφ)

        return 2*j*J(n, m) + signφ*J(n, abs(newφ), m)
    end
    φ < zero(T) && return -J(n, -φ, m)

    # Reduction of parameter
    if zero(T) < φ < T(π/2) && m*sin(φ)^2 ≤ one(T)
        nc = one(T) - n
        iszero(m) && !iszero(n) && return (FukushimaT(tan(φ), nc) - φ) / n

        iszero(m) && iszero(n) && return φ/2 - sin(2*φ)/4

        isone(m) && !isone(n) && return (atanh(sin(φ)) - FukushimaT(sin(φ), -n))/nc

        isone(m) && isone(n) && return (sin(φ)/cos(φ)^2 - atanh(sin(φ)))/2
    end

    if one(T) < m < inv(sin(φ)^2)
        φR = asin(√m*sin(φ))
        nR = n/m
        mR = inv(m)
        #return NaN
        return mR*√mR*J(nR, φR, mR)
    elseif m < zero(T)
        mc = one(T) - m
        φN = asin(sqrt(mc/(one(T)−m*sin(φ)^2)) * sin(φ) ) 
        nN = (n−m)/mc 
        mN = -m/mc
        return mN*√mN*J(nN, φN, mN)
    end
    # Reduction of Characteristics
    if zero(T) < φ < one(T) && zero(T) < m < one(T) 
        if n > one(T)
            # t = inv(x^2)
            # h = y - x
            nc = one(T) - n
            t1 = tan(φ)/sqrt(one(T)−m*sin(φ)^2) 
            h1 = nc*(n−m)/n 
            n1 = m/n
            return (-F(φ, m) + FukushimaT(t1, h1) - n1*J(n1, φ, m))/n
        elseif n < zero(T)
            mc = one(T) - m
            nc = one(T) - n
            t2 = sin(φ)*cos(φ)/sqrt(one(T)−m*sin(φ)^2) 
            h2 = -n*(m-n)/nc
            n2 = (m−n)/nc 
            #return NaN
            return(F(φ, m) - FukushimaT(t2, h2) - (mc/nc)*J(n2, φ, m))/nc
        end
    end
    
    zero(T) < φ < T(π/2) && zero(T) < m < one(T) && zero(T) < n < one(T) && return rawJ(n, φ, m)
    return T(NaN)
end

function rawJ(n::A, φ::B ,m::C) where {A,B,C}
    T = promote_type(A,B,C)
    mc = one(T) - m
    nc = one(T) - n
    ys = T == Float32 ? T(0.95) : T(0.9) # T(0.95) for single
    φs = T == Float32 ? T(1.345) : T(1.249) # T(1.345) for single
    s = sin(φ)
    c = cos(φ)
    z = c/√(mc + m*c^2)
    w = √(mc/(one(T)-m*s^2))*s
    h = n*nc*(n-m)
    ts = s*z/nc
    tc = √((one(T)-c^2)*(one(T)-w^2))/nc
    if φ < φs
        return Js(n, s, m)
    elseif z^2 < ys
        return J(n, m) - Js(n, z, m) - FukushimaT(ts, h)
    elseif c < w
        return Js(n, s, m)
    else
        return J(n, m) - Js(n, sqrt(one(T)-w^2), m) - FukushimaT(tc, h)
    end
end

"""
``J(n;\\varphi \\,|\\,m)=\\frac{\\Pi(n;\\pi/2|\\, m) - K(m)}{n}.``

Returns the associate complete elliptic integral of the third kind.

# Arguments

- `n` : Characteristic
- `m` : Elliptic modulus
"""
function J(n::A, m::B) where {A,B}
    T = promote_type(A,B)
    n > one(T) && return m/n*J(m/n, m)
    kc = √(one(T)-m)
    nc = one(T)-n
    return cel(kc, nc, zero(T), one(T))
end

#function Jc(n::T, c::T, m::T) where{T}
#    #_ybuf = @SVector [zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T)]#
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
#        #@set! _ybuf[I] = one(T) - xi
#        @set! _xbuf[I] = xi
#        ci = √xi
#        di = √(mc+m*xi)
#        xi = (ci +di)/(one(T)+di)
#    end
#    J0 = JsI(n, one(T)-xi, m)
#    I == 0 && return J0
#    Ji = J0
#    yi = one(T) - xi 
#    ti = 0
#    for i in I:-1:1
#        #yim1 = _ybuf[i]
#        #sim1 = √yim1
#        #cim1 = √(one(T)-yim1)
#        #dim1 = √(one(T)-m*yim1)
#        
#        #ti = sim1*yi / (one(T)-n*(yim1 - cim1*dim1*yi))
#        #Ji = 2Ji + FukushimaT(ti, h)
#
#        #yi = yim1
#
#        xim1 = _xbuf[i]
#        yim1 = one(T) - xim1
#        sim1 = √yim1
#        cim12 = (xim1)
#        dim12 = (one(T)-m*yim1)
#        
#        ti = sim1*yi / (one(T)-n*(yim1 - sqrt(cim1*dim1)*yi))
#        Ji = 2Ji + FukushimaT(ti, h)
#
#        xi = xim1
#
#    end
#    return Ji
#
#end

function Js(n::A, s::B, m::C) where {A,B,C}
    T = promote_type(A,B,C)
    _ybuf = @SVector [zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T)]#

    nc = one(T) - n
    h = n*nc*(n-m)

    y0 = s^2

    yi = y0 

    yB = T == Float32 ? T(0.05905) : T(0.01622)
    I = 0
    for _ in 1:10
        yi < yB && break
        I += 1
        @set! _ybuf[I] = yi

        ci = √(one(T)-yi)
        di = √(one(T)-m*yi)
        yi = yi/((1 + ci)*(1 + di))
    end

    J0 = JsI(n, yi, m)
    I == zero(T) && return J0
    Ji = J0
    ti = zero(T)
    
    for i in I:-1:1
        yim1 = _ybuf[i]
        sim1 = √yim1
        cim1 = √(one(T)-yim1)
        dim1 = √(one(T)-m*yim1)

        ti = sim1*yi / (one(T)-n*(yim1 - cim1*dim1*yi))
        Ji = Ji+Ji + FukushimaT(ti, h)

        yi = yim1
    end
    return Ji

end

function JsI(n::A, y::B, m::C) where {A,B,C}
    T = promote_type(A,B,C)
    #Series[(EllipticPi[n, ArcSin[s], m] − EllipticF[ArcSin[s], m])/n, {s, zero(T), 19}]
    return sqrt(y)*y*(
        T(0.3333333333333333) + y*(
            T(0.1) + T(0.1)*m + T(0.2)*n + y*(
                T(0.05357142857142857) + n*(T(0.07142857142857142) + T(0.14285714285714285)*n) + m*(T(0.03571428571428571) + T(0.05357142857142857)*m + T(0.07142857142857142)*n) + y*(
                    T(0.034722222222222224) + m*(T(0.020833333333333332) + m*(T(0.020833333333333332) + T(0.034722222222222224)*m + T(0.041666666666666664)*n) + n*(T(0.027777777777777776) + T(0.05555555555555555)*n)) + n*(T(0.041666666666666664) + n*(T(0.05555555555555555) + T(0.1111111111111111)*n)) + y*(
                        T(0.024857954545454544) + n*(T(0.028409090909090908) + n*(T(0.03409090909090909) + n*(T(0.045454545454545456) + T(0.09090909090909091)*n))) + m*(T(0.014204545454545454) + m*(T(0.01278409090909091) + m*(T(0.014204545454545454) + T(0.024857954545454544)*m + T(0.028409090909090908)*n) + n*(T(0.017045454545454544) + T(0.03409090909090909)*n)) + n*(T(0.017045454545454544) + n*(T(0.022727272727272728) + T(0.045454545454545456)*n))) + y*(
                            T(0.01893028846153846) + m*(T(0.010516826923076924) + n*(T(0.01201923076923077) + n*(T(0.014423076923076924) + n*(T(0.019230769230769232) + T(0.038461538461538464)*n))) + m*(T(0.009014423076923076) + m*(T(0.009014423076923076) + m*(T(0.010516826923076924) + T(0.01893028846153846)*m + T(0.021033653846153848)*n) + n*(T(0.01201923076923077) + T(0.02403846153846154)*n)) + n*(T(0.010817307692307692) + n*(T(0.014423076923076924) + T(0.028846153846153848)*n)))) + n*(T(0.021033653846153848) + n*(T(0.02403846153846154) + n*(T(0.028846153846153848) + n*(T(0.038461538461538464) + T(0.07692307692307693)*n)))) + y*(
                                T(0.0150390625) + m*(T(0.008203125) + m*(T(0.0068359375) + m*(T(0.006510416666666667) + n*(T(0.0078125) + n*(T(0.010416666666666666) + T(0.020833333333333332)*n)) + m*(T(0.0068359375) + m*(T(0.008203125) + T(0.0150390625)*m + T(0.01640625)*n) + n*(T(0.009114583333333334) + T(0.01818181818181818)*n))) + n*(T(0.0078125) + n*(T(0.009375) + n*(T(0.0125) + T(0.025)*n)))) + n*(T(0.009114583333333334) + n*(T(0.010416666666666666) + n*(T(0.0125) + n*(T(0.016666666666666666) + T(0.03333333333333333)*n))))) + n*(T(0.01640625) + n*(T(0.018229166666666668) + n*(T(0.020833333333333332) + n*(T(0.025) + n*(T(0.03333333333333333) + T(0.06666666666666667)*n))))) + y*(
                                    T(0.012321920955882353) + m*(T(0.006634880514705882) + m*(T(0.005428538602941176) + m*(T(0.0050264246323529415) + n*(T(0.0057444852941176475) + n*(T(0.006893382352941176) + n*(T(0.009191176470588236) + T(0.01838235294117647)*n))) + m*(T(0.0050264246323529415) + m*(T(0.005428538602941176) + m*(T(0.006634880514705882) + T(0.012321920955882353)*m + T(0.013269761029411764)*n) + n*(T(0.007238051470588236) + T(0.014476102941176471)*n)) + n*(T(0.00603170955882353) + n*(T(0.008042279411764705) + T(0.01608455882352941)*n)))) + n*(T(0.00603170955882353) + n*(T(0.006893382352941176) + n*(T(0.008272058823529412) + n*(T(0.011029411764705883) + T(0.022058823529411766)*n))))) + n*(T(0.007238051470588236) + n*(T(0.008042279411764705) + n*(T(0.009191176470588236) + n*(T(0.011029411764705883) + n*(T(0.014705882352941176) + T(0.029411764705882353)*n)))))) + n*(T(0.013269761029411764) + n*(T(0.014476102941176471) + n*(T(0.01608455882352941) + n*(T(0.01838235294117647) + n*(T(0.022058823529411766) + n*(T(0.029411764705882353) + T(0.058823529411764705)*n)))))) + y*(
                                        T(0.01033582185444079) + n*(T(0.011024876644736841) + n*(T(0.011872944078947368) + n*(T(0.012952302631578948) + n*(T(0.014391447368421052) + n*(T(0.01644736842105263) + n*(T(0.019736842105263157) + n*(T(0.02631578947368421) + T(0.05263157894736842)*n))))))) + m*(T(0.005512438322368421) + n*(T(0.005936472039473684) + n*(T(0.006476151315789474) + n*(T(0.007195723684210526) + n*(T(0.008223684210526315) + n*(T(0.009868421052631578) + n*(T(0.013157894736842105) + T(0.02631578947368421)*n)))))) + m*(T(0.004452354029605263) + m*(T(0.004047594572368421) + n*(T(0.004497327302631579) + n*(T(0.005139802631578948) + n*(T(0.006167763157894737) + n*(T(0.008223684210526315) + T(0.01644736842105263)*n)))) + m*(T(0.0039351613898026315) + n*(T(0.004497327302631579) + n*(T(0.005396792763157895) + n*(T(0.007195723684210526) + T(0.014391447368421052)*n))) + m*(T(0.004047594572368421) + m*(T(0.004452354029605263) + n*(T(0.005936472039473684) + T(0.011872944078947368)*n) + m*(T(0.005512438322368421) + T(0.01033582185444079)*m + T(0.011024876644736841)*n)) + n*(T(0.004857113486842105) + n*(T(0.006476151315789474) + T(0.012952302631578948)*n))))) + n*(T(0.004857113486842105) + n*(T(0.005396792763157895) + n*(T(0.006167763157894737) + n*(T(0.007401315789473684) + n*(T(0.009868421052631578) + T(0.019736842105263157)*n)))))))
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
end

function custom_atanh(a::T) where T

	arg1 = abs(one(T) + a)
	arg2 = abs(one(T) - a)

	ans = (log(arg1 / arg2)) / 2
	return ans
end

function FukushimaT(t::A, h::B) where {A,B}
    T = promote_type(A,B)
	if h > zero(T)
		return atan(t * √h) / √(h)
	elseif h == zero(T)
		return t
	else
        arg = t * √(-h)
        ans = abs(arg) < one(T) ? atanh(arg) : custom_atanh(arg)
        return ans / √(-h)
	end
end

#https://link-springer-com.ezp-prod1.hul.harvard.edu/article/T(10).1007/BF02165405
function cel(kc::A, p::B, a::C, b::D) where {A,B,C,D}
    T = promote_type(A,B,C,D)
    #ca = T(1e-6)
    ca = eps(T)
    kc = abs(kc)
    e = kc
    m = one(T)

    f, g, q = T(0), T(0), T(0)
    if p > T(0)
        p = √p
        b = b/p
    else
        f = kc^2
        q = one(T) - f
        g = one(T) - p
        f = f-p
        q = (b-a*p)*q
        p = √(f/g)
        a = (a-b)/g
        b = -q*(g^2*p) + a*p
    end
    while true
        f = a
        a = b/p + a
        g = e/p
        b = f*g + b
        b = b + b
        p = g+p
        g = m
        m = kc+m
        if abs(g-kc) < g*ca
            break
        end
        kc = √e
        kc = kc+kc
        e = kc*m
    end
    return T(π/2)*(a*m+b)/(m*(m+p))
end

#https://doi-org.ezp-prod1.hul.harvard.edu/T(10).031007/s10569-008-9177-y
#https://link.springer.com/article/T(10).1007/s10569-008-9177-y

function cel2(kc::A, a::B, b::C) where {A,B,C}
    T = promote_type(A,B,C)
	return cel(kc, one(T), a, b)
end

function _Kscreen(m::T) where T
    return T(π/2)*(one(T) + m*(T(0.25) + m*(T(0.36) + m*(T(0.09765625) + m*T(0.07476806640625)))))
end

function _ΔXNloop(u::A, m::B, n::C) where {A,B,C}
    T = promote_type(A,B,C)
    up = u * inv(T(2)^n)
    up2 = up*up
    sn = up*(up2*(up2*((m*((-(m/T(5040))-T(3/112))*m-T(3/112))-T(1/5040))*up2+(m/T(120)+T(7/60))*m+T(1/120))-m/T(6)-T(1/6))+one(T))
    cn = one(T)+up2*(-(T(1/2))+up2*(T(1/24)+m/T(6)+up2*(-(T(1/720))+(-(T(11/180))-m/T(45))*m+(-(T(1/40320))+m*(-(T(17/1680))+(-(T(19/840))-m/T(630))*m))*up2)))
    dn = m*(m*(up2^2*(inv(T(24))-(T(11)*up2)/T(180))-(m*up2^3)/T(720)) + (up2*(T(1/6)-up2/T(45))-T(0.5))*up2) + one(T) 
    Δsn = up-sn
    Δcn = one(T)-cn
    Δdn = one(T)-dn

    i = one(T)
    while i <= n
        sn2 = sn*sn
        sn4 = sn2*sn2
        
        den = inv(one(T)-m*sn4)
        Δsn = T(2.0)*(Δsn*cn*dn + up*(Δcn*(one(T) - Δdn) + Δdn  - m*sn4))*den
        Δcn = ((one(T)+cn)*Δcn + sn2-2*m*sn4)*den
        Δdn = ((one(T)+dn)*Δdn + m*(sn2-2*sn4))*den

        up += up
        sn = up - Δsn
        cn = one(T) - Δcn
        dn = one(T) - Δdn

        i += one(T)
    end
    return sn, cn, dn
end

function fold_0_25(u1::T, m::T, kp::T) where T 
    u1 == zero(T) && return zero(T), one(T), one(T)

    sn, cn, dn = _ΔXNloop(u1, m, u1 > zero(T) ? max(6+(floor(log2(u1))), one(T)) : zero(T))
    den = inv(one(T)+kp-m*sn*sn)
    return den*√(one(T)+kp)*(cn*dn-kp*sn), den*√(kp*(one(T)+kp))*(cn+sn*dn), den*√kp*((one(T)+kp)*dn+m*sn*cn)
end

function fold_0_50(u1::T, m::T, Kscreen::T, Kactual::T, kp::T) where T
    u1 == zero(T) && return zero(T), one(T), one(T)

    if u1 > T(0.25)*Kscreen 
        sn, cn, dn = fold_0_25(Kactual/2 - u1, m, kp)
    else
        sn, cn, dn = _ΔXNloop(u1, m, u1 > zero(T) ?  max(6+(floor(log2(u1))), one(T)) : zero(T))
    end
    return cn/dn, kp*sn/dn, kp/dn
end

function fold_1_00(u1::T, m::T, Kscreen::T, Kactual::T, kp::T) where T
    u1 == zero(T) && return zero(T), one(T), one(T)

    if u1 > Kscreen/2
        sn, cn, dn = fold_0_50(Kactual - u1, m, Kscreen, Kactual, kp)
    elseif u1 > Kscreen/4
        sn, cn, dn = fold_0_25(Kactual/2 - u1, m, kp)
    else
        sn, cn, dn = _ΔXNloop(u1, m, u1 > zero(T) ?  max(6+(floor(log2(u1))), one(T)) : zero(T))
    end
    return cn/dn, -kp*sn/dn, kp/dn
end

function _SN(u::T, m::T) where T
    tempK = K(m)
    if u > (4tempK)
        return rawSN(u % (4tempK), m, tempK, tempK, √(one(T)-m))
    end
    return rawSN(u, m, tempK, tempK, √(one(T)-m))
end

function rawSN(u::T, m::T, Kscreen::T, Kactual::T, kp::T) where T 
    check = u ≥ 2*Kscreen 
    sign = check ? -one(T) : one(T)
    u = check ? u - 2*Kactual : u
    u > Kscreen && return sign*fold_1_00(u - Kactual, m, Kscreen, Kactual, kp)[1]
    u == Kscreen && return one(T)
    u > Kscreen/2 && return sign*fold_0_50(Kactual - u, m, Kscreen, Kactual, kp)[1]
    u == Kscreen/2 && return inv(√(one(T)+kp))
    u ≥ Kscreen/4 && return sign*fold_0_25(Kactual/2 - u, m, kp)[1]
    return sign*_ΔXNloop(u, m, u > zero(T) ? max(6+(floor(log2(u))), one(T)) : zero(T))[1]
end

function _CN(u::T, m::T) where T
    tempK = K(m)
    if u > (4tempK)
        return rawCN(u % (4tempK), m, tempK, tempK, √(one(T)-m))
    end
    return rawCN(u, m, tempK, tempK, √(one(T)-m))
end
function rawCN(u::T, m::T, Kscreen::T, Kactual::T, kp::T) where T 
    check = u ≥ 2*Kscreen 
    sign = check ? -one(T) : one(T)
    u = check ? u - 2*Kactual : u
    u > Kscreen && return sign*fold_1_00(u - Kactual, m, Kscreen, Kactual, kp)[2]
    u == Kscreen && return zero(T)
    u > Kscreen/2 && return sign*fold_0_50(Kactual - u, m, Kscreen, Kactual, kp)[2]
    u == Kscreen/2 && return √(kp/(one(T)+kp))
    u ≥ Kscreen/4 && return sign*fold_0_25(Kactual/2 - u, m, kp)[2]
    return sign*_ΔXNloop(u, m, u > zero(T) ? max(6+(floor(log2(u))), one(T)) : zero(T))[2]
end

function _DN(u::T, m::T) where T
    tempK = K(m)
    if u > 4tempK
        return rawDN(u % (4tempK), m, tempK, tempK, √(one(T)-m))
    end
    return rawDN(u, m, tempK, tempK, √(one(T)-m))
end
function rawDN(u::T, m::T, Kscreen::T, Kactual::T, kp::T) where T 
    check = u ≥ 2*Kscreen 
    u = check ? u - 2*Kactual : u
    u > Kscreen && return fold_1_00(u - Kactual, m, Kscreen, Kactual, kp)[3]
    u == Kscreen && return kp
    u > Kscreen/2 && return fold_0_50(Kactual - u, m, Kscreen, Kactual, kp)[3]
    u == Kscreen/2 && return √(kp)
    u ≥ Kscreen/4 && return fold_0_25(Kactual/2 - u, m, kp)[3]
    return _ΔXNloop(u, m, u > zero(T) ? max(6+(floor(log2(u))), one(T)) : zero(T))[3]
end

_sc_helper(jacobituple) = jacobituple[1]/jacobituple[2]
function rawSC(u::T, m::T, Kscreen::T, Kactual::T, kp::T) where T 
    u = u ≥ 4*Kscreen ?  u % 4*Kactual : u
    u = u ≥ 2*Kscreen ? u - 2*Kactual : u
    u > Kscreen && return _sc_helper(fold_1_00(u - Kactual, m, Kscreen, Kactual, kp))
    u == Kscreen && return T(Inf)
    u > Kscreen/2 && return _sc_helper(fold_0_50(Kactual - u, m, Kscreen, Kactual, kp))
    u == Kscreen/2 && return inv(√(kp))
    u ≥ Kscreen/4 && return _sc_helper(fold_0_25(Kactual/2 - u, m, kp))
    return _sc_helper(_ΔXNloop(u, m, u > zero(T) ? max(6+(floor(log2(u))), one(T)) : zero(T)))
end

function _SC(u::T, m::T) where T 
    return rawSC(u, m, K(m), K(m), √(one(T)-m))
end

_sd_helper(jacobituple) = jacobituple[1]/jacobituple[3]
function rawSD(u::T, m::T, Kscreen::T, Kactual::T, kp::T) where T 
    u = u ≥ 4*Kscreen ?  u % 4*Kactual : u
    u = u ≥ 2*Kscreen ? u - 2*Kactual : u
    u > Kscreen && return _sd_helper(fold_1_00(u - Kactual, m, Kscreen, Kactual, kp))
    u == Kscreen && return inv(kp)
    u > Kscreen/2 && return _sd_helper(fold_0_50(Kactual - u, m, Kscreen, Kactual, kp))
    u == Kscreen/2 && return inv(√((1+kp)*kp))
    u ≥ Kscreen/4 && return _sd_helper(fold_0_25(Kactual/2 - u, m, kp))
    return _sd_helper(_ΔXNloop(u, m, u > zero(T) ? max(6+(floor(log2(u))), one(T)) : zero(T)))
end

function _SD(u::T, m::T) where T
    return rawSD(u, m, K(m), K(m), √(one(T)-m))
end

"""
``\\text{sn}(u|\\,m)=\\sin\\,\\text{am}(u,m)``, where ``\\text{am}(u|\\,m)=F^{-1}(u|\\,m)`` 
is the Jacobi amplitude.

Returns the Jacobi Elliptic sn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic modulus
"""
function sn(u::T, m::T) where T  
    signu = sign(u)
    u = abs(u)
    m < one(T) && return signu*_SN(u, m)
    sqrtm = √m
    return signu*_SN(u*sqrtm, inv(m))/sqrtm
end

"""
``\\text{cn}(u|\\,m)=\\cos\\,\\text{am}(u,m)``, where ``\\text{am}(u|\\,m)=F^{-1}(u|\\,m)`` 
is the Jacobi amplitude.

Returns the Jacobi Elliptic cn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic modulus
"""
function cn(u::T, m::T) where T  
    u = abs(u)
    m < one(T) && return _CN(u, m)
    sqrtm = √m
    return _DN(u*sqrtm, inv(m))
end

"""
``\\text{dn}(u|\\,m)=\\frac{d}{d u}\\,\\text{am}(u,m)``, where ``\\text{am}(u|\\,m)=F^{-1}(u|\\,m)`` 
is the Jacobi amplitude.

Returns the Jacobi Elliptic dn.

# Arguments

- `u` : Amplitude
- `m` : Elliptic modulus
"""
function dn(u::T, m::T) where T  
    u = abs(u)
    m < one(T) && return _DN(u, m)
    sqrtm = √m
    return _CN(u*sqrtm, inv(m))
end

function sc(u::T, m::T) where T
    signu = sign(u)
    u = abs(u)

    m < one(T) && return signu*_SC(u, m)
    sqrtm = √m
    return signu*_SD(u*sqrtm, inv(m))/sqrtm
end

function sd(u::T, m::T) where T
    signu = sign(u)
    u = abs(u)

    m < one(T) && return signu*_SD(u, m)
    sqrtm = √m
    return signu*_SC(u*sqrtm, inv(m))*inv(sqrtm)
end

xn = ((:s,:(sn(u,m))), (:c,:(cn(u,m))), (:d,:(dn(u,m))), (:n,:(1.)))
for (p,num) in xn, (q,den) in xn
    f = Symbol(p, q)
    #@eval begin
    #    """
    #        $($f)(u::Real, m::Real)

    #    Compute the Jacobi elliptic function $($f)(u | m)
    #    """
    #    ($f)(u::Real, m::Real) = ($f)(Float64(u), Float64(m))
    #end

    if (p == q)
        @eval ($f)(::A, ::B) where {A,B} = begin one(promote_type(A,B)) end
    elseif (q != :n)
        @eval ($f)(u::A, m::B) where {A,B} = begin ($num)/($den) end
    end
end

function am(u::T, m::T) where T
    asin(sn(u,m))
end

end
