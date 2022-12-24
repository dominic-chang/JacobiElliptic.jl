include("ellipe.jl")
#https://doi-org.ezp-prod1.hul.harvard.edu/10.031007/s10569-008-9177-y
const HALF_PI = 1.5707963267948966

_Kscreen(m) = HALF_PI*(1.0 + m*(0.25 + m*(0.36 + m*(0.09765625 + m*0.07476806640625))))

function _ΔSN(up, m, Kscreen, Kactual, kp)
    if up > 0.031
        u = up/2
        return 2/(1-m*_SN(u,m,Kscreen, Kactual, kp)^4)*(_CN(u,m,Kscreen,Kactual, kp)*_DN(u,m,Kscreen,Kactual, kp)*_ΔSN(u,m,Kscreen,Kactual, kp)+u*(_ΔCN(u,m,Kscreen,Kactual, kp)+_ΔDN(u,m,Kscreen,Kactual, kp)-_ΔCN(u,m,Kscreen,Kactual, kp)*_ΔDN(u,m,Kscreen,Kactual, kp)-m*_SN(u,m,Kscreen,Kactual, kp)^4))
    else
        return up - _SN(up, m, Kscreen, Kactual, kp)
    end
end
function _ΔCN(up, m, Kscreen, Kactual, kp)
    if up > 0.031
        u = up/2
        return ((1+_CN(u,m,Kscreen,Kactual, kp))*_ΔCN(u,m,Kscreen,Kactual, kp)+(1-2*m*_SN(u,m,Kscreen,Kactual, kp)^2)*_SN(u,m,Kscreen,Kactual, kp)^2)/(1-m*_SN(u,m,Kscreen,Kactual, kp)^4)
    else 
        return 1 - _CN(up, m, Kscreen, Kactual, kp)
    end
end
function _ΔDN(up, m, Kscreen, Kactual, kp) 
    if up > 0.031
        u = up/2 
        return ((1+_DN(u,m,Kscreen,Kactual, kp))*_ΔDN(u,m,Kscreen,Kactual, kp)+m*(1-2*_SN(u,m,Kscreen,Kactual, kp)^2)*_SN(u,m,Kscreen,Kactual, kp)^2)/(1-m*_SN(u,m,Kscreen,Kactual, kp)^4)
    else 
        return 1 - _CN(up, m, Kscreen, Kactual, kp)
    end
end

function _XNloop(u, m ,n)

    up = u / (2.0^n)
    sn = up*(up^2*(up^2*((m*((-(m/5040)-3/112)*m-3/112)-1/5040)*up^2+(m/120+7/60)*m+1/120)-m/6-1/6)+1)
    cn = 1 + up^2*(-(1/2)+up^2*(1/24+m/6+up^2*(-(1/720)+(-(11/180)-m/45)*m+(-(1/40320)+m*(-(17/1680)+(-(19/840)-m/630)*m))*up^2)))
    dn = m*(m*(up^4*(1/24-(11*up^2)/180)-(m*up^6)/720) + (up^2*(1/6-up^2/45)-1/2)*up^2) + 1 

    for _ in 1:n
        sn2 = sn^2
        cn2 = cn^2
        dn2 = dn^2

        den = 1/(1.0-0.5*sn2^2)
        sn = 2.0*(sn*cn*dn)*den
        cn = (cn2-sn2*dn2)*den
        dn = (dn2-m*sn2*cn2)*den
    end
    return sn, cn, dn
end
function fold_0_25(u1, m, kp) 
    u1 == 0 && return 0, 1, 1
    sn, cn, dn = _XNloop(u1, m, u1 > 0 ? max(6+Int(floor(log2(u1))), 1) : 0)
    den = 1/(1+kp-m*sn^2)
    return den*√(1+kp)*(cn*dn-kp*sn), den*√(kp*(1+kp))*(cn+sn*dn), den*√kp*((1+kp)*dn+m*sn*cn)
end

function fold_0_50(u1, m, Kscreen, Kactual, kp)
    u1 == 0 && return 0, 1, 1

    if u1 > 0.25Kscreen 
        sn, cn, dn = fold_0_25(Kactual/2 - u1, m, kp)
    else
        sn, cn, dn =  _XNloop(u1, m, u1 > 0 ?  max(6+Int(floor(log2(u1))), 1) : 0)
    end
    return cn/dn, kp*sn/dn, kp/dn
end

function fold_1_00(u1, m, Kscreen, Kactual, kp)
    u1 == 0 && return 0, 1, 1

    if u1 > 0.5Kscreen
        sn, cn, dn = fold_0_50(Kactual - u1, m, Kscreen, Kactual, kp)
    elseif u1 > 0.25Kscreen 
        sn, cn, dn = fold_0_25(Kactual/2 - u1, m, kp)
    else
        sn, cn, dn = _XNloop(u1, m, u1 > 0 ?  max(6+Int(floor(log2(u1))), 1) : 0)
    end
    return cn/dn, -kp*sn/dn, kp/dn
end



funcs = ((:SN, :_SN), (:CN, :_CN), (:DN, :_DN))
for (enum, funcpair) in enumerate(funcs)
    (func, helper) = funcpair
    @eval begin
        function $(helper)(u, m, Kscreen, Kactual, kp) 
            u = u > 4Kscreen ?  u % 4Kactual : u
            u = u > 2Kscreen ? u - 2Kactual : u
            u > Kscreen && return fold_1_00(u - Kactual, m, Kscreen, Kactual, kp)[$enum]
            u > 0.5Kscreen && return fold_0_50(Kactual - u, m, Kscreen, Kactual, kp)[$enum]
            u > 0.25Kscreen && return fold_0_25(Kactual/2 - u, m, kp)[$enum]
            return _XNloop(u, m, u > 0 ? max(6+Int(floor(log2(u))), 1) : 0)[$enum]
        end

        $(func)(u, m) = $(helper)(u, m, _Kscreen(m), _K(m), √(1-m))
    end
end

