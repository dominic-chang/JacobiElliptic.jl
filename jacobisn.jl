include("ellipe.jl")
#https://doi-org.ezp-prod1.hul.harvard.edu/10.031007/s10569-008-9177-y
HALF_PI = 1.5707963267948966

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

function _SN(u, m, Kscreen, Kactual, kp)
    u == 0.0 && return zero(u)
    m == 0.0 && return sin(u)
    m == 1.0 && return tanh(u)

    
    u > 4Kscreen &&  return _SN(u % 4Kactual, m, Kscreen, Kactual, kp)
    u > 2Kscreen && return -_SN(u - 2Kactual, m, Kscreen, Kactual, kp)
    u > Kscreen && return _CN(u - Kactual, m, Kscreen, Kactual, kp)/_DN(u - Kactual, m, Kscreen, Kactual, kp)
    u >= 0.5Kscreen && return _CN(Kactual - u, m, Kscreen, Kactual, kp)/_DN(Kactual - u, m, Kscreen, Kactual, kp)
    u > 0.25Kscreen && return (up = Kactual/2 - u;√(1+kp)*(_CN(up,m,Kscreen,Kactual, kp)*_DN(up,m,Kscreen,Kactual, kp)-kp*_SN(up,m,Kscreen,Kactual, kp))/(1+kp-m*_SN(up,m,Kscreen,Kactual, kp)^2))
    u > 0.031 && return _ΔSN(u, m, Kscreen, Kactual, kp) - u 

    return u*(u^2*(u^2*((m*((-(m/5040)-3/112)*m-3/112)-1/5040)*u^2+(m/120+7/60)*m+1/120)-m/6-1/6)+1)

end

function _CN(u, m, Kscreen, Kactual, kp)
    u == zero(u) && return one(u)
    m == 0.0 && return cos(u)
    m == 1.0 && return sech(u)

    
    u > 4Kscreen &&  return _CN(u % 4Kactual, m, Kscreen, Kactual, kp)
    u > 2Kscreen &&  return -_CN(u - 2Kactual, m, Kscreen, Kactual, kp)
    u > Kscreen && return -kp*_SN(u - Kactual, m, Kscreen, Kactual, kp)/_DN(u - Kactual, m, Kscreen, Kactual, kp)
    u >= 0.5Kscreen && return kp*_SN(Kactual - u, m, Kscreen, Kactual, kp)/_DN(Kactual - u, m, Kscreen, Kactual, kp)
    u > 0.25Kscreen && (up = Kactual/2 - u; return √(kp*(1+kp))*(_CN(up,m,Kscreen,Kactual, kp)+_SN(up,m,Kscreen,Kactual, kp)*_DN(up,m,Kscreen,Kactual, kp))/(1+kp-m*_SN(up,m,Kscreen,Kactual, kp)^2))
    u > 0.031 && return _ΔCN(u, m, Kscreen, Kactual, kp) - 1

    return u*(u^2*(u^2*((m*((-(m/5040)-3/112)*m-3/112)-1/5040)*u^2+(m/120+7/60)*m+1/120)-m/6-1/6)+1)

end


function _DN(u, m, Kscreen, Kactual, kp)
    u == zero(u) && return one(u)
    m == 0.0 && return one(u)
    m == 1.0 && return sech(u)
    
    u > 4.0Kscreen &&  return _DN(u % 4.0Kactual, m, Kscreen, Kactual, kp)
    u >= 2.0Kscreen && return _DN(u - 2.0Kactual, m, Kscreen, Kactual, kp)
    u >= Kscreen && return kp/_DN(u - Kactual, m, Kscreen, Kactual, kp)
    u >= 0.5Kscreen && return kp/_DN(Kactual - u, m, Kscreen, Kactual, kp)
    u > 0.25Kscreen && return (up = Kactual/2 - u; √kp*((1+kp)*(_DN(up,m,Kscreen,Kactual, kp)+m*_SN(up,m,Kscreen,Kactual, kp)*_CN(up,m,Kscreen,Kactual, kp)))/(1+kp-m*_SN(up,m,Kscreen,Kactual, kp)^2))
    u > 0.031 && return _ΔDN(u, m, Kscreen, Kactual, kp) - 1

    return m*(m*(u^4*(1/24-(11*u^2)/180)-(m*u^6)/720)+(u^2*(1/6-u^2/45)-1/2)*u^2)+1

end

SN(u, m) = _SN(u, m, _Kscreen(m), _K(m), √(1-m))
CN(u, m) = _CN(u, m, _Kscreen(m), _K(m), √(1-m))
DN(u, m) = _DN(u, m, _Kscreen(m), _K(m), √(1-m))