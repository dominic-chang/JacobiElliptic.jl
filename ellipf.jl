include("ellipe.jl")

#Fast computation of incomplete elliptic integral of firstkind by half argument transformation
function serf(y, m)
    return 1 + y*(1/6 + m/6 + 
        y*(3/40 + (1/20 + (3*m)/40)*m + 
        y*(5/112 + m*(3/112 + (3/112 + (5*m)/112)*m) + 
            y*(35/1152 + 
                m*(5/288 + m*(1/64 + (5/288 + (35*m)/1152)*m)) + 
                y*(63/2816 + 
                    m*(35/2816 + 
                    m*(15/1408 + 
                        m*(15/1408 + (35/2816 + (63*m)/2816)*m))) + 
                    y*(231/13312 + 
                    m*(63/6656 + 
                        m*(105/13312 + 
                        m*(25/3328 + 
                        m*(105/13312 + (63/6656 + (231*m)/13312)*m)))) + 
                    y*(143/10240 + 
                        m*(77/10240 + 
                        m*(63/10240 + 
                        m*(35/6144 + 
                        m*(35/6144 + m*(63/
                        10240 + (77/10240 + (143*m)/10240)*m))))) + 
                        y*(6435/557056 + 
                        m*(429/69632 + 
                        m*(693/139264 + 
                        m*(315/69632 + 
                        m*(1225/278528 + m*(315/69632 + 
                        m*(693/139264 + (429/69632 + (6435*m)/
                        557056)*m)))))) + (12155/1245184 + 
                        m*(6435/1245184 + 
                        m*(1287/311296 + 
                        m*(1155/311296 + m*(2205/622592 + 
                        m*(2205/622592 + m*(1155/311296 + 
                        m*(1287/311296 + (6435/1245184 + (12155*m)/
                        1245184)*m))))))))*y))))))))
 
end

function asn(s, m)
    yA = 0.04094 - 0.00652*m
    y = s * s
    if y < yA
        return s*serf(y, m)
    end
    p = 1
    for _ in 1:10
        y = y / ((1+√(1-y))*(1+√(1-m*y)))
        p *= 2
        y < yA && return p*√y*serf(y, m)
    end
end

function acn(c, mc)
    m = one(1) - mc
    x = c^2
    p = one(c)
    for _ in 1:10
        if (x > 0.5) 
            return p*asn(√(1-x), m)
        end
        d = √(mc + m * x)
        x = (√x + d)/(1+d)
        p *= 2
    end
end

function _F(φ, m)
    m == 0 && return φ
    m == 1 && return atanh(sin(φ))
    signφ = sign(φ)
    φ = abs(φ)
    
    sin(φ)^2 ≤ 0.9 && return signφ*asn(sin(φ), m)

    mc = 1 - m

    c = sin(π/2-φ)
    z = c/√(mc+m*c^2)
    z^2 ≤ 0.9 && return signφ*(_K(m) - asn(z, m))

    w = √(1-z^2)
    c > w && return signφ*acn(c, m)
    return signφ*(K(m) - acn(w, m))
end

function ellipF(φ, m)
    abs(φ) < π/2 &&  return _F(φ, m)
    j = round(φ/π)
    return 2*j*_K(m) + _F(φ - j*π, m)
end