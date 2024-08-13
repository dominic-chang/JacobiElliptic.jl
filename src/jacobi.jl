#module Jacobi
import StaticArrays

export am,
    sn, cn, dn, nn,
    sd, cd, dd, nd,
    sc, cc, dc, nc,
    ss, cs, ds, ns

# Abramowitz & Stegun, section 16.4, p571
#const _ambuf = Array{Float64}(undef, 10)
function _am(u::A, m::B, tol::C) where {A,B,C}

    T = promote_type(A, B, C)

    #Pre-allocate with SVector to avoid GPU allocations
    _ambuf = StaticArrays.@MVector[zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T)]
    u == 0 && return zero(T)

    sqrt_tol = sqrt(tol)
    if m < sqrt_tol
        # A&S 16.13.4
        return u - m*(u - sin(2*u)/2)/4
    end
    m1 = 1 - m
    if m1 < sqrt_tol
        # A&S 16.15.4
        t = tanh(u)
        return asin(t) + m1*(t - u*(1 - t^2))*cosh(u)/4
    end

    a,b,c,n = 1, sqrt(m1), sqrt(m), 0
    while abs(c) > tol
        @assert n < 10
        a,b,c,n = ((a+b)/2, sqrt(a*b), (a-b)/2, n+1)
        @inbounds _ambuf[n] = c/a
    end

    #phi = ldexp(a*u, n) # Doesn't work with Enzyme
    phi = a*u*(2^n)
    for i = n:-1:1
        @inbounds phi = (phi + asin(_ambuf[i]*sin(phi)))/2
    end
    phi
end
function _am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    return _am(u, m, eps(T))
end

"""
    am(u::Real, m::Real, [tol::Real=eps(Float64)])

Returns amplitude, φ, such that u = F(φ | m)

Landen sequence with convergence to `tol` used if `√(tol) ≤ m ≤ 1 - √(tol)`
"""
function am(u::A, m::B, tol::C) where {A, B, C}
    !(0 ≤ m ≤ 1) && throw(DomainError(m, "argument m not in [0,1]"))
    return _am(u, m, tol)
end
function am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    return am(u, m, eps(T))
end
#am(u::Real, m::Real) = am(Float64(u), Float64(m))

for (f,a,b,c) in ((:sn, :(sin(phi)),                :(sqrtmu1*s), :(sqrt(mu)*sin(phi))),
                  (:cn, :(cos(phi)),                :(cos(phi)),  :(sqrt(1 - mu*sin(phi)^2))),
                  (:dn, :(sqrt(1 - m*sin(phi)^2)), :(1),        :(cos(phi))))
    @eval begin
        function ($f)(u::A, m::B) where {A, B}
            T = promote_type(A, B)
            # Abramowitz & Stegun, section 16.10, p573
            lt0 = m < 0
            gt1 = m > 1
            if !(lt0 || gt1)
                phi = _am(u,m)
                return $a
            elseif lt0
                mu1 = inv(1 - m)
                mu = -m*mu1
                sqrtmu1 = sqrt(mu1)
                v = u/sqrtmu1
                phi = _am(v,mu)
                s = sin(phi)
                return ($b)/sqrt(1 - mu*s^2)
            elseif gt1
                mu = inv(m)
                v = u*sqrt(m)
                phi = _am(v,mu)
                return $c
            end
            return T(NaN)
        end
    end
end

xn = ((:s,:(sn(u,m))), (:c,:(cn(u,m))), (:d,:(dn(u,m))), (:n,:(1)))
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

#end # module
