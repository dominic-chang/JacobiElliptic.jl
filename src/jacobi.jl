#module Jacobi
import StaticArrays
include(joinpath(@__DIR__, "common.jl"))

export sn, cn, dn, nn, sd, cd, dd, nd, sc, cc, dc, nc, ss, cs, ds, ns
for (f, a, b, c) in (
    (:sn, :(sin(phi)), :(sqrtmu1 * s), :(_sqrt(mu) * sin(phi))),
    (:cn, :(cos(phi)), :(cos(phi)), :(_sqrt(1 - mu * sin(phi)^2))),
    (:dn, :(_sqrt(1 - m * sin(phi)^2)), :(1), :(cos(phi))),
)
    @eval begin
        function ($f)(u::A, m::B) where {A,B}
            T = promote_type(A, B)
            # Abramowitz & Stegun, section 16.10, p573
            lt0 = m < 0
            gt1 = m > 1
            if !(lt0 || gt1)
                phi = _am(u, m)
                return $a
            elseif lt0
                mu1 = inv(1 - m)
                mu = -m * mu1
                sqrtmu1 = _sqrt(mu1)
                v = u / sqrtmu1
                phi = _am(v, mu)
                s = sin(phi)
                return ($b) / _sqrt(1 - mu * s^2)
            elseif gt1
                mu = inv(m)
                v = u * _sqrt(m)
                phi = _am(v, mu)
                return $c
            end
            return T(NaN)
        end
    end
end

xn = ((:s, :(sn(u, m))), (:c, :(cn(u, m))), (:d, :(dn(u, m))), (:n, :(1)))
for (p, num) in xn, (q, den) in xn
    f = Symbol(p, q)
    #@eval begin
    #    """
    #        $($f)(u::Real, m::Real)

    #    Compute the Jacobi elliptic function $($f)(u | m)
    #    """
    #    ($f)(u::Real, m::Real) = ($f)(Float64(u), Float64(m))
    #end

    if (p == q)
        @eval ($f)(::A, ::B) where {A,B} = begin
            one(promote_type(A, B))
        end
    elseif (q != :n)
        @eval ($f)(u::A, m::B) where {A,B} = begin
            ($num) / ($den)
        end
    end
end

#end # module
