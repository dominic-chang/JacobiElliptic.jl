module JacobiEllipticMetalExt

using JacobiElliptic, Metal

# Resolve the process-wide default while constructing a broadcast on the host. The
# resulting concrete callable can then be compiled into a Metal kernel.
for name in JacobiElliptic.func_syms
    f = getfield(JacobiElliptic, name)
    @eval function Base.Broadcast.broadcasted(
        style::Metal.MtlArrayStyle,
        ::typeof($f),
        args...,
    )
        if !isempty(args) &&
           args[1] isa Base.RefValue{<:JacobiElliptic.AbstractAlgorithm}
            return Base.Broadcast.Broadcasted(style, $f, args)
        end

        algorithm =
            JacobiElliptic._effective_default($f, JacobiElliptic.default_algorithm())
        call = JacobiElliptic.AlgorithmCall($f, algorithm)
        return Base.Broadcast.broadcasted(style, call, args...)
    end
end

for alg in [JacobiElliptic.CarlsonAlg, JacobiElliptic.FukushimaAlg]
    @eval function ($alg)._am(u::T, m::T, tol::T) where {T<:Float32}
        zeroT = zero(T)

        # Use an immutable local buffer so GPU backends can keep this in registers.
        _ambuf = JacobiElliptic.StaticArrays.@SVector[
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
        u == 0 && return zeroT

        sqrt_tol = sqrt(tol)
        if m < sqrt_tol
            # A&S 16.13.4
            return u - m * (u - sin(2 * u) / 2) / 4
        end
        m1 = one(T) - m
        if m1 < sqrt_tol
            # A&S 16.15.4
            t = tanh(u)
            return asin(t) + m1 * (t - u * (one(T) - t^2)) * cosh(u) / 4
        end

        a = one(T)
        b = sqrt(m1)
        c = sqrt(m)
        n = 0
        while abs(c) > tol && n < 10
            a, b, c, n = ((a + b) / 2, sqrt(a * b), (a - b) / 2, n + 1)
            @inbounds _ambuf = Base.setindex(_ambuf, c / a, n)
        end
        abs(c) > tol && return T(NaN)

        #phi = ldexp(a*u, n) # Was slower on my benchmarks
        phi = a * u * (2^n)
        for i = n:-1:1
            @inbounds phi = (phi + asin(_ambuf[i] * sin(phi))) / 2
        end
        phi
    end

end
end
