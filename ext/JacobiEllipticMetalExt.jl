module JacobiEllipticMetalExt

using JacobiElliptic, Metal

for alg in [JacobiElliptic.CarlsonAlg, JacobiElliptic.FukushimaAlg]
    @eval function ($alg)._am(u::A, m::B, tol::C) where {A,B,C}

        T = promote_type(A, B, C)

        # Use an immutable local buffer so GPU backends can keep this in registers.
        _ambuf = JacobiElliptic.StaticArrays.@SVector[
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
        ]
        u == 0 && return zero(T)

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
