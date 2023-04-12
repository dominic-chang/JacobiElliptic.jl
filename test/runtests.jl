using Test
using FastElliptic
using ArbNumerics

ArbNumerics.elliptic_k(m::Float64) = ArbNumerics.elliptic_k(ArbNumerics.ArbFloat(m))
ArbNumerics.elliptic_e(m::Float64) = ArbNumerics.elliptic_e(ArbNumerics.ArbFloat(m))
ArbNumerics.elliptic_f(φ::Float64, m::Float64) = ArbNumerics.elliptic_f(ArbNumerics.ArbFloat(φ), ArbNumerics.ArbFloat(m))


#Test EllipticK
@test begin
    all([isapprox(FastElliptic.K(m), ArbNumerics.elliptic_k(m), rtol=1e-6) for m in range(0.0,1.0,length=1000)])
end

#Test EllipticE
@test begin
    all([isapprox(FastElliptic.E(m), ArbNumerics.elliptic_e(m), rtol=1e-6) for m in range(0.0,1.0,length=1000)])
end


#Test EllipticF
@testset begin
    #test standard domain
    @test begin
        ans = []

        for φ in range(-π/2, π/2, length=1000)
            temp = true
            for m in range(1e-3, 1-1e-3,length=1000)
                temp = isapprox(FastElliptic.F(φ, m)/ArbNumerics.elliptic_f(φ, m), 1.0, atol=1e-6)
                append!(ans, temp)
                !temp && break
            end
            !temp && break
        end
        all(ans)
    end
    #test domain for abs(φ) > π/2
    @test begin
        ans = []
        for φ in range(-10π/2, 10π/2, length=10)
            temp = true
            for m in range(1e-3, 1-1e-3,length=10)
                temp = isapprox(FastElliptic.F(φ, m)/ArbNumerics.elliptic_f(φ, m), 1.0, atol=1e-6)
                append!(ans, temp)
                !temp && break
            end
        !temp && break
        end
        all(ans)
    end

    # test domain for |msin(φ)| ≤ 1 && |m| > 1
    @test begin
        ans = []
        for m in range(-100, 100,length=10)
            temp = true
            mmax = asin(1/√abs(m))
            for i in -1.0:0.10001:1.0
                φ = i*mmax
                temp = isapprox(FastElliptic.F(φ, m)/ArbNumerics.elliptic_f(φ, m), 1.0, atol=1e-6)
                append!(ans, temp)
                !temp && break
            end
            !temp && break
        end
        all(ans)
    end
end