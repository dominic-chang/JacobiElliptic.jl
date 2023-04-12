using Test
using FastElliptic
using ArbNumerics
using DelimitedFiles: readdlm

ArbNumerics.elliptic_k(m::Float64) = ArbNumerics.elliptic_k(ArbNumerics.ArbFloat(m))
ArbNumerics.elliptic_e(m::Float64) = ArbNumerics.elliptic_e(ArbNumerics.ArbFloat(m))
ArbNumerics.elliptic_f(φ::Float64, m::Float64) = ArbNumerics.elliptic_f(ArbNumerics.ArbFloat(φ), ArbNumerics.ArbFloat(m))

#Test EllipticK
@testset "Elliptic K" begin
    @test all([isapprox(FastElliptic.K(m), ArbNumerics.elliptic_k(m), rtol=1e-6) for m in range(0.0,1.0,length=1000)])
end

#Test EllipticE
@testset "Elliptic E" begin
    @test all([isapprox(FastElliptic.E(m), ArbNumerics.elliptic_e(m), rtol=1e-6) for m in range(0.0,1.0,length=1000)])
end


#Test EllipticF
@testset "Elliptic F" begin
    @testset "Standard Domain" for φ in range(-π/2, π/2, length=100)
        for m in range(1e-3, 1-1e-3,length=100)
            @test FastElliptic.F(φ, m)/ArbNumerics.elliptic_f(φ, m) ≈ 1.0 atol=1e-6
        end
    end
    #abs(φ) > π/2: 
    @testset "φ = $φ" for φ in range(-10π/2, 10π/2, length=10)
        for m in range(1e-3, 1-1e-3,length=10)
            @test FastElliptic.F(φ, m)/ArbNumerics.elliptic_f(φ, m) ≈ 1.0 atol=1e-6
        end
    end
    #|m sin(φ)| ≤ 1 && |m| > 1
    @testset "m = $m" for m in range(-100, 100,length=10)
        mmax = asin(1/√abs(m))
        for i in -1.0:0.10001:1.0
            φ = i*mmax
            @test FastElliptic.F(φ, m)/ArbNumerics.elliptic_f(φ, m) ≈ 1.0 atol=1e-6
        end
    end
end

#Taken from Elliptic.jl
#https://github.com/nolta/Elliptic.jl/blob/master/test/jacobi_tests.jl
@testset "Jacobi Elliptic Trig" begin
    @testset "Abramowitz & Stegun Table 16.1" begin
        dataloc = "data/"
        # table 16.1
        t161, _ = readdlm(joinpath(dataloc, "table_16_1.csv"), ',', header=true)
        # table 17.2
        t172, _ = readdlm(joinpath(dataloc, "table_17_2.csv"), ',', header=true)

        K_lut = Dict(zip(t172[:, 1], t172[:, 2]))

        # data vary first by ϵ ∈ 0:5:90 then α ∈ 0:5:85
        αs = 0:5:85
        ϵs = 0:5:90
        θss = reshape(t161[:, 3], length(ϵs), length(αs))
        θns = reshape(t161[:, 4], length(ϵs), length(αs))

        @testset "α = $α" for (i, α) in enumerate(αs)
            K = K_lut[α]
            m = sind(α)^2
            denom = sqrt(secd(α))

            @testset "ϵ = $ϵ" for (j, ϵ) in enumerate(ϵs)
                j₁ = length(ϵs) - j + 1
                ϵ₁ = ϵs[j₁]
                u = ϵ * K / 90

                θs = θss[j, i]
                θn = θns[j, i]
                θc = θss[j₁, i]/denom
                θd = θns[j₁, i]/denom

                @test FastElliptic.sn(u, m) ≈ θs / θn atol=2.5e-9
                @test FastElliptic.cn(u, m) ≈ θc / θn atol=2.5e-9
                @test FastElliptic.dn(u, m) ≈ θd / θn atol=2.5e-9

                @test FastElliptic.sd(u, m) ≈ θs / θd atol=2.5e-9
                if ϵ != 90
                    # very sensitive around u = K,
                    # estimate of K(0) = pi/2 + 4e-9, so cosine causes errors
                    # also, errors build up in ϵ ≥75°, so lower tolerences for that
                    @test FastElliptic.sc(u, m) ≈ θs / θc atol=1e-7
                end
            end
        end
    end

   #@testset "u = $u" for u in -1.:0.21:1.0
   #    @test FastElliptic.sn(u,0) ≈ sin(u)
   #    @test FastElliptic.cn(u,0) ≈ cos(u)
   #    @test FastElliptic.dn(u,0) ≈ 1.
   #    @test FastElliptic.sd(u,0) ≈ sin(u)
   #    @test FastElliptic.sc(u,0) ≈ tan(u)

   #    @test FastElliptic.sn(u,1) ≈ tanh(u)
   #    @test FastElliptic.cn(u,1) ≈ sech(u)
   #    @test FastElliptic.dn(u,1) ≈ sech(u)
   #    @test FastElliptic.sd(u,1) ≈ sinh(u)
   #    @test FastElliptic.sc(u,1) ≈ sinh(u)
   #end

end