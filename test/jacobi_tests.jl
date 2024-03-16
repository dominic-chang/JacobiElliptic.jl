#
# Tests for JacobiElliptic Functions
#

@testset "JacobiElliptic" begin
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

                @test JacobiElliptic.sn(u, m) ≈ θs / θn atol=2.5e-9
                @test JacobiElliptic.cn(u, m) ≈ θc / θn atol=2.5e-9
                @test JacobiElliptic.dn(u, m) ≈ θd / θn atol=2.5e-9
                @test JacobiElliptic.nn(u, m) == 1.0

                @test JacobiElliptic.sd(u, m) ≈ θs / θd atol=2.5e-9
                @test JacobiElliptic.cd(u, m) ≈ θc / θd atol=2.5e-9
                @test JacobiElliptic.dd(u, m) == 1.0
                @test JacobiElliptic.nd(u, m) ≈ θn / θd atol=2.5e-9

                @test JacobiElliptic.cc(u, m) == 1.0
                if ϵ != 90
                    # very sensitive around u = K,
                    # estimate of K(0) = pi/2 + 4e-9, so cosine causes errors
                    # also, errors build up in ϵ ≥75°, so lower tolerences for that
                    @test JacobiElliptic.sc(u, m) ≈ θs / θc atol=1e-7
                    @test JacobiElliptic.dc(u, m) ≈ θd / θc atol=1e-7
                    @test JacobiElliptic.nc(u, m) ≈ θn / θc atol=1e-7
                end

                @test JacobiElliptic.ss(u, m) == 1.0
                @test JacobiElliptic.cs(u, m) ≈ θc / θs atol=5.5e-8
                @test JacobiElliptic.ds(u, m) ≈ θd / θs atol=5.5e-8
                @test JacobiElliptic.ns(u, m) ≈ θn / θs atol=5.5e-8

                # ellipj
                s, c, d = JacobiElliptic.CarlsonAlg.ellipj(u, m)
                @test s ≈ θs / θn atol=1e-9
                @test c ≈ θc / θn atol=1e-9
                @test d ≈ θd / θn atol=1e-9
            end
        end
    end

    @testset "u = $u" for u in -1.:0.21:1.0
        @test JacobiElliptic.am(u,0.0) ≈ u
        @test JacobiElliptic.sn(u,0.0) ≈ sin(u)
        @test JacobiElliptic.cn(u,0.0) ≈ cos(u)
        @test JacobiElliptic.dn(u,0.0) ≈ 1.
        @test JacobiElliptic.nn(u,0.0) ≈ 1.
        @test JacobiElliptic.cd(u,0.0) ≈ cos(u)
        @test JacobiElliptic.sd(u,0.0) ≈ sin(u)
        @test JacobiElliptic.nd(u,0.0) ≈ 1.
        @test JacobiElliptic.dd(u,0.0) ≈ 1.
        @test JacobiElliptic.dc(u,0.0) ≈ sec(u)
        @test JacobiElliptic.nc(u,0.0) ≈ sec(u)
        @test JacobiElliptic.sc(u,0.0) ≈ tan(u)
        @test JacobiElliptic.cc(u,0.0) ≈ 1.
        @test JacobiElliptic.ns(u,0.0) ≈ csc(u)
        @test JacobiElliptic.ds(u,0.0) ≈ csc(u)
        @test JacobiElliptic.cs(u,0.0) ≈ cot(u)
        @test JacobiElliptic.ss(u,0.0) ≈ 1.

        @test JacobiElliptic.am(u,1.0) ≈ atan(sinh(u))
        @test JacobiElliptic.sn(u,1.0) ≈ tanh(u)
        @test JacobiElliptic.cn(u,1.0) ≈ sech(u)
        @test JacobiElliptic.dn(u,1.0) ≈ sech(u)
        @test JacobiElliptic.nn(u,1.0) ≈ 1.
        @test JacobiElliptic.cd(u,1.0) ≈ 1.
        @test JacobiElliptic.sd(u,1.0) ≈ sinh(u)
        @test JacobiElliptic.nd(u,1.0) ≈ cosh(u)
        @test JacobiElliptic.dd(u,1.0) ≈ 1.
        @test JacobiElliptic.dc(u,1.0) ≈ 1.
        @test JacobiElliptic.nc(u,1.0) ≈ cosh(u)
        @test JacobiElliptic.sc(u,1.0) ≈ sinh(u)
        @test JacobiElliptic.cc(u,1.0) ≈ 1.
        @test JacobiElliptic.ns(u,1.0) ≈ coth(u)
        @test JacobiElliptic.ds(u,1.0) ≈ csch(u)
        @test JacobiElliptic.cs(u,1.0) ≈ csch(u)
        @test JacobiElliptic.ss(u,1.0) ≈ 1.
    end

    @testset "errors" begin
        u = 0.5 # random value for testing
        @test_throws DomainError JacobiElliptic.am(u, -0.1)
        @test_throws DomainError JacobiElliptic.am(u, -eps(0.0))

        @test_throws DomainError JacobiElliptic.am(u, 1 + eps(1.0))
        @test_throws DomainError JacobiElliptic.am(u, 1.1)
    end
end
