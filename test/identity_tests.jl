# Exactness and identity assertions. The ratio-against-ArbNumerics tests elsewhere compare
# a value to a reference; these constrain relationships *between* the implementation's own
# outputs, which catches cancellation and branch-consistency bugs a reference comparison
# cannot -- mixing AGM's K into a Carlson identity leaves F(0, m) at ~2e-16 instead of
# exactly 0, while still matching Arb to 1e-15.

@testset "Identities" begin
    MS_OK = (-100.0, -2.0, -0.5, -1e-8, 0.0, 1e-8, 0.01, 0.3, 0.5, 0.7, 0.9, 0.99)
    US = (0.1, 0.4, 0.8, 1.2, 1.8, 2.4, -0.7, -2.4)
    alg = JacobiElliptic.CarlsonAlg          # the default algorithm

    @testset "exact at φ = 0" begin
        @testset "m = $m" for m in MS_OK
            @test alg.F(0.0, m) == 0.0
            @test alg.E(0.0, m) == 0.0
            @test alg.Pi(0.5, 0.0, m) == 0.0
        end
    end

    @testset "complete integrals are the φ = π/2 incomplete ones" begin
        @testset "m = $m" for m in MS_OK
            @test alg.F(π / 2, m) ≈ alg.K(m) rtol = 1e-14
            @test alg.E(π / 2, m) ≈ alg.E(m) rtol = 1e-14
        end
    end

    @testset "Pythagorean identities" begin
        @testset "u = $u, m = $m" for u in US, m in MS_OK
            s, c, d = alg.sn(u, m), alg.cn(u, m), alg.dn(u, m)
            @test s^2 + c^2 ≈ 1 atol = 1e-13
            @test d^2 + m * s^2 ≈ 1 atol = 1e-13
        end
    end

    @testset "Glaisher quotients agree with their definitions" begin
        @testset "u = $u, m = $m" for u in US, m in (0.1, 0.3, 0.7, 0.9)
            s, c, d = alg.sn(u, m), alg.cn(u, m), alg.dn(u, m)
            @test alg.sd(u, m) ≈ s / d rtol = 1e-13
            @test alg.cd(u, m) ≈ c / d rtol = 1e-13
            @test alg.nd(u, m) ≈ 1 / d rtol = 1e-13
            @test alg.sc(u, m) ≈ s / c rtol = 1e-13
            @test alg.dc(u, m) ≈ d / c rtol = 1e-13
            @test alg.nc(u, m) ≈ 1 / c rtol = 1e-13
            @test alg.cs(u, m) ≈ c / s rtol = 1e-13
            @test alg.ds(u, m) ≈ d / s rtol = 1e-13
            @test alg.ns(u, m) ≈ 1 / s rtol = 1e-13
            @test alg.nn(u, m) == 1
            @test alg.ss(u, m) == 1
            @test alg.cc(u, m) == 1
            @test alg.dd(u, m) == 1
        end
    end

    @testset "am inverts F" begin
        @testset "φ = $φ, m = $m" for φ in (-1.4, -0.9, -0.3, 0.0, 0.3, 0.9, 1.4),
            m in (-2.0, -0.5, 0.0, 0.3, 0.7, 0.9, 0.99)

            @test alg.am(alg.F(φ, m), m) ≈ φ atol = 1e-13
        end
    end

    @testset "odd and even symmetry" begin
        @testset "u = $u, m = $m" for u in (0.3, 0.8, 1.6, 2.4), m in (-0.5, 0.3, 0.9)
            @test alg.sn(-u, m) == -alg.sn(u, m)
            @test alg.cn(-u, m) == alg.cn(u, m)
            @test alg.dn(-u, m) == alg.dn(u, m)
            # F and E are exactly odd only on the direct branch, where they are computed as
            # sign(φ) * _F(|φ|, m). The m < 0 transformation (newφ = π/2 - φ) and the
            # |φ| > π/2 reduction are not symmetric in φ, so those agree to ~1 ulp instead.
            if m >= 0 && abs(u) <= π / 2
                @test alg.F(-u, m) == -alg.F(u, m)
                @test alg.E(-u, m) == -alg.E(u, m)
            else
                @test alg.F(-u, m) ≈ -alg.F(u, m) rtol = 1e-14
                @test alg.E(-u, m) ≈ -alg.E(u, m) rtol = 1e-14
            end
        end
    end

    @testset "J relates to Pi and F" begin
        # J(n; φ | m) = (Π(n; φ | m) - F(φ | m)) / n
        @testset "n = $n, m = $m" for n in (0.1, 0.5, 0.9), m in (0.1, 0.5, 0.9)
            φ = 0.7
            @test alg.J(n, φ, m) ≈ (alg.Pi(n, φ, m) - alg.F(φ, m)) / n rtol = 1e-10
        end
    end

    @testset "quarter period: am(K) == π/2" begin
        @testset "m = $m" for m in (0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
            K = JacobiElliptic.K(m)
            @test JacobiElliptic.am(K, m) ≈ π / 2 atol = 1e-12
            @test JacobiElliptic.sn(K, m) ≈ 1 atol = 1e-12
            @test JacobiElliptic.cn(K, m) ≈ 0 atol = 1e-12
            @test JacobiElliptic.dn(K, m) ≈ sqrt(1 - m) atol = 1e-12
        end
    end

    @testset "documented domain behaviour" begin
        # F(φ | m > 1) is real only on the first arch |φ| ≤ asin(1/√m); outside it the
        # result is NaN rather than a throw or a wrong-branch value.
        @test isnan(JacobiElliptic.CarlsonAlg.F(1.2, 1.5))
        @test isnan(JacobiElliptic.CarlsonAlg.F(1.0, 2.0))
        @test isnan(JacobiElliptic.CarlsonAlg.F(3.0, 1.5))
        # inside the arch it is an ordinary real value
        @test JacobiElliptic.CarlsonAlg.F(0.5, 1.5) isa Float64
        @test !isnan(JacobiElliptic.CarlsonAlg.F(0.5, 1.5))
    end
end
