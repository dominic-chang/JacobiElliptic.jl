# Parameter regions the sweeps in default_tests.jl / carlson_tests.jl never reach: those
# sample m ∈ (1e-3, 0.999) and n ∈ (1e-3, 0.999) only. Every defect found while auditing
# the Reactant extension lived outside that box, so these pin the boundaries down.
#
# References come from ArbNumerics at 200 bits.

_arb(f, xs...) = try
    Float64(f(map(ArbNumerics.ArbFloat, xs)...))
catch
    NaN
end

_try(f, xs...) = try
    f(xs...)
catch
    NaN
end

@testset "Edge cases" begin
    setprecision(ArbNumerics.ArbFloat, 200)
    RTOL = 1e-11
    alg = JacobiElliptic.CarlsonAlg      # the default algorithm for all of these

    @testset "m < 0" begin
        # the imaginary-modulus transformation, including |φ| > π/2, which must be reduced
        # before the transformation is applied (its asin cannot represent the outer branches)
        @testset "φ = $φ, m = $m" for φ in (0.0, 0.7, 1.5, 2.0, 3.0, -1.6, -3.0),
            m in (-100.0, -2.0, -0.5, -1e-8)

            @test alg.F(φ, m) ≈ _arb(ArbNumerics.elliptic_f, φ, m) rtol = RTOL
            @test alg.E(φ, m) ≈ _arb(ArbNumerics.elliptic_e, φ, m) rtol = RTOL
            @testset "n = $n" for n in (-2.0, 0.0, 0.5, 0.9)
                @test alg.Pi(n, φ, m) ≈ _arb(ArbNumerics.elliptic_pi, n, φ, m) rtol = RTOL
            end
        end
        @test alg.F(0.0, -0.5) == 0.0          # exact, by symmetry
    end

    @testset "n == m < 0 (0/0 in the transformation)" begin
        @testset "n = m = $m, φ = $φ" for m in (-100.0, -2.0, -0.5, -0.1),
            φ in (0.5, 1.0, 1.6, -1.2)

            @test alg.Pi(m, φ, m) ≈ _arb(ArbNumerics.elliptic_pi, m, φ, m) rtol = RTOL
        end
    end

    @testset "m == 0 and m == 1 exactly" begin
        # Π(n|0) = π/(2√(1-n)) is finite -- only m == 1 diverges
        @testset "n = $n" for n in (-2.0, -0.5, 0.0, 0.1, 0.5, 0.9)
            @test alg.Pi(n, 0.0) ≈ π / (2 * sqrt(1 - n)) rtol = RTOL
        end
        @test isinf(alg.Pi(0.5, 1.0))
        @test isinf(alg.K(1.0))
        @test alg.E(1.0) == 1.0
        # at m == 1 the Jacobi functions degenerate to the hyperbolic ones
        @testset "u = $u" for u in (0.0, 0.3, 0.8, 1.5, 2.4, 5.0)
            @test alg.am(u, 1.0) ≈ asin(tanh(u)) rtol = 1e-12
            @test alg.sn(u, 1.0) ≈ tanh(u) rtol = 1e-12
            @test alg.dn(u, 1.0) ≈ sech(u) rtol = 1e-12
        end
    end

    @testset "m > 1, within the real domain |φ| ≤ asin(1/√m)" begin
        @testset "m = $m" for m in (1.0001, 1.1, 1.5, 2.0, 10.0)
            arch = asin(1 / sqrt(m))
            @testset "φ = $φ" for φ in (0.0, 0.3 * arch, 0.7 * arch, 0.99 * arch)
                @test alg.F(φ, m) ≈ _arb(ArbNumerics.elliptic_f, φ, m) rtol = RTOL
                @test alg.E(φ, m) ≈ _arb(ArbNumerics.elliptic_e, φ, m) rtol = RTOL
            end
        end
    end

    @testset "n outside (0, 1)" begin
        @testset "n = $n, m = $m" for n in (-2.0, -0.5, 0.0, 1.0001, 1.7, 3.0),
            m in (0.1, 0.5, 0.9)

            @test alg.Pi(n, 0.5, m) ≈ _arb(ArbNumerics.elliptic_pi, n, 0.5, m) rtol = RTOL
        end
        @test isinf(alg.Pi(1.0, 0.5))   # the 1/(1 - n sin²θ) pole at θ = π/2
    end

    # ----------------------------------------------------------------------------------
    # Known defects, reproducible and marked broken so that a fix flips them to Pass
    # rather than going unnoticed.
    # ----------------------------------------------------------------------------------
    @testset "m > 1 outside the real arch returns NaN" begin
        # F(φ|m) and Π(n;φ|m) are real only on |φ| ≤ asin(1/√m); beyond it the integrand
        # turns negative somewhere in [0, φ]. Previously this threw DomainError for some φ
        # and silently returned a wrong-branch value for others (F(3.0, 1.5) gave 0.1423).
        @testset "m = $m" for m in (1.0001, 1.1, 1.5, 2.0, 10.0)
            arch = asin(1 / sqrt(m))
            @testset "φ = $φ" for φ in (1.01 * arch, 1.5 * arch, π / 2 + 0.1, 3.0, -3.0, 5.0)
                @test isnan(alg.F(φ, m))
                @test isnan(alg.Pi(0.5, φ, m))
            end
        end
    end

    @testset "n > 1 and m > 1 (reciprocal modulus)" begin
        # n₁ = m/n is an exact 2-cycle when both exceed 1, so the n > 1 reduction used to
        # recurse forever. Π is now reduced through 1/m first.
        @testset "n = $n, φ = $φ, m = $m" for (n, φ, m) in (
            (1.0001, 0.5, 1.5),
            (1.7, 0.5, 2.0),
            (3.0, 0.1, 10.0),
            (1.3, 0.6, 1.2),
            (-2.0, 0.4, 2.0),
        )
            @test alg.Pi(n, φ, m) ≈ _arb(ArbNumerics.elliptic_pi, n, φ, m) rtol = RTOL
        end
    end

    @testset "known defects: FukushimaAlg" begin
        # FukushimaAlg still carries the three defects that were fixed in CarlsonAlg,
        # plus a large-|m| error in the incomplete E. It is never exercised outside
        # m ∈ (0, 1) by the other suites, which is why these went unnoticed.
        FA = JacobiElliptic.FukushimaAlg

        # 1. incomplete E at large |m|: E(0.7, -100) returns 2.35e12 for a true 2.51
        @test_broken isapprox(
            _try(FA.E, 0.7, -100.0),
            _arb(ArbNumerics.elliptic_e, 0.7, -100.0),
            rtol = 1e-11,
        )

        # 2. |φ| > π/2 with m < 0: the transformation is applied before the argument is
        #    reduced, so the half-period is lost (fixed in CarlsonAlg)
        @testset "|φ| > π/2, m < 0: $((φ, m))" for (φ, m) in
                                                   ((2.0, -2.0), (3.0, -0.5), (-3.0, -2.0))
            @test_broken isapprox(
                _try(FA.Pi, 0.5, φ, m),
                _arb(ArbNumerics.elliptic_pi, 0.5, φ, m),
                rtol = 1e-11,
            )
        end

        # 3. n == m < 0 is 0/0 in the transformation (fixed in CarlsonAlg)
        @testset "n == m = $m" for m in (-2.0, -0.5)
            @test_broken isapprox(
                _try(FA.Pi, m, 0.5, m),
                _arb(ArbNumerics.elliptic_pi, m, 0.5, m),
                rtol = 1e-11,
            )
        end

        # 4. n == 1 is a pole; CarlsonAlg returns Inf, FukushimaAlg returns a finite value
        @test_broken isinf(_try(FA.Pi, 1.0, 0.5))
    end
end
