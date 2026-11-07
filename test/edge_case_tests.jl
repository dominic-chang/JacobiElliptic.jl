# Parameter regions the sweeps in default_tests.jl / carlson_tests.jl never reach: those
# sample m ∈ (1e-3, 0.999) and n ∈ (1e-3, 0.999) only. Every defect found while auditing
# the Reactant extension lived outside that box, so these pin the boundaries down.
#
# References come from ArbNumerics at 200 bits.

_arb(f, xs...) =
    try
        Float64(f(map(ArbNumerics.ArbFloat, xs)...))
    catch
        NaN
    end

@testset "Edge cases" begin
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
            @testset "φ = $φ" for φ in
                                  (1.01 * arch, 1.5 * arch, π / 2 + 0.1, 3.0, -3.0, 5.0)
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

    @testset "FukushimaAlg: ported fixes" begin
        # The |φ| > π/2 reduction, the n == m closed form and the n == 1 pole were fixed
        # in CarlsonAlg first and then ported here; these guard the port.
        FA = JacobiElliptic.FukushimaAlg

        @testset "|φ| > π/2 with m < 0: $((φ, m))" for (φ, m) in (
            (2.0, -2.0),
            (3.0, -0.5),
            (-3.0, -2.0),
            (1.6, -2.0),
            (5.0, -0.5),
        )
            @test FA.Pi(0.5, φ, m) ≈ _arb(ArbNumerics.elliptic_pi, 0.5, φ, m) rtol = RTOL
        end

        @testset "n == m = $m" for m in (-2.0, -0.5, -0.1)
            @test FA.Pi(m, 0.5, m) ≈ _arb(ArbNumerics.elliptic_pi, m, 0.5, m) rtol = RTOL
        end

        @test isinf(FA.Pi(1.0, 0.5))          # the 1/(1 - n sin²θ) pole at θ = π/2
        @testset "Pi(n, 0) stays finite: n = $n" for n in (-2.0, -0.5, 0.0, 0.5, 0.9)
            @test FA.Pi(n, 0.0) ≈ π / (2 * sqrt(1 - n)) rtol = RTOL
        end
    end

    @testset "incomplete E at large negative m" begin
        # `_Ds` used to stop its descent too early for m < 0, leaving `_Dsi`'s series to
        # diverge: E(0.7, -100) came back as 2.35e12 against a true 2.51, and Float32 was
        # off by 4e31. The descent cut-off now scales with 1 - m.
        FA = JacobiElliptic.FukushimaAlg

        @testset "φ = $φ, m = $m" for φ in (0.3, 0.7, 1.2, 1.5, 3.0, -0.7),
            m in (-5.0, -10.0, -50.0, -100.0, -1.0e4, -1.0e8, -1.0e15)

            @test FA.E(φ, m) ≈ _arb(ArbNumerics.elliptic_e, φ, m) rtol = 1e-11
        end

        # Float32 stops at |m| = 1e6, where `_Dsi` overflows: it evaluates its m-polynomial
        # by Horner before scaling by s^2k, so the m^7 coefficient passes 3.4e38. That
        # limit predates this fix -- Float32 returned NaN there before it too.
        @testset "Float32: φ = $φ, m = $m" for φ in (0.3f0, 0.7f0, 1.5f0),
            m in (-10.0f0, -100.0f0, -1.0f3, -1.0f4, -1.0f5)

            @test FA.E(φ, m) ≈ _arb(ArbNumerics.elliptic_e, φ, m) rtol = 1.0f-5
        end
    end
end

@testset "GPU-hostile error paths" begin
    # Neither of these is reachable through the public API on the CPU -- the argument
    # reduction in F/K/ellipke keeps the base helpers in range, and `_am`'s AGM counter is
    # never negative. They matter because GPU backends compile the dead branches anyway,
    # and Metal faults on throws it cannot lower (an opaque device exception, no name).

    @testset "_am raises 2^n in float, not Int" begin
        # An Int base sends `2^n` through `Base.power_by_squaring`, which carries a
        # `DomainError` branch for negative exponents. `n` only counts up from 0, so the
        # branch is dead -- but it still compiles into GPU kernels, and Metal traps on it
        # (`throw_domerr_powbysq at kernels.jl:49` in a device backtrace). A float base has
        # no such branch. Checked on the unoptimized typed IR: in a small kernel inlining
        # constant-folds the branch away and hides the regression.
        @testset "$alg, $T" for alg in
                                (JacobiElliptic.CarlsonAlg, JacobiElliptic.FukushimaAlg),
            T in (Float64, Float32)

            ir = sprint(
                io -> InteractiveUtils.code_warntype(
                    io,
                    alg._am,
                    (T, T, T);
                    optimize = false,
                ),
            )
            powers = filter(l -> occursin(" ^ ", l), split(ir, '\n'))
            @test !isempty(powers)                        # the probe still finds the call
            @test all(l -> endswith(l, "::$T"), powers)   # ...and it is raised in T
        end
    end

    @testset "Carlson base kernels return NaN out of domain" begin
        # DRF/DRFD report out-of-domain input with a nonzero error code rather than
        # throwing, so the base kernels must turn that into NaN. The `@assert`s these
        # replace turned it into an `AssertionError`, which a GPU cannot lower -- and a
        # SIMD group runs both arms of a divergent branch, so lanes the argument reduction
        # already rejected still reach these calls and tripped the assert for the group.
        CA = JacobiElliptic.CarlsonAlg

        @test isnan(CA._rawF(2.0, 0.5))                 # |sinφ| > 1
        @test isnan(CA._rawF(0.5f0, 5.0f0))             # 1 - m sin²φ < 0
        @test isnan(first(CA._K_base(2.0)))             # m > 1
        @test isnan(first(CA._K_base(2.0f0)))
        @test all(isnan, CA._ellipke_base(2.0))
        @test all(isnan, CA._ellipke_base(2.0f0))

        # and the in-domain results are untouched
        @test CA._rawF(sin(0.5), 0.5) ≈ _arb(ArbNumerics.elliptic_f, 0.5, 0.5) rtol = 1e-11
        @test first(CA._K_base(0.5)) ≈ _arb(ArbNumerics.elliptic_k, 0.5) rtol = 1e-11
    end
    @testset "sn, cn for m near 1 at |u| up to 2K(m)" begin
        # A&S 16.15.4 (first order in 1 - m) holds only for (1 - m) cosh(u)^2 << 1; the
        # cases reach u ≈ 1.9 K(m). (u, m, sn, cn) from mpmath at 50 digits.
        cases = (
            (5.8739636482105215, 0.999999999, 0.9999841889869129, 0.005623324300283354),
            (17.621890944631566, 0.999999999, 0.9999841889869129, -0.0056233243002833545),
            (22.32106186319998, 0.999999999, 0.82580254414494, -0.5639593585390212),
            (22.80272397010518, 0.999999999999, 0.9999995000059054, -0.0009999939695327612),
            (28.883450362133228, 0.999999999999, 0.9087293000336927, -0.4173859835455366),
            (27.984123535935506, 0.999999999999999, 0.9999999841949322, -0.00017779239405545436),
            (35.44655647885164, 0.999999999999999, 0.9531943021914329, -0.30235843343585994),
            (-19.314222997786473, 0.9999999917019833, -0.9683816763736394, -0.24947330290790634),
        )
        @testset "u = $u, m = $m" for (u, m, sn_ref, cn_ref) in cases
            @test JacobiElliptic.sn(u, m) ≈ sn_ref atol = 1e-12
            @test JacobiElliptic.cn(u, m) ≈ cn_ref atol = 1e-12
        end
        @test JacobiElliptic.cn(3.0, 1.0) ≈ sech(3.0) rtol = 1e-14     # m = 1: am = gd(u)
        @test JacobiElliptic.sn(40.0, 1.0) ≈ tanh(40.0) rtol = 1e-14
    end
end
