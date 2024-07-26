using Zygote
using ForwardDiff
using Enzyme
using SpecialFunctions

@testset "alg:$alg" for alg in [JacobiElliptic.CarlsonAlg, ]# JacobiElliptic.FukushimaAlg]
    @testset "Zygote and ForwardDiff" begin

        num_trials = 1
        ks = rand(num_trials)
        ϕs = rand(num_trials) .* 2π
        ns = rand(num_trials)

        @show ks, ϕs, ns

        # Test several known derivative identities across a wide range of values of ϕ and m
        # in order to verify that derivatives work correctly
        @testset for (k, ϕ, n) in zip(ks, ϕs, ns)
            m = k^2
            dk_dm = inv(2k)

            # I. Tests for complete integrals, from https://en.wikipedia.org/wiki/Elliptic_integral
            # 1. K'(m) == K'(k) * dk/dm == E(k) / (k * (1 - k^2)) - K(k)/k
            @test Zygote.gradient(alg.K, m)[1] ≈ (alg.E(m) / (k * (1 - k^2)) - alg.K(m) / k) * dk_dm
            @test ForwardDiff.derivative(alg.K, m) ≈ (alg.E(m) / (k * (1 - k^2)) - alg.K(m) / k) * dk_dm
            @test Enzyme.autodiff(Reverse, alg.K, Active, Active(m))[1][1] ≈ (alg.E(m) / (k * (1 - k^2)) - alg.K(m) / k) * dk_dm

            # 2. E'(m) == E'(k) * dk/dm == (E(k) - K(k))/k
            @test Zygote.gradient(alg.E, m)[1] ≈ (alg.E(m) - alg.K(m))/k * dk_dm
            @test ForwardDiff.derivative(alg.E, m) ≈ (alg.E(m) - alg.K(m))/k * dk_dm
            @test Enzyme.autodiff(Reverse, alg.E, Active, Active(m))[1][1] ≈ (alg.E(m) - alg.K(m))/k * dk_dm

            # II. Tests for incomplete integrals, from https://functions.wolfram.com/EllipticIntegrals/EllipticF/introductions/IncompleteEllipticIntegrals/ShowAll.html
            # 3. ∂ϕ(F(ϕ, m)) == 1 / √(1 - m*sin(ϕ)^2)
            _F = alg.F
            @test Zygote.gradient(ϕ -> _F(ϕ, m), ϕ)[1] ≈ 1 / √(1 - m*sin(ϕ)^2) atol=1e-5
            @test ForwardDiff.derivative(ϕ -> _F(ϕ, m), ϕ) ≈ 1 / √(1 - m*sin(ϕ)^2) atol=1e-5
            @test Enzyme.autodiff(Reverse, ϕ -> _F(ϕ, m), Active, Active(ϕ))[1][1] ≈ 1 / √(1 - m*sin(ϕ)^2) atol=1e-5

            # 4. ∂m(F(ϕ, m)) == E(ϕ, m) / (2 * m * (1 - m)) - F(ϕ, m) / 2m - sin(2ϕ) / (4 * (1-m) * √(1 - m * sin(ϕ)^2))
            @test Zygote.gradient(m -> _F(ϕ, m), m)[1] ≈
                alg.E(ϕ, m) / (2 * m * (1 - m)) -
                alg.F(ϕ, m) / 2 / m -
                sin(2*ϕ) / (4 * (1 - m) * √(1 - m * sin(ϕ)^2)) atol=1e-5
            @test ForwardDiff.derivative(m -> _F(ϕ, m), m) ≈
                alg.E(ϕ, m) / (2 * m * (1 - m)) -
                alg.F(ϕ, m) / 2 / m -
                sin(2*ϕ) / (4 * (1 - m) * √(1 - m * sin(ϕ)^2)) atol=1e-5
            @test Enzyme.autodiff(Reverse, m -> _F(ϕ, m), Active, Active(m))[1][1] ≈
                alg.E(ϕ, m) / (2 * m * (1 - m)) -
                alg.F(ϕ, m) / 2 / m -
                sin(2*ϕ) / (4 * (1 - m) * √(1 - m * sin(ϕ)^2))  atol=1e-5

            _E = alg.E
            # 5. ∂ϕ(E(ϕ, m)) == √(1 - m * sin(ϕ)^2)
            @test Zygote.gradient(ϕ -> _E(ϕ, m), ϕ)[1] ≈ √(1 - m * sin(ϕ)^2) atol=1e-5
            @test ForwardDiff.derivative(ϕ -> _E(ϕ, m), ϕ) ≈ √(1 - m * sin(ϕ)^2) atol=1e-5
            @test Enzyme.autodiff(Reverse, ϕ -> _E(ϕ, m), Active, Active(ϕ))[1][1] ≈ √(1 - m * sin(ϕ)^2) atol=1e-5

            # 6. ∂m(E(ϕ, m)) == (E(ϕ, m) - F(ϕ, m)) / 2m
            @test Zygote.gradient(m -> _E(ϕ, m), m)[1] ≈ (alg.E(ϕ, m) - alg.F(ϕ, m)) / 2m
            @test ForwardDiff.derivative(m -> _E(ϕ, m), m) ≈ (alg.E(ϕ, m) - alg.F(ϕ, m)) / 2m
            @test Enzyme.autodiff(Reverse, m -> _E(ϕ, m), Active, Active(m))[1][1] ≈ (alg.E(ϕ, m) - alg.F(ϕ, m)) / 2m atol=1e-5

        end
    end

    @testset "SpecialCases" begin
        _E = alg.E
        @test ForwardDiff.gradient(x -> _E(x[1], x[2]), [π/2.0, 0.0]) |> collect ≈  [1.0, -π/8]
        @test Enzyme.autodiff(Reverse, _E, Active, Active(π/2.0), Active(0.0))[1] |> collect ≈  [1.0, -π/8]
        @test Zygote.gradient(x -> _E(x[1], x[2]), [π/2.0, 0.0])[1] ≈  [1.0, -π/8]

    end

end