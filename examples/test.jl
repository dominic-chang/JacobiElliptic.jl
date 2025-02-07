using Revise
using JacobiElliptic
using BenchmarkTools
#@profview begin
#    for _ in 1:10_000_000
#        T = Float64
#        JacobiElliptic.Pi(JacobiElliptic.Fukushima(), rand(T), rand(T), rand(T))
#    end
#    nothing
#end
#
#@profview begin
#    for _ in 1:10_000_000
#        T = Float64
#        JacobiElliptic.Pi(JacobiElliptic.Carlson(), rand(T), rand(T), rand(T))
#    end
#    nothing
#end
#
@benchmark JacobiElliptic.Pi(JacobiElliptic.Fukushima(), rand(), rand(), rand())
@benchmark JacobiElliptic.Pi(JacobiElliptic.Carlson(), rand(), rand(), rand())

@benchmark JacobiElliptic.K(JacobiElliptic.Fukushima(), rand())
@benchmark JacobiElliptic.K(JacobiElliptic.Carlson(), rand())

@benchmark JacobiElliptic.F(JacobiElliptic.Fukushima(), rand() * 0.1, rand() + 1)
@benchmark JacobiElliptic.F(JacobiElliptic.Carlson(), rand() * 0.1, rand() + 1)

#@profview begin
#    for _ in 1:10_000_000
#        JacobiElliptic.F(JacobiElliptic.Fukushima(), rand(), rand())
#    end
#    nothing
#end
#

@benchmark JacobiElliptic.E(JacobiElliptic.Fukushima(), rand())
@benchmark JacobiElliptic.E(JacobiElliptic.Carlson(), rand())

#@profview begin
#    for _ in 1:10_000_000
#        JacobiElliptic.E(JacobiElliptic.Fukushima(), rand(Float32), rand(Float32))
#    end
#    nothing
#end

@btime JacobiElliptic.E(JacobiElliptic.Fukushima(), rand(Float32), rand(Float32))

@benchmark JacobiElliptic.E(JacobiElliptic.Fukushima(), rand(Float32), rand(Float32))
@benchmark JacobiElliptic.E(JacobiElliptic.Carlson(), rand(Float32), rand(Float32))

@btime JacobiElliptic.K(JacobiElliptic.Fukushima(), rand())
@btime JacobiElliptic.K(JacobiElliptic.Carlson(), rand())

using Enzyme
Enzyme.Compiler.RunAttributor[] = false
Enzyme.API.printall!(false)
φ = 0.11106641f0
m = 0.999f0

JacobiElliptic.CarlsonAlg.ellipke(m)#, m)

(k, φ, n) = (0.6276449453871871, 0.1435247713165507, 0.7888794230982612)
Enzyme.autodiff(Reverse, m -> JacobiElliptic.CarlsonAlg.F(φ, m), Active, Active(m))[1][1]
m = k^2
dk_dm = 0.5 / k
@benchmark autodiff(
    Reverse,
    m -> JacobiElliptic.CarlsonAlg.Pi(rand(), rand(), rand()),
    Active,
    Active(m),
)
@benchmark autodiff(
    Forward,
    m -> JacobiElliptic.CarlsonAlg.Pi(rand(), rand(), rand()),
    Duplicated,
    Duplicated(m, 1.0),
)
@benchmark autodiff(
    Reverse,
    m -> JacobiElliptic.FukushimaAlg.Pi(rand(), rand(), rand()),
    Active,
    Active(m),
)
@benchmark autodiff(
    Forward,
    m -> JacobiElliptic.FukushimaAlg.Pi(rand(), rand(), rand()),
    Duplicated,
    Duplicated(m, 1.0),
)

@benchmark autodiff(
    Reverse,
    m -> JacobiElliptic.CarlsonAlg.F(rand(), rand()),
    Active,
    Active(m),
)
@benchmark autodiff(
    Forward,
    m -> JacobiElliptic.CarlsonAlg.F(rand(), rand()),
    Duplicated,
    Duplicated(m, 1.0),
)
@benchmark autodiff(
    Reverse,
    m -> JacobiElliptic.FukushimaAlg.F(rand(), rand()),
    Active,
    Active(m),
)
@benchmark autodiff(
    Forward,
    m -> JacobiElliptic.FukushimaAlg.F(rand(), rand()),
    Duplicated,
    Duplicated(m, 1.0),
)


@btime autodiff(
    Reverse,
    m -> JacobiElliptic.CarlsonAlg.E(rand(), rand()),
    Active,
    Active(m),
)
@btime autodiff(
    Forward,
    m -> JacobiElliptic.CarlsonAlg.E(rand(), rand()),
    Duplicated,
    Duplicated(m, 1.0),
)
@time autodiff(Reverse, m -> JacobiElliptic.FukushimaAlg.E(φ, m), Active, Active(m))
@time autodiff(
    Forward,
    m -> JacobiElliptic.FukushimaAlg.E(φ, m),
    Duplicated,
    Duplicated(m, 1.0),
)



@time ForwardDiff.derivative(m -> Elliptic2.Pi(φ, m, n), m)# ≈ (E(m) - K(m))/k * dk_dm
@time ForwardDiff.derivative(E, m)# ≈ (E(m) - K(m))/k * dk_dm


test = x -> Elliptic2.SLATEC.DRF(0.0, 1.0 - x, 1.0)[1]
ForwardDiff.derivative(test, m)# ≈ (E(m) - K(m))/k * dk_dm


@btime Pi(rand(), rand(), rand())
@time ∂φF, ∂mF = autodiff(Reverse, F, Active, Active(φ), Active(m))[1]
#@test ∂φF ≈ typ(1)#/(1/sqrt(1 - m*sin(φ)^2)) ≈ typ(1) 
∂mF / (
    -(
        (
            E((φ), (m)) + (-1 + m) * F((φ), (m)) -
            (m * cos((φ)) * sin((φ))) / sqrt(1 - m * sin((φ))^2)
        ) / (2 * (-1 + m) * m)
    )
) ≈ 1.0
