
# Testing Algorithm uality by comparing Carlson Elliptic Integral Algorithm implementation in Arbnumerics
function ArbNumerics.elliptic_k(m::T) where T 
    return ArbNumerics.elliptic_k(ArbNumerics.ArbFloat(m))
end
function ArbNumerics.elliptic_k(m::Float64) 
    return ArbNumerics.elliptic_k(ArbNumerics.ArbFloat(m))
end
function ArbNumerics.elliptic_e(m::T) where T 
    return ArbNumerics.elliptic_e(ArbNumerics.ArbFloat(m))
end
function ArbNumerics.elliptic_e(m::Float64) 
    return ArbNumerics.elliptic_e(ArbNumerics.ArbFloat(m))
end
function ArbNumerics.elliptic_f(φ::T, m::T) where T 
    return ArbNumerics.elliptic_f(ArbNumerics.ArbFloat(φ), ArbNumerics.ArbFloat(m))
end
function ArbNumerics.elliptic_f(φ::Float64, m::Float64)
    return ArbNumerics.elliptic_f(ArbNumerics.ArbFloat(φ), ArbNumerics.ArbFloat(m))
end
function ArbNumerics.elliptic_e(φ::T, m::T) where T 
    return ArbNumerics.elliptic_e(ArbNumerics.ArbFloat(φ), ArbNumerics.ArbFloat(m))
end
function ArbNumerics.elliptic_e(φ::Float64, m::Float64)
    return ArbNumerics.elliptic_e(ArbNumerics.ArbFloat(φ), ArbNumerics.ArbFloat(m))
end
ArbNumerics.elliptic_pi(n::AbstractFloat, φ::AbstractFloat, m::AbstractFloat) = ArbNumerics.elliptic_pi(ArbNumerics.ArbFloat(n), ArbNumerics.ArbFloat(φ), ArbNumerics.ArbFloat(m))
function ArbNumerics.elliptic_pi(n::Float64, m::Float64)
    return ArbNumerics.elliptic_pi(ArbNumerics.ArbFloat(n), ArbNumerics.ArbFloat(m))
end
function ArbNumerics.elliptic_pi(n::Float32, m::Float32)
    return ArbNumerics.elliptic_pi(ArbNumerics.ArbFloat(n), ArbNumerics.ArbFloat(m))
end

@testset "alg:$alg" for alg in [JacobiElliptic.CarlsonAlg, JacobiElliptic.FukushimaAlg]
    @testset "Elliptic K"  begin
        @testset for typ in [Float32, Float64]
            @testset "$typ $m" for m in [0.001, 0.01, 0.1, 0.2, 0.3, 0.5, 0.5, 0.8, 0.9, 0.99, 0.999]
                @test alg.K(typ(m)) ≈ typ(ArbNumerics.elliptic_k(typ(m)))
            end
        end
    end

    @testset "Complete Elliptic E"  begin
        @testset "Type: $typ" for typ in [Float32, Float64]
            @testset  "$(typ(m))" for m in [0.001, 0.01, 0.1, 0.19, 0.29, 0.39, 0.49, 0.59, 0.69, 0.79, 0.89, 0.99, 0.999]
                @test alg.E(typ(m)) ≈ ArbNumerics.elliptic_e(typ(m)) 
            end
        end
    end

    #Test EllipticF
    @testset "Elliptic F" begin
        @testset "$typ" for typ in [Float32, Float64]
            @testset "Standard Domain, φ:$(typ(φ))" for φ in range(-π/2, π/2, length=100)
                @testset "m=$(typ(m))" for m in [0.001, 0.01, 0.1, 0.2, 0.3, 0.5, 0.5, 0.8, 0.9, 0.99, 0.999]
                    @test alg.F(typ(φ), typ(m)) / typ(ArbNumerics.elliptic_f(typ(φ), typ(m))) ≈ one(typ) atol=1e-4
                end
            end
            @testset "abs(φ) > π/2:" begin
                @testset "φ = $(typ(φ))" for φ in range(-4π/2, 4π/2, length=100)
                    @testset "m=$(typ(m))" for m in range(1e-3, 1-1e-3,length=10)
                        @test alg.F(typ(φ), typ(m)) / ArbNumerics.elliptic_f(typ(φ), typ(m)) ≈ one(typ) atol=1e-4
                    end
                end
            end
            @testset "|m sin(φ)| ≤ 1 && |m| > 1" begin
                @testset "m = $(typ(m))" for m in [-100, -10, -5, -2, -1.1, 1.1, 2, 5, 10, 100]
                    mmax = asin(1/√abs(m))
                    @testset "φ = $(i*mmax)" for i in -1:0.1:1.0
                        φ = i*(mmax*(1-1e-3))
                        if φ != zero(typ)
                            @test alg.F(typ(φ), typ(m)) / ArbNumerics.elliptic_f(typ(φ), typ(m)) ≈ one(typ) atol=1e-4
                        end
                    end
                end
            end
        end
    end

    #Test EllipticE
    @testset "Elliptic E" begin
        @testset "$typ" for typ in [Float32, Float64]
            @testset "Standard Domain, φ:$(typ(φ))" for φ in range(-π/2, π/2, length=100)
                @testset "m=$(typ(m))" for m in [0.001, 0.01, 0.1, 0.2, 0.3, 0.5, 0.5, 0.8, 0.9, 0.99, 0.999]
                    @test alg.E(typ(φ), typ(m)) / typ(ArbNumerics.elliptic_e(typ(φ), typ(m))) ≈ one(typ) atol=1e-4
                end
            end
            @testset "abs(φ) > π/2:" begin
                @testset "φ = $(typ(φ))" for φ in range(-4π/2, 4π/2, length=100)
                    @testset "m=$(typ(m))" for m in range(1e-3, 1-1e-3,length=10)
                        @test alg.E(typ(φ), typ(m)) / ArbNumerics.elliptic_e(typ(φ), typ(m)) ≈ one(typ) atol=1e-4
                    end
                end
            end
            @testset "|m sin(φ)| ≤ 1 && |m| > 1" begin
                @testset "m = $(typ(m))" for m in [1.1, 2, 5, 10, 100]
                    mmax = asin(1/√abs(m))
                    @testset "φ = $(i*mmax)" for i in -1:0.1:1.0
                        φ = i*(mmax*(0.999))
                        if φ != zero(typ)
                            @test alg.E(typ(φ), typ(m)) / ArbNumerics.elliptic_e(typ(φ), typ(m)) ≈ one(typ) atol=1e-4
                        end
                    end
                end
            end
        end
    end
    #Test EllipticPI
    @testset "Elliptic Pi" begin
        @testset "$typ" for typ in [Float32, Float64]
            @testset "Standard Domain" begin
                @testset "φ : $φ" for φ in range(1e-3, π/2, length=100)
                    @testset "n : $n" for n in range(1e-3,typ == Float32 ? 0.99 : 0.999, length=10)
                        @testset "Standard m" begin
                            @testset "m : $m" for m in range(0.1, typ == Float32 ? 0.99 : 0.999,length=20)
                                @test alg.Pi(typ(n), typ(φ), typ(m))/ArbNumerics.elliptic_pi(typ(n), typ(φ), typ(m)) ≈ typ(1.0) atol=1e-4
                            end
                        end
                        @testset "small m" begin
                            @testset "m :$m" for m in [3.874e-4, 4.25e-3, 6.83e-2]
                                @test alg.Pi(typ(n), typ(φ), typ(m))/ArbNumerics.elliptic_pi(typ(n), typ(φ), typ(m)) ≈ typ(1.0) atol=1e-4
                            end
                        end
                        @testset "large m" begin
                            @testset "m :$m" for m in [1-3.874e-4, 1-4.25e-3, 1-6.83e-2]
                                @test alg.Pi(typ(n), typ(φ), typ(m))/ArbNumerics.elliptic_pi(typ(n), typ(φ), typ(m)) ≈ typ(1.0) atol=1e-4
                            end
                        end
                        @testset "negative m" begin
                            @testset "m :$m" for m in [-3.874e-4, -4.25e-3, -6.83e-2]
                                @test alg.Pi(typ(n), typ(φ), typ(m))/ArbNumerics.elliptic_pi(typ(n), typ(φ), typ(m)) ≈ typ(1.0) atol=1e-4
                            end
                        end
                    end
                end
            end
            @testset "|m sin(φ)| ≤ 1 && |m| > 1" begin
                @testset "m = $(typ(m))" for m in [-100, -10, -5, -2, -1.1, 1.1, 2, 5, 10, 100]
                    @testset "n : $n" for n in range(1e-3, typ == Float32 ? 0.99 : 0.999, length=10)
                        mmax = asin(1/√abs(m))
                        @testset "φ = $(i*mmax)" for i in -1:0.1:1.0
                            φ = i*(mmax*(1-1e-3))
                            if φ != zero(typ)
                                @test alg.Pi(typ(n), typ(φ), typ(m)) / ArbNumerics.elliptic_pi(typ(n), typ(φ), typ(m)) ≈ one(typ) atol=1e-4
                            end
                        end
                    end
                end
            end
            @testset "φ > π/2" begin
                @testset "φ : $φ" for φ in range(π/2, 3π, length=10)
                    @testset "n : $n" for n in range(1e-3,typ == Float32 ? 0.99 : 0.999, length=10)
                        @testset "m : $m" for m in range(1e-3, typ == Float32 ? 0.99 : 0.999,length=20)
                            @test alg.Pi(typ(n), typ(φ), typ(m))/ArbNumerics.elliptic_pi(typ(n), typ(φ), typ(m)) ≈ typ(1.0) atol=1e-5
                        end
                        @testset "large m" begin
                            @testset "m :$m)" for m in [1-1e-5, 1-3.874e-4, 1-4.25e-3, 1-6.83e-2, 1-1.32e-1]
                                @test alg.Pi(typ(n), typ(φ), typ(m))/ArbNumerics.elliptic_pi(typ(n), typ(φ), typ(m)) ≈ typ(1.0) atol=1e-5
                            end
                        end
                    end
                end
            end
        end
    end
    
    #Test EllipticPI
    @testset "Complete Elliptic Pi" begin
        @testset "$typ" for typ in [Float32, Float64]
            @testset "n : $n" for n in range(1e-3,10, length=10)
                @testset "m : $m" for m in range(0.1, 0.999,length=20)
                    @test alg.Pi(typ(n), typ(m))/ArbNumerics.elliptic_pi(typ(n), typ(m)) ≈ typ(1.0) atol=1e-4
                end
            end
        end
    end
    
    #Taken from Elliptic.jl
    #https://github.com/nolta/Elliptic.jl/blob/main/test/jacobi_tests.jl
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
                
                    @test alg.sn(u, m) ≈ θs / θn atol=2.5e-6
                    @test alg.cn(u, m) ≈ θc / θn atol=2.5e-6
                    @test alg.dn(u, m) ≈ θd / θn atol=2.5e-6
                
                    @test alg.sd(u, m) ≈ θs / θd atol=2.5e-6
                    if ϵ != 90
                        # very sensitive around u = K,
                        # estimate of K(0) = pi/2 + 4e-9, so cosine causes errors
                        # also, errors build up in ϵ ≥75°, so lower tolerences for that
                        @test alg.sc(u, m) ≈ θs / θc atol=1e-7
                    end
                end
            end
        end
       @testset "u = $u" for u in -1.:0.21:1.0
           @test alg.sn(u,0.) ≈ sin(u)
           @test alg.cn(u,0.) ≈ cos(u)
           @test alg.dn(u,0.) ≈ 1.
           @test alg.sd(u,0.) ≈ sin(u)
           @test alg.sc(u,0.) ≈ tan(u)
    
           @test alg.sn(u,1.) ≈ tanh(u)
           @test alg.cn(u,1.) ≈ sech(u)
           @test alg.dn(u,1.) ≈ sech(u)
           @test alg.sd(u,1.) ≈ sinh(u)
           @test alg.sc(u,1.) ≈ sinh(u)
       end
   
    end

end