using Metal
using CUDA

@testset "Elliptic K"  begin
    @testset "GPU" begin
        arr = [0.001f0, 0.01f0, 0.1f0, 0.2f0, 0.3f0, 0.5f0, 0.5f0, 0.8f0, 0.9f0, 0.99f0, 0.999f0]
        if Sys.isapple()
            @test length(Metal.devices()) > 0
            mtlarr = MtlArray(arr)
            gpuvals= Array(K.(mtlarr))

            @testset "$val" for (iter, val) in enumerate(arr)
                @test gpuvals[iter] ≈ elliptic_k(val)
            end
        else
            cuarr = CuArray(arr)
            gpuvals = Array(K.(cuarr))

            @testset "$val" for (iter, val) in enumerate(arr)
                @test gpuvals[iter] ≈ elliptic_k(val)
            end
        end
    end

end

@testset "Elliptic E"  begin
    @testset "GPU" begin
        arr = [0.001f0, 0.01f0, 0.1f0, 0.2f0, 0.3f0, 0.5f0, 0.5f0, 0.8f0, 0.9f0, 0.99f0, 0.999f0]
        if Sys.isapple()
            @test length(Metal.devices()) > 0
            mtlarr = MtlArray(arr)
            gpuvals= Array(E.(mtlarr))

            @testset "$val" for (iter, val) in enumerate(arr)
                @test gpuvals[iter] ≈ elliptic_e(val)
            end
        else
            cuarr = CuArray(arr)
            gpuvals= Array(E.(cuarr))

            @testset "$val" for (iter, val) in enumerate(arr)
                @test gpuvals[iter] ≈ elliptic_e(val)
            end
        end
    end
end

#Test EllipticF
@testset "Elliptic F" begin
    @testset "GPU" begin
        arr = [0.001f0, 0.01f0, 0.1f0, 0.2f0, 0.3f0, 0.5f0, 0.5f0, 0.8f0, 0.9f0, 0.99f0, 0.999f0]
        if Sys.isapple()
            @test length(Metal.devices()) > 0
            mtlarr = MtlArray(arr)
            @testset "φ: $φ" for φ in range(0f0, Float32(π/2), length=10)
                gpuvals= Array(F.(φ, mtlarr))

                @testset "$val" for (iter, val) in enumerate(arr)
                    @test gpuvals[iter] ≈ elliptic_f(φ, val)
                end
            end
        else
            cuarr = CuArray(arr)
            @testset "φ: $φ" for φ in range(0f0, Float32(π/2), length=10)
                gpuvals= Array(F.(φ, cuarr))

                @testset "$val" for (iter, val) in enumerate(arr)
                    @test gpuvals[iter] ≈ elliptic_f(φ, val)
                end
            end
        end
    end
end

#Test EllipticPI
@testset "Elliptic Pi" begin
    @testset "$typ" for typ in [Float32, Float64]
        @testset "Standard Domain" for φ in range(1e-3, π/2, length=100)
            @testset "φ :$(round(φ*180/pi, digits=2))" for n in range(1e-3,1-1e-3, length=100)
                @testset "m :$m)" for m in range(1e-3, 0.95,length=100)
                    @test FastElliptic.Pi(typ(n), typ(φ), typ(m))/ArbNumerics.elliptic_pi(typ(n), typ(φ), typ(m)) ≈ typ(1.0) atol=1e-4
                end
            end
        end
    end
end


#Taken from Elliptic.jl
#https://github.com/nolta/Elliptic.jl/blob/master/test/jacobi_tests.jl
@testset "Jacobi Elliptic Trig" begin
    @testset "GPU" begin
        arr = range(0.0f0, 3.14f0, length=10)
        if Sys.isapple()
            @test length(Metal.devices()) > 0

            mtlarr = MtlArray(arr)
            m = 0.1f0

            @test typeof(FastElliptic.sn.(mtlarr, m)) == MtlVector{Float32}
            @test typeof(FastElliptic.cn.(mtlarr, m)) == MtlVector{Float32}
            @test typeof(FastElliptic.dn.(mtlarr, m)) == MtlVector{Float32}
            @test typeof(FastElliptic.sd.(mtlarr, m)) == MtlVector{Float32}
            @test typeof(FastElliptic.sc.(mtlarr, m)) == MtlVector{Float32}
        else
            cuarr = CuArray(arr)
            m = 0.1f0

            @test typeof(FastElliptic.sn.(cuarr, m)) == MtlVector{Float32}
            @test typeof(FastElliptic.cn.(cuarr, m)) == MtlVector{Float32}
            @test typeof(FastElliptic.dn.(cuarr, m)) == MtlVector{Float32}
            @test typeof(FastElliptic.sd.(cuarr, m)) == MtlVector{Float32}
            @test typeof(FastElliptic.sc.(cuarr, m)) == MtlVector{Float32}

        end

    end
end