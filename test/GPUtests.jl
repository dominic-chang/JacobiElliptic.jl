using Metal
using CUDA

@testset "Elliptic K" begin
    @testset "GPU" begin
        arr = [
            0.001f0,
            0.01f0,
            0.1f0,
            0.2f0,
            0.3f0,
            0.5f0,
            0.5f0,
            0.8f0,
            0.9f0,
            0.99f0,
            0.999f0,
        ]
        mtlarr = MtlArray(arr)
        gpuvals = Array(K.(mtlarr))

        @testset "$val" for (iter, val) in enumerate(arr)
            @test gpuvals[iter] ≈ elliptic_k(val)
        end
    end

end

@testset "Elliptic E" begin
    @testset "GPU" begin
        arr = [
            0.001f0,
            0.01f0,
            0.1f0,
            0.2f0,
            0.3f0,
            0.5f0,
            0.5f0,
            0.8f0,
            0.9f0,
            0.99f0,
            0.999f0,
        ]
        mtlarr = MtlArray(arr)
        gpuvals = Array(E.(mtlarr))

        @testset "$val" for (iter, val) in enumerate(arr)
            @test gpuvals[iter] ≈ elliptic_e(val)
        end
    end
end

#Test EllipticF
@testset "Elliptic F" begin
    @testset "GPU" begin
        arr = [
            0.001f0,
            0.01f0,
            0.1f0,
            0.2f0,
            0.3f0,
            0.5f0,
            0.5f0,
            0.8f0,
            0.9f0,
            0.99f0,
            0.999f0,
        ]
        mtlarr = MtlArray(arr)
        @testset "φ: $φ" for φ in range(0.0f0, Float32(π / 2), length = 10)
            gpuvals = Array(F.(φ, mtlarr))

            @testset "$val" for (iter, val) in enumerate(arr)
                @test gpuvals[iter] ≈ elliptic_f(φ, val)
            end
        end
    end
end

#Test Complete EllipticPI
@testset "Complete Elliptic Pi" begin
    typ = Float32
    @testset "GPU" begin
        mtlarr = MtlArray([typ(n) for n in range(1e-3, 1e-3, length = 10)])
        @testset "Standard m" begin
            @testset "m :$m)" for m in range(zero(typ), one(typ), length = 10)
                JacobiElliptic.Pi.(mtlarr, typ(m))
            end
        end
    end
end

#Test EllipticPI
@testset "Elliptic Pi" begin
    typ = Float32
    @testset "GPU" begin
        mtlarr = MtlArray([φ for φ in range(1.0f-3, Float32(π / 2), length = 2)])
        @testset "n :$(round(φ*180/pi, digits=2))" for n in range(1e-3, 1e-3, length = 10)
            @testset "Standard m" begin
                @testset "m :$m)" for m in range(zero(typ), one(typ), length = 10)
                    gpuvals = JacobiElliptic.Pi.(typ(n), mtlarr, typ(m))
                end
            end
        end
    end
end


#Taken from Elliptic.jl
#https://github.com/nolta/Elliptic.jl/blob/main/test/jacobi_tests.jl
@testset "Jacobi Elliptic Trig" begin
    @testset "GPU" begin
        arr = range(0.0f0, 3.14f0, length = 10)
        if Sys.isapple()
            @test length(Metal.devices()) > 0

            mtlarr = MtlArray(arr)
            m = 0.1f0

            @test typeof(JacobiElliptic.sn.(mtlarr, m)) == MtlVector{Float32}
            @test typeof(JacobiElliptic.cn.(mtlarr, m)) == MtlVector{Float32}
            @test typeof(JacobiElliptic.dn.(mtlarr, m)) == MtlVector{Float32}
            @test typeof(JacobiElliptic.sd.(mtlarr, m)) == MtlVector{Float32}
            @test typeof(JacobiElliptic.sc.(mtlarr, m)) == MtlVector{Float32}
        else
            cuarr = CuArray(arr)
            m = 0.1f0

            @test typeof(JacobiElliptic.sn.(cuarr, m)) == MtlVector{Float32}
            @test typeof(JacobiElliptic.cn.(cuarr, m)) == MtlVector{Float32}
            @test typeof(JacobiElliptic.dn.(cuarr, m)) == MtlVector{Float32}
            @test typeof(JacobiElliptic.sd.(cuarr, m)) == MtlVector{Float32}
            @test typeof(JacobiElliptic.sc.(cuarr, m)) == MtlVector{Float32}

        end

    end
end
