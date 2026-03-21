@testset "Metal" begin
    metal_loaded = false

    @static if Sys.isapple()
        try
            @eval using Metal
            metal_loaded = true
        catch _
            metal_loaded = false
        end
    end

    if !Sys.isapple()
        @info "Skipping Metal tests on non-Apple platforms."
        @test true
    elseif !metal_loaded
        @info "Skipping Metal tests because Metal.jl is not available in the test environment."
        @test true
    elseif isempty(Metal.devices())
        @info "Skipping Metal tests because no Metal device is available."
        @test true
    else
        ms = Float32[0.01, 0.1, 0.5, 0.9]
        us = Float32[0.1, 0.4, 0.8, 1.2]
        m_gpu = Metal.MtlArray(ms)
        u_gpu = Metal.MtlArray(us)

        @testset "Default API" begin
            k_gpu = K.(m_gpu)
            @test k_gpu isa Metal.MtlArray
            e_gpu = E.(m_gpu)
            @test e_gpu isa Metal.MtlArray
            f_gpu = F.(u_gpu, 0.5f0)
            @test f_gpu isa Metal.MtlArray
            pi_gpu = Pi.(0.2f0, u_gpu, 0.5f0)
            @test pi_gpu isa Metal.MtlArray

            @test Array(k_gpu) ≈ K.(ms)
            @test Array(e_gpu) ≈ E.(ms)
            @test Array(f_gpu) ≈ F.(us, 0.5f0)
            @test Array(pi_gpu) ≈ Pi.(0.2f0, us, 0.5f0)
        end

        @testset "Jacobi API" begin
            sn_gpu = sn.(u_gpu, 0.5f0)
            @test sn_gpu isa Metal.MtlArray
            cn_gpu = cn.(u_gpu, 0.5f0)
            @test cn_gpu isa Metal.MtlArray
            dn_gpu = dn.(u_gpu, 0.5f0)
            @test dn_gpu isa Metal.MtlArray
            sc_gpu = sc.(u_gpu, 0.5f0)
            @test sc_gpu isa Metal.MtlArray
            sd_gpu = sd.(u_gpu, 0.5f0)
            @test sd_gpu isa Metal.MtlArray

            @test Array(sn_gpu) ≈ sn.(us, 0.5f0)
            @test Array(cn_gpu) ≈ cn.(us, 0.5f0)
            @test Array(dn_gpu) ≈ dn.(us, 0.5f0)
            @test Array(sc_gpu) ≈ sc.(us, 0.5f0)
            @test Array(sd_gpu) ≈ sd.(us, 0.5f0)
        end

        @testset "Algorithm Dispatch" begin
            fukushima_k_gpu = K.(Ref(Fukushima()), m_gpu)
            @test fukushima_k_gpu isa Metal.MtlArray
            carlson_f_gpu = F.(Ref(Carlson()), u_gpu, 0.5f0)
            @test carlson_f_gpu isa Metal.MtlArray
            fukushima_sn_gpu = sn.(Ref(Fukushima()), u_gpu, 0.5f0)
            @test fukushima_sn_gpu isa Metal.MtlArray

            @test Array(fukushima_k_gpu) ≈ K.(Ref(Fukushima()), ms)
            @test Array(carlson_f_gpu) ≈ F.(Ref(Carlson()), us, 0.5f0)
            @test Array(fukushima_sn_gpu) ≈ sn.(Ref(Fukushima()), us, 0.5f0)
        end
    end
end
