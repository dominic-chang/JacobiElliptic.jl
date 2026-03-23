
@testset "Reactant" begin

    reactant_jit(expr::Expr) = Core.eval(@__MODULE__, :(Reactant.@jit $expr))

    ms = Float64[0.01, 0.1, 0.5, 0.9]
    k_edge_ms = Float64[-0.5, 0.01, 1.5]
    us = Float64[0.1, 0.4, 0.8, 1.2]
    gt_one_us = Float64[0.1, 0.4, 0.8]
    m_r = Reactant.to_rarray(ms)
    k_edge_m_r = Reactant.to_rarray(k_edge_ms)
    u_r = Reactant.to_rarray(us)
    gt_one_u_r = Reactant.to_rarray(gt_one_us)

    @testset "Default API" begin
        k_r = @jit K.(m_r)
        k_edge_r = @jit K.(k_edge_m_r)
        e_r = @jit E.(m_r)
        e_edge_r = @jit E.(k_edge_m_r)
        f_r = reactant_jit(:(F.($u_r, 0.5)))
        f_gt_one_r = reactant_jit(:(F.($gt_one_u_r, 1.5)))
        f_neg_r = reactant_jit(:(F.($u_r, -0.5)))
        f_m_r = reactant_jit(:(F.($u_r, $m_r[3])))
        f_m_gt_one_r = reactant_jit(:(F.($gt_one_u_r, $k_edge_m_r[3])))
        f_m_neg_r = reactant_jit(:(F.($u_r, $k_edge_m_r[1])))
        am_r = reactant_jit(:(am.($u_r, 0.5)))
        am_gt_one_r = reactant_jit(:(am.($gt_one_u_r, 1.5)))
        am_neg_r = reactant_jit(:(am.($u_r, -0.5)))
        am_m_r = reactant_jit(:(am.($u_r, $m_r[3])))
        am_m_gt_one_r = reactant_jit(:(am.($gt_one_u_r, $k_edge_m_r[3])))
        am_m_neg_r = reactant_jit(:(am.($u_r, $k_edge_m_r[1])))
        #pi_r = reactant_jit(:(Pi.(0.2, $u_r, 0.5)))
        #pi_scalar_r = reactant_jit(:(Pi(0.2, $u_r[3], 0.5)))
        #pi_neg_m_r = reactant_jit(:(Pi(0.2, $u_r, -0.5)))
        #pi_gt_one_m_r = reactant_jit(:(Pi(0.2, $u_r, 1.5)))

        @test Array(k_r) ≈ K.(ms)
        @test Array(k_edge_r) ≈ K.(k_edge_ms)
        @test Array(e_r) ≈ E.(ms)
        @test Array(e_edge_r) ≈ E.(k_edge_ms)
        @test Array(f_r) ≈ F.(us, 0.5)
        @test Array(f_gt_one_r) ≈ F.(gt_one_us, 1.5)
        @test Array(f_neg_r) ≈ F.(us, -0.5)
        @test Array(f_m_r) ≈ F.(us, ms[3])
        @test Array(f_m_gt_one_r) ≈ F.(gt_one_us, k_edge_ms[3])
        @test Array(f_m_neg_r) ≈ F.(us, k_edge_ms[1])
        @test Array(am_r) ≈ am.(us, 0.5)
        @test Array(am_gt_one_r) ≈ am.(gt_one_us, 1.5)
        @test Array(am_neg_r) ≈ am.(us, -0.5)
        @test Array(am_m_r) ≈ am.(us, ms[3])
        @test Array(am_m_gt_one_r) ≈ am.(gt_one_us, k_edge_ms[3])
        @test Array(am_m_neg_r) ≈ am.(us, k_edge_ms[1])
        #@test Array(pi_r) ≈ Pi.(0.2, us, 0.5)
        #@test pi_scalar_r ≈ Pi(0.2, us[3], 0.5)
        #@test Array(pi_neg_m_r) ≈ Pi.(0.2, us, -0.5)
        #pi_gt_one_expected = Pi.(0.2, us, 1.5)
        #@test all(
        #    ((isnan(a) && isnan(b)) || isapprox(a, b)) for
        #    (a, b) in zip(Array(pi_gt_one_m_r), pi_gt_one_expected)
        #)
    end

    #@testset "Jacobi API" begin
    #    sn_r = reactant_jit(:(sn.($u_r, 0.5)))
    #    cn_r = reactant_jit(:(cn.($u_r, 0.5)))
    #    dn_r = reactant_jit(:(dn.($u_r, 0.5)))
    #    sc_r = reactant_jit(:(sc.($u_r, 0.5)))
    #    sd_r = reactant_jit(:(sd.($u_r, 0.5)))

    #    @test Array(sn_r) ≈ sn.(us, 0.5)
    #    @test Array(cn_r) ≈ cn.(us, 0.5)
    #    @test Array(dn_r) ≈ dn.(us, 0.5)
    #    @test Array(sc_r) ≈ sc.(us, 0.5)
    #    @test Array(sd_r) ≈ sd.(us, 0.5)
    #end

    #@testset "Algorithm Dispatch" begin
    #    #fukushima_k_r = reactant_jit(:(K.(Ref(Fukushima()), $m_r)))
    #    carlson_f_r = reactant_jit(:(F.(Ref(Carlson()), $u_r, 0.5)))
    #    carlson_f_gt_one_r = reactant_jit(:(F.(Ref(Carlson()), $u_r, 1.5)))
    #    carlson_f_neg_r = reactant_jit(:(F.(Ref(Carlson()), $u_r, -0.5)))
    #    carlson_pi_r = reactant_jit(:(Pi.(Ref(Carlson()), 0.2, $u_r, 0.5)))
    #    #fukushima_sn_r = reactant_jit(:(sn.(Ref(Fukushima()), $u_r, 0.5)))

    #    #@test Array(fukushima_k_r) ≈ K.(Ref(Fukushima()), ms)
    #    @test Array(carlson_f_r) ≈ F.(Ref(Carlson()), us, 0.5)
    #    @test Array(carlson_f_gt_one_r) ≈ F.(Ref(Carlson()), us, 1.5)
    #    @test Array(carlson_f_neg_r) ≈ F.(Ref(Carlson()), us, -0.5)
    #    @test Array(carlson_pi_r) ≈ Pi.(Ref(Carlson()), 0.2, us, 0.5)
    #    #@test Array(fukushima_sn_r) ≈ sn.(Ref(Fukushima()), us, 0.5)
    #end
end
