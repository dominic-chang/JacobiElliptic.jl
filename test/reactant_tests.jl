
@testset "Reactant" begin

    reactant_jit(expr::Expr) = Core.eval(@__MODULE__, :(Reactant.@jit $expr))
    reactant_ext = Base.get_extension(JacobiElliptic, :JacobiEllipticReactantExt)

    ms = Float64[0.01, 0.1, 0.5, 0.9]
    pi_ns = Float64[0.0, 0.2, 0.8, 1.7]
    k_edge_ms = Float64[-0.5, 0.01, 1.5]
    us = Float64[0.1, 0.4, 0.8, 1.2]
    gt_one_us = Float64[0.1, 0.4, 0.8]
    wide_us = Float64[1.8, 2.4]
    m_r = Reactant.to_rarray(ms)
    n_r = Reactant.to_rarray(pi_ns)
    k_edge_m_r = Reactant.to_rarray(k_edge_ms)
    u_r = Reactant.to_rarray(us)
    gt_one_u_r = Reactant.to_rarray(gt_one_us)
    wide_u_r = Reactant.to_rarray(wide_us)
    @testset "Default API" begin
        k_r = @jit K.(m_r)
        k_edge_r = @jit K.(k_edge_m_r)
        e_r = @jit E.(m_r)
        e_edge_r = @jit E.(k_edge_m_r)
        incomplete_e_r = [reactant_jit(:(E($(u), 0.5))) for u in us]
        incomplete_e_gt_one_r = [reactant_jit(:(E($(u), 1.5))) for u in gt_one_us]
        incomplete_e_neg_r = [reactant_jit(:(E($(u), -0.5))) for u in us]
        incomplete_e_m_r = [reactant_jit(:(E($(u), $(ms[3])))) for u in us]
        incomplete_e_m_gt_one_r = [reactant_jit(:(E($(u), $(k_edge_ms[3])))) for u in gt_one_us]
        incomplete_e_m_neg_r = [reactant_jit(:(E($(u), $(k_edge_ms[1])))) for u in us]
        incomplete_e_r = reactant_jit(:(E.($(u_r), 0.5)))
        incomplete_e_gt_one_r = reactant_jit(:(E.($(gt_one_u_r), 1.5)))
        incomplete_e_neg_r = reactant_jit(:(E.($(u_r), -0.5)))
        incomplete_e_m_r = reactant_jit(:(E.($(u_r), $(m_r[3]))))
        incomplete_e_m_gt_one_r = reactant_jit(:(E.($(gt_one_u_r), $(k_edge_m_r[3]))))
        incomplete_e_m_neg_r = reactant_jit(:(E.($(u_r), $(k_edge_m_r[1]))))
        incomplete_e_wide_r = reactant_jit(:(E.($(wide_u_r), 0.5)))
        incomplete_e_wide_gt_one_r = reactant_jit(:(E.($(wide_u_r), 1.5)))
        incomplete_e_wide_neg_r = reactant_jit(:(E.($(wide_u_r), -0.5)))
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
        complete_pi_r = reactant_jit(:(Pi.($n_r, 0.5)))
        complete_pi_scalar_r = reactant_jit(:(Pi($(pi_ns[2]), $(ms[3]))))
        complete_pi_zero_n_r = reactant_jit(:(Pi($(pi_ns[1]), $(ms[3]))))
        complete_pi_neg_m_r = reactant_jit(:(Pi.($n_r, $(k_edge_m_r[1]))))
        complete_pi_gt_one_m_r = reactant_jit(:(Pi.($n_r, $(k_edge_m_r[3]))))

        @test Array(k_r) ≈ K.(ms)
        @test Array(k_edge_r) ≈ K.(k_edge_ms)
        @test Array(e_r) ≈ E.(ms)
        @test Array(e_edge_r) ≈ E.(k_edge_ms)
        @test incomplete_e_r ≈ E.(us, 0.5)
        @test incomplete_e_gt_one_r ≈ E.(gt_one_us, 1.5)
        @test incomplete_e_neg_r ≈ E.(us, -0.5)
        @test incomplete_e_m_r ≈ E.(us, ms[3])
        @test incomplete_e_m_gt_one_r ≈ E.(gt_one_us, k_edge_ms[3])
        @test incomplete_e_m_neg_r ≈ E.(us, k_edge_ms[1])
        @test incomplete_e_wide_r ≈ E.(wide_us, 0.5)
        wide_gt_one_expected = E.(wide_us, 1.5)
        @test all(
            ((isnan(a) && isnan(b)) || isapprox(a, b)) for
            (a, b) in zip(incomplete_e_wide_gt_one_r, wide_gt_one_expected)
        )
        @test incomplete_e_wide_neg_r ≈ E.(wide_us, -0.5)
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
        @test Array(complete_pi_r) ≈ Pi.(pi_ns, 0.5)
        @test complete_pi_scalar_r ≈ Pi(pi_ns[2], ms[3])
        @test complete_pi_zero_n_r ≈ Pi(pi_ns[1], ms[3])
        @test Array(complete_pi_neg_m_r) ≈ Pi.(pi_ns, k_edge_ms[1])
        complete_pi_gt_one_expected = Pi.(pi_ns, k_edge_ms[3])
        @test all(
            ((isnan(a) && isnan(b)) || isapprox(a, b)) for
            (a, b) in zip(Array(complete_pi_gt_one_m_r), complete_pi_gt_one_expected)
        )
    end

    @testset "Jacobi API" begin
        sn_r = reactant_jit(:(sn.($u_r, 0.5)))
        cn_r = reactant_jit(:(cn.($u_r, 0.5)))
        dn_r = reactant_jit(:(dn.($u_r, 0.5)))
        sn_neg_r = reactant_jit(:(sn.($u_r, -0.5)))
        cn_neg_r = reactant_jit(:(cn.($u_r, -0.5)))
        dn_neg_r = reactant_jit(:(dn.($u_r, -0.5)))
        sn_gt_one_r = reactant_jit(:(sn.($gt_one_u_r, 1.5)))
        cn_gt_one_r = reactant_jit(:(cn.($gt_one_u_r, 1.5)))
        dn_gt_one_r = reactant_jit(:(dn.($gt_one_u_r, 1.5)))
        sd_r = reactant_jit(:(sd.($u_r, 0.5)))
        cd_r = reactant_jit(:(JacobiElliptic.cd.($u_r, 0.5)))
        #dd_r = reactant_jit(:(dd.($u_r, 0.5)))
        nd_r = reactant_jit(:(nd.($u_r, 0.5)))
        sc_r = reactant_jit(:(sc.($u_r, 0.5)))
        #cc_r = reactant_jit(:(cc.($u_r, 0.5)))
        dc_r = reactant_jit(:(dc.($u_r, 0.5)))
        nc_r = reactant_jit(:(nc.($u_r, 0.5)))
        #ss_r = reactant_jit(:(ss.($u_r, 0.5)))
        cs_r = reactant_jit(:(cs.($u_r, 0.5)))
        ds_r = reactant_jit(:(ds.($u_r, 0.5)))
        ns_r = reactant_jit(:(ns.($u_r, 0.5)))
        #nn_r = reactant_jit(:(nn.($u_r, 0.5)))

        @test Array(sn_r) ≈ sn.(us, 0.5)
        @test Array(cn_r) ≈ cn.(us, 0.5)
        @test Array(dn_r) ≈ dn.(us, 0.5)
        @test Array(sn_neg_r) ≈ sn.(us, -0.5)
        @test Array(cn_neg_r) ≈ cn.(us, -0.5)
        @test Array(dn_neg_r) ≈ dn.(us, -0.5)
        @test Array(sn_gt_one_r) ≈ sn.(gt_one_us, 1.5)
        @test Array(cn_gt_one_r) ≈ cn.(gt_one_us, 1.5)
        @test Array(dn_gt_one_r) ≈ dn.(gt_one_us, 1.5)
        @test Array(sd_r) ≈ sd.(us, 0.5)
        @test Array(cd_r) ≈ JacobiElliptic.cd.(us, 0.5)
        #@test Array(dd_r) ≈ dd.(us, 0.5)
        @test Array(nd_r) ≈ nd.(us, 0.5)
        @test Array(sc_r) ≈ sc.(us, 0.5)
        #@test Array(cc_r) ≈ cc.(us, 0.5)
        @test Array(dc_r) ≈ dc.(us, 0.5)
        @test Array(nc_r) ≈ nc.(us, 0.5)
        #@test Array(ss_r) ≈ ss.(us, 0.5)
        @test Array(cs_r) ≈ cs.(us, 0.5)
        @test Array(ds_r) ≈ ds.(us, 0.5)
        @test Array(ns_r) ≈ ns.(us, 0.5)
        #@test Array(nn_r) ≈ nn.(us, 0.5)
    end

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
