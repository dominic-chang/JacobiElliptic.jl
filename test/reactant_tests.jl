@testset "Reactant" begin
    ms = Float64[0.01, 0.1, 0.5, 0.9]
    edge_ms = Float64[-0.5, 0.01, 1.5]
    pi_ns = Float64[0.0, 0.2, 0.8, 1.7]
    us = Float64[0.1, 0.4, 0.8, 1.2]
    short_us = Float64[0.1, 0.4, 0.8]
    wide_us = Float64[1.8, 2.4]

    m_r = Reactant.to_rarray(ms)
    edge_m_r = Reactant.to_rarray(edge_ms)
    n_r = Reactant.to_rarray(pi_ns)
    u_r = Reactant.to_rarray(us)
    short_u_r = Reactant.to_rarray(short_us)
    wide_u_r = Reactant.to_rarray(wide_us)

    scalar_u_r = Reactant.to_rarray(us[3:3])
    scalar_short_u_r = Reactant.to_rarray(short_us[2:2])
    scalar_n_r = Reactant.to_rarray(pi_ns[2:2])
    zero_n_r = Reactant.to_rarray(pi_ns[1:1])
    scalar_m_r = Reactant.to_rarray(ms[3:3])

    isapprox_with_nan =
        (actual, expected) -> all(
            pair -> (isnan(pair[1]) && isnan(pair[2])) || isapprox(pair[1], pair[2]),
            zip(Array(actual), expected),
        )

    @testset "Complete elliptic integrals" begin
        k_r = @jit K.(m_r)
        @test Array(k_r) ≈ K.(ms)

        k_edge_r = @jit K.(edge_m_r)
        @test Array(k_edge_r) ≈ K.(edge_ms)

        e_r = @jit E.(m_r)
        @test Array(e_r) ≈ E.(ms)

        e_edge_r = @jit E.(edge_m_r)
        @test Array(e_edge_r) ≈ E.(edge_ms)
    end

    @testset "Incomplete elliptic integrals" begin
        @testset "E" begin
            e_r = @jit E.(u_r, 0.5)
            @test Array(e_r) ≈ E.(us, 0.5)

            e_gt_one_r = @jit E.(short_u_r, 1.5)
            @test Array(e_gt_one_r) ≈ E.(short_us, 1.5)

            e_negative_r = @jit E.(u_r, -0.5)
            @test Array(e_negative_r) ≈ E.(us, -0.5)

            e_traced_m_r = @jit E.(u_r, m_r[3])
            @test Array(e_traced_m_r) ≈ E.(us, ms[3])

            e_traced_gt_one_m_r = @jit E.(short_u_r, edge_m_r[3])
            @test Array(e_traced_gt_one_m_r) ≈ E.(short_us, edge_ms[3])

            e_traced_negative_m_r = @jit E.(u_r, edge_m_r[1])
            @test Array(e_traced_negative_m_r) ≈ E.(us, edge_ms[1])

            e_wide_r = @jit E.(wide_u_r, 0.5)
            @test Array(e_wide_r) ≈ E.(wide_us, 0.5)

            e_wide_gt_one_r = @jit E.(wide_u_r, 1.5)
            @test isapprox_with_nan(e_wide_gt_one_r, E.(wide_us, 1.5))

            e_wide_negative_r = @jit E.(wide_u_r, -0.5)
            @test Array(e_wide_negative_r) ≈ E.(wide_us, -0.5)
        end

        @testset "F" begin
            f_r = @jit F.(u_r, 0.5)
            @test Array(f_r) ≈ F.(us, 0.5)

            f_gt_one_r = @jit F.(short_u_r, 1.5)
            @test Array(f_gt_one_r) ≈ F.(short_us, 1.5)

            f_negative_r = @jit F.(u_r, -0.5)
            @test Array(f_negative_r) ≈ F.(us, -0.5)

            f_traced_m_r = @jit F.(u_r, m_r[3])
            @test Array(f_traced_m_r) ≈ F.(us, ms[3])

            f_traced_gt_one_m_r = @jit F.(short_u_r, edge_m_r[3])
            @test Array(f_traced_gt_one_m_r) ≈ F.(short_us, edge_ms[3])

            f_traced_negative_m_r = @jit F.(u_r, edge_m_r[1])
            @test Array(f_traced_negative_m_r) ≈ F.(us, edge_ms[1])
        end

        @testset "Pi" begin
            pi_r = @jit Pi.(0.2, u_r, 0.5)
            @test Array(pi_r) ≈ Pi.(0.2, us, 0.5)

            pi_scalar_r = @jit Pi.(0.2, scalar_u_r, 0.5)
            @test Array(pi_scalar_r) ≈ Pi.(0.2, us[3:3], 0.5)

            pi_negative_m_r = @jit Pi.(0.2, u_r, edge_m_r[1])
            @test Array(pi_negative_m_r) ≈ Pi.(0.2, us, edge_ms[1])

            pi_n_gt_one_r = @jit Pi.(1.7, scalar_short_u_r, 0.5)
            @test Array(pi_n_gt_one_r) ≈ Pi.(1.7, short_us[2:2], 0.5)

            pi_wide_r = @jit Pi.(0.2, wide_u_r, 0.5)
            @test Array(pi_wide_r) ≈ Pi.(0.2, wide_us, 0.5)

            pi_traced_n_r = @jit Pi.(n_r, us[2], 0.5)
            @test Array(pi_traced_n_r) ≈ Pi.(pi_ns, us[2], 0.5)
        end
    end

    @testset "Complete elliptic integral Pi" begin
        pi_r = @jit Pi.(n_r, 0.5)
        @test Array(pi_r) ≈ Pi.(pi_ns, 0.5)

        pi_scalar_r = @jit Pi.(scalar_n_r, scalar_m_r)
        @test Array(pi_scalar_r) ≈ Pi.(pi_ns[2:2], ms[3:3])

        pi_zero_n_r = @jit Pi.(zero_n_r, scalar_m_r)
        @test Array(pi_zero_n_r) ≈ Pi.(pi_ns[1:1], ms[3:3])

        pi_negative_m_r = @jit Pi.(n_r, edge_m_r[1])
        @test Array(pi_negative_m_r) ≈ Pi.(pi_ns, edge_ms[1])

        pi_gt_one_m_r = @jit Pi.(n_r, edge_m_r[3])
        @test isapprox_with_nan(pi_gt_one_m_r, Pi.(pi_ns, edge_ms[3]))
    end

    @testset "Jacobi amplitude" begin
        am_r = @jit am.(u_r, 0.5)
        @test Array(am_r) ≈ am.(us, 0.5)

        am_gt_one_r = @jit am.(short_u_r, 1.5)
        @test Array(am_gt_one_r) ≈ am.(short_us, 1.5)

        am_negative_r = @jit am.(u_r, -0.5)
        @test Array(am_negative_r) ≈ am.(us, -0.5)

        am_traced_m_r = @jit am.(u_r, m_r[3])
        @test Array(am_traced_m_r) ≈ am.(us, ms[3])

        am_traced_gt_one_m_r = @jit am.(short_u_r, edge_m_r[3])
        @test Array(am_traced_gt_one_m_r) ≈ am.(short_us, edge_ms[3])

        am_traced_negative_m_r = @jit am.(u_r, edge_m_r[1])
        @test Array(am_traced_negative_m_r) ≈ am.(us, edge_ms[1])
    end

    @testset "Jacobi elliptic functions" begin
        @testset "Primary functions" begin
            sn_r = @jit sn.(u_r, 0.5)
            @test Array(sn_r) ≈ sn.(us, 0.5)

            cn_r = @jit cn.(u_r, 0.5)
            @test Array(cn_r) ≈ cn.(us, 0.5)

            dn_r = @jit dn.(u_r, 0.5)
            @test Array(dn_r) ≈ dn.(us, 0.5)

            sn_negative_r = @jit sn.(u_r, -0.5)
            @test Array(sn_negative_r) ≈ sn.(us, -0.5)

            cn_negative_r = @jit cn.(u_r, -0.5)
            @test Array(cn_negative_r) ≈ cn.(us, -0.5)

            dn_negative_r = @jit dn.(u_r, -0.5)
            @test Array(dn_negative_r) ≈ dn.(us, -0.5)

            sn_gt_one_r = @jit sn.(short_u_r, 1.5)
            @test Array(sn_gt_one_r) ≈ sn.(short_us, 1.5)

            cn_gt_one_r = @jit cn.(short_u_r, 1.5)
            @test Array(cn_gt_one_r) ≈ cn.(short_us, 1.5)

            dn_gt_one_r = @jit dn.(short_u_r, 1.5)
            @test Array(dn_gt_one_r) ≈ dn.(short_us, 1.5)
        end

        @testset "Quotient functions" begin
            sd_r = @jit sd.(u_r, 0.5)
            @test Array(sd_r) ≈ sd.(us, 0.5)

            cd_r = @jit JacobiElliptic.cd.(u_r, 0.5)
            @test Array(cd_r) ≈ JacobiElliptic.cd.(us, 0.5)

            nd_r = @jit nd.(u_r, 0.5)
            @test Array(nd_r) ≈ nd.(us, 0.5)

            sc_r = @jit sc.(u_r, 0.5)
            @test Array(sc_r) ≈ sc.(us, 0.5)

            dc_r = @jit dc.(u_r, 0.5)
            @test Array(dc_r) ≈ dc.(us, 0.5)

            nc_r = @jit nc.(u_r, 0.5)
            @test Array(nc_r) ≈ nc.(us, 0.5)

            cs_r = @jit cs.(u_r, 0.5)
            @test Array(cs_r) ≈ cs.(us, 0.5)

            ds_r = @jit ds.(u_r, 0.5)
            @test Array(ds_r) ≈ ds.(us, 0.5)

            ns_r = @jit ns.(u_r, 0.5)
            @test Array(ns_r) ≈ ns.(us, 0.5)
        end
    end

    # TODO: Add Reactant coverage for the constant quotient functions
    # (nn, dd, cc, and ss) and explicit algorithm dispatch.
end
