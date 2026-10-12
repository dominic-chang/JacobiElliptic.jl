@testset "Default algorithm selection" begin
    original = default_algorithm()

    try
        @test original isa AGM

        @test set_default_algorithm!(Fukushima()) isa Fukushima
        @test default_algorithm() isa Fukushima
        @test K(0.4) == K(Fukushima(), 0.4)
        @test E(0.7, 0.4) == E(Fukushima(), 0.7, 0.4)
        @test F(0.7, 0.4) == F(Fukushima(), 0.7, 0.4)
        @test Pi(0.2, 0.7, 0.4) == Pi(Fukushima(), 0.2, 0.7, 0.4)
        @test J(0.2, 0.7, 0.4) == J(Fukushima(), 0.2, 0.7, 0.4)
        @test sn(0.7, 0.4) == sn(Fukushima(), 0.7, 0.4)
        @test am(0.7, 0.4) == am(Fukushima(), 0.7, 0.4)

        @test set_default_algorithm!(Carlson) isa Carlson
        @test default_algorithm() isa Carlson
        @test K(0.4) == K(Carlson(), 0.4)
        @test F(0.7, 0.4) == F(Carlson(), 0.7, 0.4)
        @test sn(0.7, 0.4) == sn(Carlson(), 0.7, 0.4)
        @test am(0.7, 0.4) == am(Carlson(), 0.7, 0.4)

        @test set_default_algorithm!(AGM()) isa AGM
        @test default_algorithm() isa AGM
        @test K(0.4) == K(AGM(), 0.4)
        @test E(0.4) == E(AGM(), 0.4)
        @test F(0.7, 0.4) == F(Carlson(), 0.7, 0.4)
        @test sn(0.7, 0.4) == sn(Carlson(), 0.7, 0.4)
        @test am(0.7, 0.4) == am(Carlson(), 0.7, 0.4)
    finally
        set_default_algorithm!(original)
    end
end
