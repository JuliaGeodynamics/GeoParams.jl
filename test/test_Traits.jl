using Test, GeoParams

@testset "Traits" begin
    # Define a range of rheological components
    v1 = DiffusionCreep()
    v2 = DislocationCreep()
    v3 = LinearViscous()
    el = ConstantElasticity()           # elasticity
    pl = DruckerPrager()        # plasticity which ends up with the same yield stress as pl3
    # CompositeRheologies
    c = CompositeRheology(v1, v2, v3)
    r1 = (
        SetMaterialParams(;
            CompositeRheology = CompositeRheology(v1, v2, v3),
        ),
        SetMaterialParams(;
            CompositeRheology = CompositeRheology(el, v3),
        ),
    )
    r2 = (
        SetMaterialParams(;
            CompositeRheology = CompositeRheology(el, v3),
        ),
        SetMaterialParams(;
            CompositeRheology = CompositeRheology(el, v3),
        ),
    )

    # test basic cases
    for r in (v1, v2, pl)
        @test islinear(r) isa NonLinearRheology
    end
    for r in (v3, el)
        @test islinear(r) isa LinearRheology
    end
    @test_throws ArgumentError islinear("potato")

    # test composite cases
    @test islinear(v1, v2) isa NonLinearRheology
    @test islinear(v1, v3) isa NonLinearRheology
    @test islinear(v3, el) isa LinearRheology

    @test islinear(tuple(v1, v2)) isa NonLinearRheology
    @test islinear(tuple(v1, v3)) isa NonLinearRheology
    @test islinear(tuple(v3, el)) isa LinearRheology

    @test islinear(CompositeRheology(v1, v2, v3)) isa NonLinearRheology
    @test islinear(CompositeRheology(el, v3)) isa LinearRheology

    # test MaterialParams cases
    @test islinear(r1[1]) isa NonLinearRheology
    @test islinear(r1[2]) isa LinearRheology
    @test islinear(r1)    isa NonLinearRheology

    @test islinear(r2[1]) isa LinearRheology
    @test islinear(r2[2]) isa LinearRheology
    @test islinear(r2)    isa LinearRheology
end
