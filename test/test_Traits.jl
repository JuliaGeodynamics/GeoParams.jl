using Test, GeoParams

@testset "Rheology traits" begin
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
            CompositeRheology = CompositeRheology(v2, el),
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
        SetMaterialParams(;
            CompositeRheology = CompositeRheology(el, v3),
        ),
    )

    # test basic cases
    for r in (v1, v2, pl)
        @test islinear(r) isa NonLinearRheologyTrait
    end
    for r in (v3, el)
        @test islinear(r) isa LinearRheologyTrait
    end
    @test_throws ArgumentError islinear("potato")

    # test composite cases
    @test islinear(v1, v2) isa NonLinearRheologyTrait
    @test islinear(v1, v3) isa NonLinearRheologyTrait
    @test islinear(v3, el) isa LinearRheologyTrait

    @test islinear(tuple(v1, v2)) isa NonLinearRheologyTrait
    @test islinear(tuple(v1, v3)) isa NonLinearRheologyTrait
    @test islinear(tuple(v3, el)) isa LinearRheologyTrait

    @test islinear(CompositeRheology(v1, v2, v3)) isa NonLinearRheologyTrait
    @test islinear(CompositeRheology(el, v3))     isa LinearRheologyTrait

    # # test MaterialParams cases
    @test islinear(r1[1]) isa NonLinearRheologyTrait
    @test islinear(r1[2]) isa NonLinearRheologyTrait
    @test islinear(r1[3]) isa LinearRheologyTrait
    @test islinear(r1)    isa NonLinearRheologyTrait

    @test islinear(r2[1]) isa LinearRheologyTrait
    @test islinear(r2[2]) isa LinearRheologyTrait
    @test islinear(r2[3]) isa LinearRheologyTrait
    @test islinear(r2)    isa LinearRheologyTrait
end

@testset "Density traits" begin
    v1 = ConstantDensity()
    v2 = PT_Density()
    d1 = SetMaterialParams(; Density = v1)
    d2 = SetMaterialParams(; Density = v2)
    r1 = (d1, d2)
    r2 = (d1, d1)

    @test isconstant(v1) isa ConstantDensityTrait
    @test isconstant(v2) isa NonConstantDensityTrait
    @test isconstant(d1) isa ConstantDensityTrait
    @test isconstant(d2) isa NonConstantDensityTrait
    @test isconstant(r1) isa NonConstantDensityTrait
    @test isconstant(r2) isa ConstantDensityTrait
end