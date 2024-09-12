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
            CompositeRheology = CompositeRheology(pl, el, v3),
        ),
        SetMaterialParams(;
            CompositeRheology = CompositeRheology(v2, el),
        ),
        SetMaterialParams(;
            CompositeRheology = CompositeRheology(el, v3),
        ),
        SetMaterialParams(;
            CompositeRheology = CompositeRheology(v1, v2),
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

    # viscoelastic rheology traits
    # test basic cases
    for r in (v1, v2, v3, pl)
        @test isviscoelastic(r) isa NonElasticRheologyTrait
    end
    for r in (el,)
        @test isviscoelastic(r) isa ElasticRheologyTrait
    end
    @test_throws ArgumentError isviscoelastic("potato")

    # test composite cases
    @test isviscoelastic(v1, v2) isa NonElasticRheologyTrait
    @test isviscoelastic(v1, v3) isa NonElasticRheologyTrait
    @test isviscoelastic(v3, el) isa ElasticRheologyTrait

    @test isviscoelastic(tuple(v1, v2)) isa NonElasticRheologyTrait
    @test isviscoelastic(tuple(v1, v3)) isa NonElasticRheologyTrait
    @test isviscoelastic(tuple(v3, el)) isa ElasticRheologyTrait

    @test isviscoelastic(CompositeRheology(v1, v2, v3, pl)) isa NonElasticRheologyTrait
    @test isviscoelastic(CompositeRheology(el, v3, pl))     isa ElasticRheologyTrait

    # test MaterialParams cases
    @test isviscoelastic(r1[1]) isa ElasticRheologyTrait
    @test isviscoelastic(r1[2]) isa ElasticRheologyTrait
    @test isviscoelastic(r1[3]) isa ElasticRheologyTrait
    @test isviscoelastic(r1[4]) isa NonElasticRheologyTrait

    # test get_G and get_Kb
    r = CompositeRheology(v1, v2, v3)
    @test GeoParams.get_G(isviscoelastic(r), r) == 0
    @test GeoParams.get_Kb(isviscoelastic(r), r) == Inf

    ## linear rheology traits
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

    # test MaterialParams cases
    @test islinear(r1[1]) isa NonLinearRheologyTrait
    @test islinear(r1[2]) isa NonLinearRheologyTrait
    @test islinear(r1[3]) isa LinearRheologyTrait
    @test islinear(r1)    isa NonLinearRheologyTrait

    @test islinear(r2[1]) isa LinearRheologyTrait
    @test islinear(r2[2]) isa LinearRheologyTrait
    @test islinear(r2[3]) isa LinearRheologyTrait
    @test islinear(r2)    isa LinearRheologyTrait

    ## plastic rheology traits
    # test basic cases
    for r in (v1, v2, v3, el)
        @test isplasticity(r) isa NonPlasticRheologyTrait
    end
    @test isplasticity(pl) isa PlasticRheologyTrait
    @test_throws ArgumentError isplasticity("potato")

    # test composite cases
    @test isplasticity(v1, v2) isa NonPlasticRheologyTrait
    @test isplasticity(pl, el) isa PlasticRheologyTrait
    @test isplasticity(tuple(v1, v2)) isa NonPlasticRheologyTrait
    @test isplasticity(tuple(pl, el)) isa PlasticRheologyTrait
    
    @test isplasticity(CompositeRheology(v1, v2, v3)) isa NonPlasticRheologyTrait
    @test isplasticity(CompositeRheology(el, pl, v1)) isa PlasticRheologyTrait

    # test MaterialParams cases
    @test isplasticity(r1) isa PlasticRheologyTrait
    @test isplasticity(r2) isa NonPlasticRheologyTrait
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
    @test_throws ArgumentError isconstant("potato")
end

