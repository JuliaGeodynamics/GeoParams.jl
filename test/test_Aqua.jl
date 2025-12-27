if VERSION < v"1.13"
    using Aqua, Test, GeoParams

    ## Failing tests: hard to fix in the current state of the pkg
    # Aqua.test_unbound_args(GeoParams)
    # Aqua.test_ambiguities(GeoParams, color = true)
    # Aqua.test_piracies(GeoParams)

    @testset "Project extras" begin
        @test Aqua.test_project_extras(GeoParams).value
    end

    @testset "Undefined exports" begin
        @test Aqua.test_undefined_exports(GeoParams.MeltingParam).value
        @test Aqua.test_undefined_exports(GeoParams.Units).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.GravitationalAcceleration).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.PhaseDiagrams).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.ChemicalDiffusion).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.HeatCapacity).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.Conductivity).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.LatentHeat).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.RadioactiveHeat).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.Shearheating).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.Permeability).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.SeismicVelocity).value
        @test Aqua.test_undefined_exports(GeoParams.MaterialParameters.Density).value
        @test Aqua.test_undefined_exports(GeoParams.TASclassification).value
        @test Aqua.test_undefined_exports(GeoParams.ZirconAges).value
        @test Aqua.test_undefined_exports(GeoParams.Dislocation).value
        @test Aqua.test_undefined_exports(GeoParams.Diffusion).value
        @test Aqua.test_undefined_exports(GeoParams.GBS).value
        @test Aqua.test_undefined_exports(GeoParams.Peierls).value
        @test Aqua.test_undefined_exports(GeoParams.NonLinearPeierls).value
        @test Aqua.test_undefined_exports(GeoParams.Tables).value
        @test Aqua.test_undefined_exports(GeoParams.ChemicalDiffusion).value
    end

    @testset "Compats" begin
        @test !Aqua.test_deps_compat(
            GeoParams;
            check_julia = true,
            check_extras = false,
            check_weakdeps = true,
        ).anynonpass
        @test Aqua.test_stale_deps(GeoParams).value
    end

    @testset "Persistent tasks" begin
        Aqua.test_persistent_tasks(GeoParams)
    end
end
