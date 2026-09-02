using Test, GeoParams

@testset "Computations.jl" begin

    # Parameters whose value depends on arguments have no single-argument
    # `compute_*` method, so an argument-free call has nowhere to dispatch.
    @testset "argument-free calls raise MethodError" begin
        for (f, x) in (
                (compute_density, PT_Density()),
                (compute_density, Compressible_Density()),
                (compute_conductivity, TP_Conductivity()),
                (compute_meltfraction, MeltingParam_4thOrder()),
                (compute_permeability, HazenPermeability()),
                (compute_permeability, PowerLawPermeability()),
                (compute_permeability, CarmanKozenyPermeability()),
            )
            @test_throws MethodError f(x)
        end
    end

    # Parameters with an argument-free value keep their single-argument method.
    @testset "constant parameters still evaluate without arguments" begin
        @test compute_density(ConstantDensity()) == 2900.0
        @test compute_conductivity(ConstantConductivity()) == 3.0
        @test compute_permeability(ConstantPermeability()) == 1.0e-12
        @test compute_heatcapacity(ConstantHeatCapacity()) == 1050.0
    end

    @testset "argument-carrying calls are unaffected" begin
        args = (; P = 1.0e9, T = 1000.0, ϕ = 0.1)

        @test compute_density(PT_Density(), args) ≈ 5713.0
        @test compute_conductivity(TP_Conductivity(), args) ≈ 1.6201114206128133
        @test compute_permeability(CarmanKozenyPermeability(), args) ≈ 1000.0

        rheologies = (
            SetMaterialParams(; Phase = 1, Density = PT_Density()),
            SetMaterialParams(; Phase = 2, Density = ConstantDensity()),
        )

        @test compute_density(rheologies, 1, args) ≈ 5713.0
        @test compute_density(rheologies, 2, args) ≈ 2900.0
        @test compute_density(rheologies, (0.5, 0.5), args) ≈ 4306.5

        ρ = zeros(3)
        compute_density!(
            ρ, rheologies, [1, 2, 1], (; P = fill(1.0e9, 3), T = fill(1000.0, 3))
        )
        @test ρ ≈ [5713.0, 2900.0, 5713.0]
    end

end
