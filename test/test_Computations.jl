using Test, GeoParams
import ForwardDiff

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

        @test compute_density(PT_Density(), args) ≈ 5736.764050000001
        @test compute_conductivity(TP_Conductivity(), args) ≈ 1.6201114206128133
        @test compute_permeability(CarmanKozenyPermeability(), args) ≈ 1000.0

        rheologies = (
            SetMaterialParams(; Phase = 1, Density = PT_Density()),
            SetMaterialParams(; Phase = 2, Density = ConstantDensity()),
        )

        @test compute_density(rheologies, 1, args) ≈ 5736.764050000001
        @test compute_density(rheologies, 2, args) ≈ 2900.0
        @test compute_density(rheologies, (0.5, 0.5), args) ≈ 4318.382025000001

        ρ = zeros(3)
        compute_density!(
            ρ, rheologies, [1, 2, 1], (; P = fill(1.0e9, 3), T = fill(1000.0, 3))
        )
        @test ρ ≈ [5736.764050000001, 2900.0, 5736.764050000001]
    end

    @testset "array phase lookup preserves material identity" begin
        for F in (Float64, Float32), count in (1, 4, 16), ids in (
                    ntuple(identity, count), ntuple(i -> count - i + 1, count),
                    ntuple(i -> 10i, count), ntuple(i -> i == 1 ? 2 : i, count),
                )
            materials = ntuple(
                i -> SetMaterialParams(;
                    Phase = ids[i], Density = ConstantDensity(ρ = F(2700 + 20i) * kg / m^3),
                    Conductivity = ConstantConductivity(k = F(i) * Watt / m / K)
                ), count
            )
            phases = reshape([ids[mod(i - 1, count) + 1] for i in 1:(4count)], count, 4)
            args = (; T = fill(F(1000), size(phases)), P = F(1.0e9))
            out = zeros(F, size(phases))
            # Scalar lookup is the reference, including first-match duplicate IDs.
            expected = map(p -> compute_density(materials, p, (; T = F(1000), P = F(1.0e9))), phases)
            compute_density!(out, materials, phases, args)
            @test out ≈ expected
            compute_conductivity!(out, materials, phases, args)
            @test out ≈ map(p -> compute_conductivity(materials, p, (; T = F(1000))), phases)
            destination = view(out, :, 1:2)
            phase_view = view(phases, :, 1:2)
            compute_density!(destination, materials, phase_view, (; T = F(1000), P = F(1.0e9)))
            @test destination ≈ expected[:, 1:2]
            @test_throws ArgumentError compute_density!(out, materials, fill(Int64(0), size(phases)), args)
            @test_throws ArgumentError compute_density!(out, materials, fill(Int64(maximum(ids) + 1), size(phases)), args)
            @test compute_density!(F[], materials, Int64[], (; T = F(1000))) === nothing
        end
        mixed = (
            SetMaterialParams(; Phase = 1, Density = PT_Density()),
            SetMaterialParams(; Phase = 2, Density = ConstantDensity()),
        )
        out = zeros(4)
        compute_density!(out, mixed, [1, 2, 2, 1], (; T = 1000.0, P = 1.0e9))
        @test out ≈ [5736.764050000001, 2900.0, 2900.0, 5736.764050000001]
    end

    @testset "array phase lookup supports automatic differentiation" begin
        materials = ntuple(i -> SetMaterialParams(; Phase = i, Density = PT_Density()), 16)
        phases = collect(Int64, 1:16)
        temperatures = collect(1000.0:1015.0)
        total_density = function (temperature)
            out = similar(temperature)
            compute_density!(out, materials, phases, (; T = temperature, P = 1.0e9))
            return sum(out)
        end
        expected = [
            ForwardDiff.derivative(t -> compute_density(materials, phase, (; T = t, P = 1.0e9)), temperature)
                for (phase, temperature) in zip(phases, temperatures)
        ]
        @test ForwardDiff.gradient(total_density, temperatures) ≈ expected
    end

end
