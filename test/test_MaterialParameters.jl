using Test
using GeoParams

@testset "MaterialParameters" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 1000km)

    # Define a struct for a first phase
    Phase1 = SetMaterialParams(;
        Name = "test1",
        Phase = 22,
        CreepLaws = (PowerlawViscous(), LinearViscous(; η = 1.0e21Pa * s)),
        Gravity = ConstantGravity(; g = 11.0m / s^2),
        Density = ConstantDensity(),
        CharDim = CharUnits_GEO,
    )

    @test Phase1.Density[1].ρ.val ≈ 2.9e-16 rtol = 1.0e-6
    @test Phase1.Gravity[1].g * 1.0 ≈ 1.1e19 rtol = 1.0e-6
    @test Phase1.CreepLaws[1].η0 * 1.0 == 0.1
    @test Phase1.CreepLaws[2].η * 1.0 == 100.0

    Phase2 = SetMaterialParams(;
        Name = "test1",
        Phase = 2,
        CreepLaws = (LinearViscous(), LinearViscous(; η = 1.0e21Pa * s)),
        Density = ConstantDensity(),
    )

    # Non-dimensionalize all values:
    Phase2 = nondimensionalize(Phase2, CharUnits_GEO)
    @test Phase2.CreepLaws[2].η.val ≈ 100.0
    @test Phase2.Density[1].ρ.val ≈ 2.9e-16

    # Dimensionalize all values again:
    Phase2 = dimensionalize(Phase2, CharUnits_GEO)
    @test Phase2.Density[1].ρ.val ≈ 2900
    @test Phase2.CreepLaws[2].η.val ≈ 1.0e21

    # Create array with several phases
    MatParam = Vector{MaterialParams}(undef, 2)
    Phase = 1
    MatParam[Phase] = SetMaterialParams(;
        Name = "Upper Crust",
        Phase = Phase,
        CreepLaws = (PowerlawViscous(), LinearViscous(; η = 1.0e23Pa * s)),
        Density = ConstantDensity(; ρ = 2900kg / m^3),
    )
    Phase = 2
    MatParam[Phase] = SetMaterialParams(;
        Name = "Lower Crust",
        Phase = Phase,
        CreepLaws = (PowerlawViscous(; n = 5.0), LinearViscous(; η = 1.0e21Pa * s)),
        Density = PT_Density(; ρ0 = 3000kg / m^3),
    )

    @test MatParam[2].Density[1].α.val ≈ 3.0e-5
    @test MatParam[2].CreepLaws[1].n.val == 5.0

    # test adding phase Diagrams
    MatParam = Array{MaterialParams, 1}(undef, 1)
    Phase = 1
    MatParam[Phase] = SetMaterialParams(;
        Name = "Mantle",
        Phase = Phase,
        CreepLaws = (PowerlawViscous(), LinearViscous(; η = 1.0e23Pa * s)),
        Density = PerpleX_LaMEM_Diagram("test_data/Peridotite_dry.in"),
        CharDim = GEO_units(),
    )

    @test MatParam[1].Density[1].Rho(1, 10) ≈ 3.180692806182894e-18

    MatParam = Array{MaterialParams, 1}(undef, 1)
    Phase = 1
    MatParam[Phase] = SetMaterialParams(;
        Name = "Mantle",
        Phase = Phase,
        CreepLaws = (PowerlawViscous(), LinearViscous(; η = 1.0e23Pa * s)),
        Density = PT_Density(; ρ0 = 3000kg / m^3),
        CharDim = nothing,
    )

    @test isdimensional(MatParam[1].Density[1])

    # Base.show for NTuple{N, MaterialParams}
    buf = IOBuffer()
    show(buf, (Phase1, Phase2))
    @test !isempty(String(take!(buf)))

    # ConvField with multiple creep laws (NTuple path)
    Phase_multi = SetMaterialParams(;
        Name = "multi_creep", Phase = 1,
        CreepLaws = (LinearViscous(), PowerlawViscous()),
    )
    @test length(Phase_multi.CreepLaws) == 2

    # ConvField error: maxAllowedFields exceeded for density
    @test_throws Exception SetMaterialParams(;
        Name = "bad_density", Phase = 1, Density = (ConstantDensity(), ConstantDensity())
    )

    # set_gravity: explicit Gravity supplied (not auto-set from Density)
    Phase_grav = SetMaterialParams(;
        Name = "with_gravity", Phase = 1,
        Density = ConstantDensity(), Gravity = ConstantGravity(; g = 11.0m / s^2),
    )
    @test Phase_grav.Gravity[1].g.val ≈ 11.0

    # nondimensionalize_phase: invalid CharDim triggers an error
    @test_throws Exception SetMaterialParams(; Name = "bad_chardim", Phase = 1, CharDim = "invalid")
end
