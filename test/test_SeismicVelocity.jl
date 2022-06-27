using Test
using GeoParams
@testset "SeismicVelocity.jl" begin

    #Make sure structure is isbits
    x = ConstantSeismicVelocity()
    isbits(x)

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=1000km)

    # Define constant velocities 
    x1 = ConstantSeismicVelocity(; Vp=8.05km / s, Vs=3.5km / s)
    @test Value(x1.Vp) == 8.05km / s
    @test Value(x1.Vs) == 3.5km / s

    x1 = nondimensionalize(x1, CharUnits_GEO)
    @test NumValue(x1.Vp) ≈ 8.050000000000001e9
    @test NumValue(x1.Vs) ≈ 3.5e9

    # Compute
    @test compute_pwave_velocity(1.0, 1.0, x1) ≈ 8.050000000000001e9
    @test compute_swave_velocity(1.0, 1.0, x1) ≈ 3.5e9

    # Read Phase diagram interpolation object
    fname = "./test_data/Peridotite.in"
    PD_data = PerpleX_LaMEM_Diagram(fname)
    @test PD_data.Vp(1500, 1e7) ≈ 6.5290725233303935
    @test PD_data.Vs(1500, 1e7) ≈ 2.4874400647487658

    @test compute_pwave_velocity(1e7, 1500, PD_data) ≈ 6.5290725233303935
    @test compute_swave_velocity(1e7, 1500, PD_data) ≈ 2.4874400647487658

    # Do the same but non-dimensionalize the result
    CharDim = GEO_units()
    PD_data1 = PerpleX_LaMEM_Diagram(fname; CharDim=CharDim)

    rho_ND = PD_data1.Rho(
        nondimensionalize(1500K, CharDim), nondimensionalize(1e8 * Pa, CharDim)
    )
    Vp_ND = PD_data1.Vp(
        nondimensionalize(1500K, CharDim), nondimensionalize(1e8 * Pa, CharDim)
    )
    Vs_ND = PD_data1.Vs(
        nondimensionalize(1500K, CharDim), nondimensionalize(1e8 * Pa, CharDim)
    )

    # redimensionalize and check with value from original structure that did not use non-dimensionalization 
    @test ustrip(dimensionalize(rho_ND, kg / m^3, CharDim)) ≈ PD_data.Rho(1500, 1e8)
    @test ustrip(dimensionalize(Vp_ND, km / s, CharDim)) ≈ PD_data.Vp(1500, 1e8)
    @test ustrip(dimensionalize(Vs_ND, km / s, CharDim)) ≈ PD_data.Vs(1500, 1e8)

    # Test computation of velocity for the whole computational domain, using arrays 
    MatParam = Array{MaterialParams,1}(undef, 3)
    MatParam[1] = SetMaterialParams(;
        Name="Mantle",
        Phase=0,
        SeismicVelocity=PerpleX_LaMEM_Diagram("test_data/Peridotite.in"),
    )

    MatParam[2] = SetMaterialParams(;
        Name="Crust", Phase=1, SeismicVelocity=ConstantSeismicVelocity()
    )

    MatParam[3] = SetMaterialParams(;
        Name="UpperCrust",
        Phase=2,
        SeismicVelocity=ConstantSeismicVelocity(; Vp=10km / s, Vs=3km / s),
    )

    # test computing material properties
    Phases = ones(Int64, 400, 400) * 0
    Phases[:, 20:end] .= 1
    Phases[:, 300:end] .= 2

    Vp = zeros(size(Phases))
    Vs = zeros(size(Phases))
    T = ones(size(Phases))
    P = ones(size(Phases)) * 10

    compute_pwave_velocity!(Vp, Phases, P, T, MatParam)
    @test sum(Vp) / 400^2 ≈ 8.541562850000005

    # test computing material properties when we have PhaseRatios, instead of Phase numbers
    PhaseRatio = zeros(size(Phases)..., length(MatParam))
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz + 1)
        PhaseRatio[I] = 1.0
    end

    compute_swave_velocity!(Vs, PhaseRatio, P, T, MatParam)
    @test sum(Vs) / 400^2 ≈ 4.0739837

    Vp_cor, Vs_cor = melt_correction(
        26.0, 94.5, 61.0, 2802.0, 3198.0, 7.4, 4.36, 0.01, 0.15
    )
    @test [Vp_cor, Vs_cor] ≈ [7.336238790906285, 4.314027804335563]

    Vs_anel = anelastic_correction(0, 4.36734, 5.0, 1250.0)
    @test Vs_anel ≈ 4.1182815519599325
end
