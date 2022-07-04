using Test
using GeoParams

@testset "SeismicVelocity.jl" begin
    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=10km)

    # Constant seismic velocity capacity
    x = ConstantSeismicVelocity()
    @test isbits(x) == true
    info = param_info(x)

    x_nd = x
    x_nd = nondimensionalize(x_nd, CharUnits_GEO)

    @test Value(x.Vp) ≈ 8.1km / s
    @test Value(x.Vs) ≈ 4.5km / s
    @test UnitValue(x_nd.Vp) ≈ 8.1e11
    @test UnitValue(x_nd.Vs) ≈ 4.5e11

    @test UnitValue(compute_wave_velocity(x_nd, (; wave=:Vp))) ≈ 8.1e11
    @test UnitValue(compute_wave_velocity(x_nd, (; wave=:Vs))) ≈ 4.5e11
    @test UnitValue(compute_wave_velocity(x_nd, (; wave=:VpVs))) ≈ 1.8

    # Check that it works if we give a phase array
    MatParam = Array{MaterialParams,1}(undef, 2)
    MatParam[1] = SetMaterialParams(;
        Name="Mantle", Phase=1, SeismicVelocity=ConstantSeismicVelocity()
    )

    MatParam[2] = SetMaterialParams(;
        Name="Crust",
        Phase=2,
        SeismicVelocity=PerpleX_LaMEM_Diagram("test_data/Peridotite_dry.in"),
    )

    Mat_tup = Tuple(MatParam)

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    Phases[:, :, 20:end] .= 2

    Vp = zeros(size(Phases))
    Vs = zeros(size(Phases))
    VpVs = zeros(size(Phases))
    T = ones(size(Phases)) * 1500
    P = zeros(size(Phases))

    args = (; T=T, P=P, wave=:Vp)
    compute_wave_velocity!(Vp, Mat_tup, Phases, args)

    args = (; T=T, P=P, wave=:Vs)
    compute_wave_velocity!(Vs, Mat_tup, Phases, args)

    args = (; T=T, P=P, wave=:VpVs)
    compute_wave_velocity!(VpVs, Mat_tup, Phases, args)

    @test Vp[1] == 8.1
    @test Vp[1, 1, end] ≈ 5.500887338991992
    @test Vs[1] == 4.5
    @test Vs[1, 1, end] ≈ 2.68 

    @test VpVs[1] ≈ 1.8
    @test VpVs[1, 1, end] ≈ 2.05

    Vp_cor, Vs_cor = melt_correction(
        26.0, 94.5, 61.0, 2802.0, 3198.0, 7.4, 4.36, 0.01, 0.15
    )
    @test [Vp_cor, Vs_cor] ≈  [7.331657177397843, 4.314027804335563]

    Vs_cor = porosity_correction(
         94.5, 61.0, 1000.0, 3198.0, 4.36, 0.25, 0.25
    )
    @test [Vs_cor] ≈ [2.226167083352012]


    Vs_anel = anelastic_correction(0, 4.36734, 5.0, 1250.0)
    @test Vs_anel ≈ 4.343623758644558


    # Correct phase diagrams for melt, anelasticity and porosity

    PD= correct_wavevelocities_phasediagrams(PerpleX_LaMEM_Diagram("test_data/Peridotite_dry.in"), apply_melt_correction=false, apply_porosity_correction=false, apply_anelasticity_correction=false)
    @test PD.Vs(1500,5e8) ≈ 4.329112661008007
    @test PD.Vp(1500,5e8) ≈ 7.376899430314534

    PD= correct_wavevelocities_phasediagrams(PerpleX_LaMEM_Diagram("test_data/Peridotite_dry.in"), apply_melt_correction=false, apply_porosity_correction=false, apply_anelasticity_correction=true)
    @test PD.Vs(1500,5e8) ≈ 2.6516022088436753
    @test PD.Vp(1500,5e8) ≈ 7.376899430314534

    PD= correct_wavevelocities_phasediagrams(PerpleX_LaMEM_Diagram("test_data/Peridotite_dry.in"), apply_melt_correction=false, apply_porosity_correction=true, apply_anelasticity_correction=true)
    @test PD.Vs(1500,5e8) ≈ 2.647426620329417
    @test PD.Vp(1500,5e8) ≈ 7.376899430314534

    PD= correct_wavevelocities_phasediagrams(PerpleX_LaMEM_Diagram("test_data/Peridotite_dry.in"), apply_melt_correction=true, apply_porosity_correction=true, apply_anelasticity_correction=true)
    @test PD.Vs(1500,1e8) ≈ 1.4709017701696432
    @test PD.Vp(1500,5e8) ≈ 6.774831213555857
end
