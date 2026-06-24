using Test
using GeoParams

@testset "SeismicVelocity.jl" begin
    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 10km)

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

    @test UnitValue(compute_wave_velocity(x_nd, (; wave = :Vp))) ≈ 8.1e11
    @test UnitValue(compute_wave_velocity(x_nd, (; wave = :Vs))) ≈ 4.5e11
    @test UnitValue(compute_wave_velocity(x_nd, (; wave = :VpVs))) ≈ 1.8

    # Check that it works if we give a phase array
    MatParam = Array{MaterialParams, 1}(undef, 2)
    MatParam[1] = SetMaterialParams(;
        Name = "Mantle", Phase = 1, SeismicVelocity = ConstantSeismicVelocity()
    )

    MatParam[2] = SetMaterialParams(;
        Name = "Crust",
        Phase = 2,
        SeismicVelocity = PerpleX_LaMEM_Diagram("test_data/Peridotite_dry.in"),
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

    args = (; T = T, P = P, wave = :Vp)
    compute_wave_velocity!(Vp, Mat_tup, Phases, args)

    args = (; T = T, P = P, wave = :Vs)
    compute_wave_velocity!(Vs, Mat_tup, Phases, args)

    args = (; T = T, P = P, wave = :VpVs)
    compute_wave_velocity!(VpVs, Mat_tup, Phases, args)

    @test Vp[1] == 8.1e3
    @test Vp[1, 1, end] ≈ 5.500887338991992
    @test Vs[1] == 4.5e3
    @test Vs[1, 1, end] ≈ 2.68

    @test VpVs[1] ≈ 1.8
    @test VpVs[1, 1, end] ≈ 2.05

    # NOTE: This will be made obsolete by melt_correction_Takei
    #    Vp_cor, Vs_cor = melt_correction(
    #        26.0, 94.5, 61.0, 2802.0, 3198.0, 7.4, 4.36, 0.01, 0.15
    #    )
    #    @test [Vp_cor, Vs_cor] ≈  [7.331657177397843, 4.314027804335563]

    # NOTE: This will be made obsolete by melt_correction_Takei
    #  Vs_cor = porosity_correction(
    #       94.5, 61.0, 1000.0, 3198.0, 4.36, 0.25, 0.25
    #  )
    #  @test [Vs_cor] ≈ [2.226167083352012]

    Vs_anel = anelastic_correction(0, 4.36734, 5.0, 1250.0)
    @test Vs_anel ≈ 4.343623758644558

    # testing the new seismic velocity correction for partial melt
    ρL = 2000.0
    ρS = 3300.0
    Vs0 = 3000.0
    Vp0 = 6000.0
    α = 0.4
    ϕ = 0.7
    Kb_S = 250.0
    Ks_S = 162.0
    Kb_L = 200.0
    R = 0.1

    melt_correction_Takei(Kb_L, Kb_S, Ks_S, ρL, ρS, Vp0, Vs0, ϕ, α)

    ϕ_vec = 0:0.01:1
    Vs_new = zero(ϕ_vec)
    Vp_new = zero(ϕ_vec)

    for i in eachindex(ϕ_vec)
        Vs_new[i], Vp_new[i] = melt_correction_Takei(
            Kb_L, Kb_S, Ks_S, ρL, ρS, Vp0, Vs0, ϕ_vec[i], α
        )
    end

    @test Vs_new[10] ≈ 2750.713407744307
    @test Vp_new[10] ≈ 5792.937134183798

    # ConstantSeismicVelocity vararg constructor
    x_vararg = ConstantSeismicVelocity(8100m / s, 4500m / s)
    @test Value(x_vararg.Vp) ≈ 8100m / s

    # anelastic_correction: water = 1 (damp) and water = 2 (wet)
    Vs_dry = anelastic_correction(0, 4.36734, 5.0, 1250.0)
    Vs_damp = anelastic_correction(1, 4.36734, 5.0, 1250.0)
    Vs_wet = anelastic_correction(2, 4.36734, 5.0, 1250.0)
    @test Vs_damp < 4.36734   # correction reduces velocity
    @test Vs_wet < 4.36734
    @test Vs_dry > Vs_damp    # more water -> larger correction

    # melt_correction_Takei: ϕ = 0 -> velocities unchanged
    Vs_nomelt, Vp_nomelt = melt_correction_Takei(Kb_L, Kb_S, Ks_S, ρL, ρS, Vp0, Vs0, 0.0, α)
    @test Vs_nomelt ≈ Vs0 atol = 1.0
    @test Vp_nomelt ≈ Vp0 atol = 1.0

    # melt_correction_Takei: α = 1.0 -> R_func is NaN -> correction is skipped (else branch)
    Vs_nan, Vp_nan = melt_correction_Takei(Kb_L, Kb_S, Ks_S, ρL, ρS, Vp0, Vs0, 0.7, 1.0)
    @test Vs_nan == Vs0
    @test Vp_nan == Vp0

    # melt_correction_Takei: extreme contrast clamps negative velocities to zero
    Vs_clamp, Vp_clamp = melt_correction_Takei(1.0, 250.0, 162.0, 2000.0, 3300.0, 6000.0, 3000.0, 0.8, 0.1)
    @test Vs_clamp == 0.0
    @test Vp_clamp == 0.0

    # melt_correction (Takei 1998, restored): matches the published reference values
    Vp_cor, Vs_cor = melt_correction(26.0, 94.5, 61.0, 2802.0, 3198.0, 7.4, 4.36, 0.01, 0.15)
    @test [Vp_cor, Vs_cor] ≈ [7.331657177397843, 4.314027804335563]

    # porosity_correction (restored): guarded against unphysical negative velocities
    Vs_poro = porosity_correction(94.5, 61.0, 1000.0, 3198.0, 4.36, 0.25, 0.25)
    @test Vs_poro ≥ 0.0

    # anelastic_correction: invalid water mode throws a descriptive error
    @test_throws ArgumentError anelastic_correction(3, 4.36734, 5.0, 1250.0)

    # compute_wave_velocity for a phase without a SeismicVelocity parametrization -> 0
    mat_noVs = SetMaterialParams(; Name = "noVs", Phase = 1)
    @test compute_wave_velocity(mat_noVs, (1.0, 2.0)) == 0.0

    # correct_wavevelocities_phasediagrams: full pipeline on a real lookup table
    PD = PerpleX_LaMEM_Diagram("test_data/Peridotite_dry.in")

    # default options (anelasticity + Takei melt + guarded porosity)
    PD_def = correct_wavevelocities_phasediagrams(PD)
    @test PD_def isa GeoParams.MaterialParameters.PhaseDiagrams.PhaseDiagram_LookupTable
    @test PD_def.Vp_uncorrected !== nothing
    @test PD_def.Vs_uncorrected !== nothing

    # exercise the legacy (non-Takei) melt_correction branch as well
    PD_legacy = correct_wavevelocities_phasediagrams(
        PD; apply_porosity_correction = false, melt_correction_takei = false, water = 2
    )
    @test PD_legacy isa GeoParams.MaterialParameters.PhaseDiagrams.PhaseDiagram_LookupTable

    # corrected velocities are stored as evaluatable interpolation objects
    T0, P0 = PD.solid_Vs.T0, PD.solid_Vs.P0
    @test PD_def.Vp(T0, P0) isa Real
    @test PD_def.Vs(T0, P0) isa Real
end
