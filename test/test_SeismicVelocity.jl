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

    @test compute_wave_velocity(x_nd, (; wave = :Vp)) ≈ 8.1e11
    @test compute_wave_velocity(x_nd, (; wave = :Vs)) ≈ 4.5e11
    @test compute_wave_velocity(x_nd, (; wave = :VpVs)) ≈ 1.8
    @test compute_wave_velocity(x, (; T = 1.0f3, wave = :Vp)) === 8.1f3
    @test_throws "`wave` must be :Vp, :Vs or :VpVs" compute_wave_velocity(x, (; wave = :Vq))

    # Check that it works if we give a phase array
    MatParam = Array{MaterialParams, 1}(undef, 2)
    MatParam[1] = SetMaterialParams(;
        Name = "Mantle", Phase = 1, SeismicVelocity = ConstantSeismicVelocity()
    )

    MatParam[2] = SetMaterialParams(;
        Name = "Crust",
        Phase = 2,
        SeismicVelocity = PerpleX_LaMEM_Diagram(joinpath(@__DIR__, "test_data", "Peridotite_dry.in")),
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

    # reference values from an independent implementation of the same inclusion model
    # (solid: K = 79.44, G = 41.67, ρ = 3033, Vs = 3.71, Vp = 6.67; aspect ratio 0.5)
    @test [melt_correction_Takei(12.96, 79.44, 41.67, 2220.0, 3033.0, 6.67, 3.71, 0.2, 0.5)...] ≈ [3.0326, 5.6144] atol = 1.0e-4
    @test [melt_correction_Takei(12.96, 79.44, 41.67, 2220.0, 3033.0, 6.67, 3.71, 0.4, 0.5)...] ≈ [2.3215, 4.6277] atol = 1.0e-4
    @test [melt_correction_Takei(2.49, 79.44, 41.67, 1000.0, 3033.0, 6.67, 3.71, 0.2, 0.5)...] ≈ [3.1819, 5.6347] atol = 1.0e-4
    @test [melt_correction_Takei(12.96f0, 79.44f0, 41.67f0, 2220.0f0, 3033.0f0, 6.67f0, 3.71f0, 0.2f0, 0.5f0)...] ≈ [3.0326f0, 5.6144f0] atol = 1.0f-3
    @test melt_correction_Takei(12.96f0, 79.44f0, 41.67f0, 2220.0f0, 3033.0f0, 6.67f0, 3.71f0, 0.2f0, 0.5f0) isa NTuple{2, Float32}
    @test melt_correction(12.96f0, 79.44f0, 41.67f0, 2220.0f0, 3033.0f0, 6.67f0, 3.71f0, 0.2f0, 0.5f0) isa NTuple{2, Float32}

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

    # melt_correction (Takei 1998 framework moduli)
    Vp_cor, Vs_cor = melt_correction(26.0, 94.5, 61.0, 2802.0, 3198.0, 7.4, 4.36, 0.01, 0.84)
    @test [Vp_cor, Vs_cor] ≈ [7.2883533354949295, 4.265842133195837]
    # with equal densities ΔVs/Vs = ΛG ϕ/2; ΛG = (1 - μ_framework/G)/ϕ ≈ 4.44 for contiguity 0.84
    Vs_eqρ = melt_correction(26.0, 94.5, 61.0, 3000.0, 3000.0, 7.4, 4.36, 0.01, 0.84)[2]
    @test (1 - Vs_eqρ / 4.36) / 0.005 ≈ 4.443 atol = 1.0e-3
    @test melt_correction(26.0, 94.5, 61.0, 2802.0, 3198.0, 7.4, 4.36, 0.0, 0.84) == (7.4, 4.36)

    # porosity_correction: the shear part of the inclusion model of melt_correction_Takei
    ϕ_pore = 0.474 / (1 + 0.25 * 0.071)^5.989
    Vs_poro = porosity_correction(94.5, 61.0, 1000.0, 3198.0, 4.36, 0.25, 0.25)
    @test Vs_poro ≈ melt_correction_Takei(2.49, 94.5, 61.0, 1000.0, 3198.0, 7.4, 4.36, ϕ_pore, 0.25)[1]
    @test porosity_correction(94.5f0, 61.0f0, 1000.0f0, 3198.0f0, 4.36f0, 0.25f0, 0.25f0) isa Float32

    # anelastic_correction: invalid water mode throws a descriptive error
    @test_throws ArgumentError anelastic_correction(3, 4.36734, 5.0, 1250.0)

    # compute_wave_velocity for a phase without a SeismicVelocity parametrization -> 0
    mat_noVs = SetMaterialParams(; Name = "noVs", Phase = 1)
    @test compute_wave_velocity(mat_noVs, (1.0, 2.0)) == 0.0

    # correct_wavevelocities_phasediagrams: full pipeline on a real lookup table
    PD = PerpleX_LaMEM_Diagram(joinpath(@__DIR__, "test_data", "Peridotite_dry.in"))

    Vs_solid = copy(PD.solid_Vs.coefs)
    Vp_solid = copy(PD.solid_Vp.coefs)

    # default options (porosity + melt)
    PD_def = correct_wavevelocities_phasediagrams(PD)
    @test PD_def isa GeoParams.MaterialParameters.PhaseDiagrams.PhaseDiagram_LookupTable
    @test PD_def.Vp_uncorrected === PD.Vp
    @test PD_def.Vs_uncorrected === PD.Vs
    @test PD.solid_Vs.coefs == Vs_solid     # the input diagram is left untouched
    @test PD.solid_Vp.coefs == Vp_solid

    # melt only: P and S velocities are not swapped, and melt only reduces them
    PD_melt = correct_wavevelocities_phasediagrams(PD; apply_porosity_correction = false)
    @test all(PD_melt.Vp.coefs .> PD_melt.Vs.coefs)
    @test all(PD_melt.Vs.coefs .<= Vs_solid)
    @test all(PD_melt.Vp.coefs .<= Vp_solid)
    molten = PD.meltFrac.coefs .> 0
    @test PD_melt.Vs.coefs[.!molten] == Vs_solid[.!molten]
    i = findfirst(molten)
    @test (PD_melt.Vs.coefs[i], PD_melt.Vp.coefs[i]) == melt_correction_Takei(
        PD.melt_bulkModulus.coefs[i], PD.solid_bulkModulus.coefs[i], PD.solid_shearModulus.coefs[i],
        PD.meltRho.coefs[i], PD.rockRho.coefs[i], Vp_solid[i], Vs_solid[i], PD.meltFrac.coefs[i], 0.1
    )

    # weighted combination: porosity and melt corrections both start from the solid velocities
    PD_w = correct_wavevelocities_phasediagrams(PD; combine = :weighted)
    PD_poro = correct_wavevelocities_phasediagrams(PD; apply_melt_correction = false)
    @test PD_w.Vs.coefs[.!molten] == PD_poro.Vs.coefs[.!molten]
    @test PD_w.Vp.coefs[.!molten] == Vp_solid[.!molten]
    @test PD_w.Vs.coefs[molten] != PD_def.Vs.coefs[molten]
    @test_throws "`combine` must be :sequential or :weighted" correct_wavevelocities_phasediagrams(PD; combine = :sum)

    # anelasticity reduces the S-wave velocity of the solid
    PD_anel = correct_wavevelocities_phasediagrams(
        PD; apply_porosity_correction = false, apply_melt_correction = false, apply_anelasticity_correction = true
    )
    @test all(PD_anel.Vs.coefs .<= Vs_solid)
    @test PD_anel.Vs.coefs[end, end] < Vs_solid[end, end]

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
