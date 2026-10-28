using Test
using GeoParams
import GeoParams.NonLinearPeierls

@testset "NonLinearPeierlsCreepLaws" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1.0Pa * s, length = 1.0m)

    # Define a linear viscous creep law ---------------------------------
    x1 = NonLinearPeierlsCreep()
    @test isbits(x1)
    str = sprint(show, x1)
    @test occursin("NonLinearPeierlsCreep:", str)
    @test Value(x1.n) == 2.0
    @test Value(x1.q) == 1.0
    @test Value(x1.o) == 0.5
    @test Value(x1.A) == 5.7e11MPa^(-2.0) * s^(-1.0)
    @test repr("text/plain", x1) isa String

    # perform a computation with the peierls creep laws
    # Calculate EpsII, using a set of pre-defined values
    CharDim = GEO_units(;
        length = 1000.0km, viscosity = 1.0e19Pa * s, stress = 100.0MPa, temperature = 1000.0C
    )
    EpsII = GeoUnit(1.0s^-1.0)
    TauII = GeoUnit(1.0e3MPa)
    T = GeoUnit(473.0C)

    # compute a pure non linear peierls creep rheology
    p = SetNonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010; n = 3.1)
    @test p.n.val == 3.1

    # compute a pure non linear peierls creep rheology
    p = SetNonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010)
    T = 200.0 + 273.15
    args = (; T = T)
    TauII = 1.0e9
    ε = compute_εII(p, TauII, args)
    @test ε ≈ 2.3708323119196878e-17

    # same but with removing the tensor correction
    ε_notensor = compute_εII(remove_tensor_correction(p), TauII, args)
    @test ε_notensor ≈ 2.3087681699131527e-22

    # test with arrays
    τII_array = ones(10) * 1.0e9
    ε_array = similar(τII_array)
    T_array = ones(size(τII_array)) * (200.0 + 273.15)
    args_array = (; T = T_array)
    compute_εII!(ε_array, p, τII_array, args_array)
    @test ε_array[1] ≈ ε

    # compute when args are scalars
    compute_εII!(ε_array, p, τII_array, args)
    @test ε_array[1] ≈ ε

    # wet olivine, stress-strainrate curve
    # p    = SetNonLinearPeierlsCreep("Dry Olivine | Mei et al. (2010)")
    εII = exp10.(-22:0.5:-12)
    τII = exp10.(9:0.05:10)               # preallocate array
    T = 200.0 + 273.15
    args = (; T = T)
    εII = zero(τII)
    compute_εII!(εII, p, τII, args)
    eta_array1 = @. 0.5 * τII / εII

    # test overriding the default values
    # a = SetNonLinearPeierlsCreep("Dry Olivine | Mei et al. (2010)"; E=475.0kJ / mol)
    # @test Value(a.E) == 475.0kJ / mol

    # Do some basic checks on all creeplaws in the DB
    CharDim = GEO_units()
    creeplaw_list = nonlinearpeierls_law_list() # all creeplaws in database
    for fun in creeplaw_list
        p = SetNonLinearPeierlsCreep(fun)                # original creep law
        p_nd = nondimensionalize(p, CharDim)                # non-dimensionalized
        @test isfinite(NumValue(p_nd.A))         # extreme n can over- or underflow the scaling
        @test p_nd == SetNonLinearPeierlsCreep(fun, CharDim) # check that the non-dimensionalized version is the same as the original
        p_dim = dimensionalize(p, CharDim)                   # dimensionalized

        # Check that values are the same after non-dimensionalisation & dimensionalisation
        for field in fieldnames(typeof(p_dim))
            val_original = getfield(p, field)
            val_final = getfield(p_dim, field)
            if isa(val_original, GeoUnit)
                @test Value(val_original) == Value(val_final)
            end
        end

        # Perform computations with the rheology
        args = (T = 273 + 200.0, d = 100.0e-6, τII_old = 1.0e6)
        ε = 1.0e-15
        τ_ini = 0.3 * NumValue(p.TauP) / p.FT               # initial guess for stress iterations, below the Peierls stress
        τ = Peierls_stress_iterations(p, τ_ini, ε, args)
        ε_test = compute_εII(p, τ; args...)
        @test ε ≈ ε_test
    end

    # check that derivative calculation is correct
    p = SetNonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010)
    TauII = 1.0e9
    args = (; T = 273 + 200.0)
    depsdtau = dεII_dτII(p, TauII; args...)
    @test depsdtau ≈ 5.633292966430465e-25

    # param_info resolves the builder from the stored Name (find_creep_law)
    info_named = param_info(SetNonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010))
    @test info_named isa MaterialParamsInfo
    @test !isempty(info_named.BibTex_Reference)
    @test isempty(param_info(NonLinearPeierlsCreep()).BibTex_Reference)                # empty-name early return
    @test isempty(param_info(NonLinearPeierlsCreep(; Name = "nope_xyz")).BibTex_Reference)  # unknown -> fallback

    # nondimensional unit transform path
    @test GeoParams.Transform_NonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010, GEO_units()) isa AbstractCreepLaw

    # Mei et al. (2010), eq. (3): ε̇ = A_P σ² exp(-E/RT (1 - √(σ/σP))), A_P = 1.4e-7 MPa⁻² s⁻¹, σ the
    # differential stress in MPa (σ = √3 τII) and ε̇ the axial strain rate (ε̇II = ε̇ / (2/√3)).
    σ, T = 1.0e3, 873.0
    a = SetNonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010)
    @test compute_εII(a, σ * 1.0e6 / sqrt(3); T) ≈
        1.4e-7 * σ^2 * exp(-320.0e3 / (NumValue(a.R) * T) * (1 - sqrt(σ / 5.9e3))) / (2 / sqrt(3)) rtol = 1.0e-6
end

# The (σ/μ)² laws are stored as NonLinearPeierlsCreep with n = 2 and A/μ² (σ in MPa)
@testset "NonLinearPeierls database vs. printed laws" begin
    σ, T = 1000.0, 700.0
    # Burdette & Hirth (2022): differential stress, μ = 35 GPa
    p = SetNonLinearPeierlsCreep(NonLinearPeierls.antigorite_Burdette_2022)
    ε = compute_εII(p, σ * 1.0e6 / p.FT; T) * p.FE
    @test ε ≈ exp(-0.624) * (σ / 35.0e3)^2 * exp(-86.3e3 / (NumValue(p.R) * T) * (1 - σ / 2420)^1.18) rtol = 1.0e-6

    # Hufford et al. (2026): σ = 2τ, ε̇ = γ̇/√3, μ = 60.6 GPa
    p = SetNonLinearPeierlsCreep(NonLinearPeierls.glaucophane_glide_Hufford_2026)
    @test p.Apparatus == SimpleShear
    ε = compute_εII(p, σ * 1.0e6 / 2; T) * 2 / sqrt(3)
    @test ε ≈ 10^38.0181 * (σ / 60.6e3)^2 * exp(-800.0e3 / (NumValue(p.R) * T) * (1 - σ / 22.0e3)^2) rtol = 1.0e-6

    # Above the Peierls stress the barrier vanishes (non-integer q = 1.18 must not throw)
    p = SetNonLinearPeierlsCreep(NonLinearPeierls.antigorite_Burdette_2022)
    σ_high = 3000.0
    @test compute_εII(p, σ_high * 1.0e6 / p.FT; T) * p.FE ≈ exp(-0.624) * (σ_high / 35.0e3)^2 rtol = 1.0e-6

    # Nondimensional Float32: the prefactor (≈ 5e43) exceeds floatmax(Float32); the rate is recomputed in Float64
    CharDim = GEO_units()
    p = SetNonLinearPeierlsCreep(NonLinearPeierls.glaucophane_glide_Hufford_2026, CharDim)
    τn, Tn = nondimensionalize(1.0e9Pa, CharDim), nondimensionalize(950.0K, CharDim)
    ε32 = compute_εII(p, Float32(τn); T = Float32(Tn))
    @test ε32 isa Float32
    @test ε32 ≈ compute_εII(p, τn; T = Tn) rtol = 1.0e-4
end
