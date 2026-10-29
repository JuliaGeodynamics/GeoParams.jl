using Test
using GeoParams
using LaTeXStrings
import GeoParams.Dislocation

@testset "DislocationCreepLaws" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 1000km)

    # Define a linear viscous creep law ---------------------------------
    x1 = DislocationCreep()
    @test isbits(x1)
    @test param_info(x1).Equation === L"$\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}$"
    str = sprint(show, x1)
    @test occursin("DislocationCreep:", str)
    @test x1.n.val == 1.0
    @test x1.A.val == 1.5e6

    x2 = DislocationCreep(; n = 3)
    @test x2.A.val == 1.5e6

    # perform a computation with the dislocation creep laws
    # Calculate EpsII, using a set of pre-defined values
    CharDim = GEO_units()
    EpsII = GeoUnit(0s^-1)
    EpsII = nondimensionalize(EpsII, CharDim)
    TauII = GeoUnit(100MPa)
    TauII = nondimensionalize(TauII, CharDim)
    P = GeoUnit(100MPa)
    P = nondimensionalize(P, CharDim)
    T = GeoUnit(500C)
    T = nondimensionalize(T, CharDim)
    f = GeoUnit(50MPa)
    f = nondimensionalize(f, CharDim)
    args = (; P = 100.0, T = 500.0, f = 50.0e6)
    Phase = SetMaterialParams(;
        Name = "Viscous Matrix",
        Phase = 2,
        Density = ConstantDensity(),
        CreepLaws = DislocationCreep(; n = 3NoUnits, r = 1NoUnits),
        CharDim = CharDim,
    )
    TauII = 1.0e6
    ε = compute_εII(x1, TauII, args)

    # Test some of the preset rheologies
    p = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003; n = 3.1)
    @test p.n.val == 3.1

    # Test some of the preset rheologies
    p = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
    TauII = 0.3e6Pa
    args = (; T = 1673.0K, P = 0.0Pa)
    ε = compute_εII(p, TauII, args)
    @test ε ≈ 2.7319532582144474e-13 / s

    # same but while removing the tensor correction
    ε_notensor = compute_εII(remove_tensor_correction(p), TauII, args)
    @test ε_notensor ≈ 4.612967949163285e-14 / s

    # test with arrays
    τII_array = ones(10) * 1.0e9
    ε_array = similar(τII_array)
    T_array = ones(size(τII_array)) * (500 + 273)

    args_array = (; T = T_array, P = 100.0, f = 5.0e7)

    compute_εII!(ε_array, x1, τII_array, args_array)
    @test ε_array[1] ≈ 1.5408816790109496e-17

    # EXERCISE 6.1 of the Gerya textbook (as implemented in the matlab exercise)
    Ad = 2.5e-17Pa^-3.5 * s^-1
    n = 3.5NoUnits
    Ea = 532.0e3J / mol
    V = 0m^3 / mol
    R = 8.3145J / K / mol
    T = (1200 + 273.15)K
    ε = 1.0e-14 / s
    F2 = 1 / 2^((n - 1) / n) / 3^((n + 1) / 2 / n)
    η_book = F2 / Ad^(1 / n) / ε^((n - 1) / n) * exp(Ea / n / R / T)
    τ_book = 2 * η_book * ε

    # Same but using GeoParams & dimensional values:
    p = SetDislocationCreep(Dislocation.dry_olivine_Gerya_2019)
    args = (; T = T)
    τ = compute_τII(p, ε, args)        # compute stress
    @test τ ≈ τ_book

    ε1 = compute_εII(p, τ, args)
    @test ε1 ≈ ε

    # Do the same, but using Floats as input
    args1 = (; T = ustrip(T))
    τ1 = compute_τII(p, 1.0e-14, args1)    # only using Floats
    @test ustrip(τ) ≈ τ1

    # compute devatoric stress & viscosity profiles (as in book)
    Tvec = (400 + 273.15):10:(1200 + 273.15)
    εvec = ones(size(Tvec)) * 1.0e-14
    τvec = zero(εvec)
    args2 = (; T = Tvec)
    compute_τII!(τvec, p, εvec, args2)    # only using Floats
    ηvec = τvec ./ (2 * εvec)
    @test sum(ηvec) / length(ηvec) ≈ 4.124658696991946e24

    p = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)

    args = (; T = (650 + 273.15))
    ε_vec = exp10.(-22:-12)
    τ_vec = zero(ε_vec)
    compute_τII!(τ_vec, p, ε_vec, args)
    η_vec = τ_vec ./ (2 * ε_vec)

    # Do some basic checks on all creeplaws in the DB
    CharDim = GEO_units()
    creeplaw_list = dislocation_law_list()
    for fun in creeplaw_list
        # @show key
        p = SetDislocationCreep(fun)        # original creep law
        p_nd = nondimensionalize(p, CharDim)    # non-dimensionalized
        @test isfinite(NumValue(p_nd.A))         # extreme n can over- or underflow the scaling
        p_dim = dimensionalize(p, CharDim)       # dimensionalized

        # Check that values are the same after non-dimensionalisation & dimensionalisation
        for field in fieldnames(typeof(p_dim))
            val_original = getfield(p, field)
            val_final = getfield(p_dim, field)
            if isa(val_original, GeoUnit)
                @test Value(val_original) == Value(val_final)
            end
        end

        # Perform computations with the rheology
        args = (T = 900.0, d = 100.0e-6, τII_old = 1.0e6)
        ε = 1.0e-15
        τ = compute_τII(p, ε, args)
        ε_test = compute_εII(p, τ, args)
        @test ε ≈ ε_test

        # the same law in Float32: in SI units the prefactor reaches 1e-55 and
        # τ^n reaches 1e39, neither of which Float32 holds, while the stresses
        # and strain rates they relate do
        τ32 = compute_τII(p, Float32(ε), args)
        @test Float64(τ32) ≈ τ rtol = 1.0e-3
        ε32 = compute_εII(p, Float32(τ), args)
        @test Float64(ε32) ≈ ε rtol = 1.0e-3

        # test overriding the default values
        # a = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)"; V=1e-6m^3 / mol)
        # @test Value(a.V) == 1e-6m^3 / mol
    end

    # ---- param_info name->function resolution (find_creep_law) ----
    # default (empty-name) law takes the early-return branch -> Equation only, no Comment/BibTex
    info_default = param_info(DislocationCreep())
    @test info_default isa MaterialParamsInfo
    @test isempty(info_default.BibTex_Reference)

    # named law: param_info resolves the builder from the stored Name and pulls the
    # database Comment + BibTex. (Previously threw `String not callable`.)
    v_named = SetDislocationCreep(Dislocation.dry_anorthite_Rybacki_2006)
    info_named = param_info(v_named)
    @test info_named isa MaterialParamsInfo
    @test !isempty(info_named.BibTex_Reference)
    # the resolved info must match calling the builder's database_info directly
    @test info_named.BibTex_Reference ==
        GeoParams.MaterialParameters.ConstitutiveRelationships.dislocation_database_info(Dislocation.dry_anorthite_Rybacki_2006).BibTex_Reference

    # name not present in the database -> graceful fallback to Equation only (isnothing(f) branch)
    v_unknown = DislocationCreep(; Name = "totally_made_up_law_xyz")
    info_unknown = param_info(v_unknown)
    @test info_unknown isa MaterialParamsInfo
    @test isempty(info_unknown.BibTex_Reference)

    # dimensional compute path: Quantity strain rate -> @unpack_units branch
    τ_dim = compute_τII(v_named, 1.0e-15 / s, (; T = 1273.15K, P = 1.0e9Pa))
    @test τ_dim isa Quantity
end

# Lusk et al. (2021), Table 2: ε = A σⁿ f^r exp(-(Q + PV)/RT) with σ, f, P in MPa;
# GeoParams takes the second invariant, so σ and ε go through FT and FE.
@testset "Lusk et al. (2021) parameters" begin
    σ, T, P, fH2O = 100.0, 623.15, 300.0, 42.0
    for (law, logA, n, r, Q, V) in (
            (Dislocation.lowP_wet_quartzite_Lusk_2021, -9.3, 3.5, 0.49, 118.0e3, 2.59),
            (Dislocation.wet_quartzite_Lusk_2021, -6.36, 2.1, 0.2, 94.0e3, 1.44),
            (Dislocation.highP_wet_quartzite_Lusk_2021, -7.9, 2.0, 0.49, 77.0e3, 2.59),
        )
        p = SetDislocationCreep(law)
        ε = compute_εII(p, σ * 1.0e6 / p.FT; T, P = P * 1.0e6, f = fH2O * 1.0e6) * p.FE
        @test ε ≈ 10^logA * σ^n * fH2O^r * exp(-(Q + V * P) / (NumValue(p.R) * T)) rtol = 1.0e-6
    end
end

# Corrected database entries, each against its flow law written in the paper's units:
# ε = A σⁿ f^r exp(-(Q + PV)/RT) with σ (differential stress) and fH2O in MPa, V in cm³/mol, P in MPa.
# GeoParams takes τII, P and fH2O in Pa; σ and ε go through FT and FE.
@testset "Dislocation creep database vs. paper units" begin
    σ, T, P, fH2O = 100.0, 1273.0, 300.0, 200.0
    for (law, A, n, r, Q, V) in (
            (Dislocation.wet_olivine_Hirth_2003, 1600.0, 3.5, 1.2, 520.0e3, 22.0),     # HK03 Table 1
            (Dislocation.wet_anorthite_Rybacki_2006, 10^0.2, 3.0, 1.0, 345.0e3, 38.0), # Rybacki 2006 Table 5
            (Dislocation.wet_clinopyroxene_Chen_2006, 10^6.7, 2.7, 3.0, 670.0e3, 0.0),
            (Dislocation.wet_quartzite_Hirth_2001, 10^-11.2, 4.0, 1.0, 135.0e3, 0.0),
            (Dislocation.wet_quartzite_Rutter_2004, 10^-4.93, 2.97, 1.0, 242.0e3, 0.0),
            (Dislocation.wet_quartzite_Tokle_2019, 1.1e-12, 3.0, 1.2, 115.0e3, 0.0),           # Tokle 2019 Table 1
            (Dislocation.lowstress_wet_quartzite_Tokle_2019, 1.75e-12, 4.0, 1.0, 125.0e3, 0.0),
            (Dislocation.wet_olivine_Mei_2000b, 1.5e3, 3.0, 0.98, 470.0e3, 20.0),               # Mei & Kohlstedt 2000b eq. (2)
            (Dislocation.dry_olivine_Karato_2003, 10^6.1, 3.0, 0.0, 510.0e3, 14.0),             # Karato & Jung 2003 Table 3
            (Dislocation.wet_olivine_Karato_2003, 10^2.9, 3.0, 1.2, 470.0e3, 24.0),
            (Dislocation.wet_quartzite_Lu_2019, 6.0e-15, 4.0, 2.7, 132.0e3, 35.3),              # Lu & Jiang 2019 eq. (6)
            (Dislocation.wet_quartzite_Lusk_2021, 10^-6.36, 2.1, 0.2, 94.0e3, 1.44),            # Lusk 2021 Table 2
            (Dislocation.lowP_wet_quartzite_Lusk_2021, 10^-9.3, 3.5, 0.49, 118.0e3, 2.59),
            (Dislocation.highP_wet_quartzite_Lusk_2021, 10^-7.9, 2.0, 0.49, 77.0e3, 2.59),
            (Dislocation.high_stress_wet_dunite_Chopra_1981, 10^3.98, 3.35, 0.0, 444.0e3, 0.0), # Chopra & Paterson 1981 Table IV
            (Dislocation.low_stress_wet_dunite_Chopra_1981, 10^3.88, 2.44, 0.0, 386.0e3, 0.0),
            (Dislocation.dry_dunite_Chopra_1984, 10^4.46, 3.6, 0.0, 535.0e3, 0.0),              # Chopra & Paterson 1984 Table 3
            (Dislocation.diabase_Caristan_1982, 6.12e-2, 3.05, 0.0, 276.0e3, 0.0),              # Caristan 1982 p. 6789
            (Dislocation.dry_clinopyroxene_Bystricky_Mackwell_2001, 10^9.8, 4.7, 0.0, 760.0e3, 0.0), # eq. (2)
            (Dislocation.strong_diabase_Mackwell_1998, 8.0, 4.7, 0.0, 485.0e3, 0.0),            # Mackwell 1998 p. 980
            (Dislocation.mafic_granulite_Wilks_1990, 1.4e4, 4.2, 0.0, 445.0e3, 0.0),            # Wilks & Carter 1990 Fig. 6B
            (Dislocation.plagioclase_An75_Ji_1993, 3.27e-4, 3.2, 0.0, 238.0e3, 0.0),            # Ji & Zhao 1993 Table 2
            (Dislocation.wet_omphacite_Zhang_2006, 10^-2.0, 3.5, 0.0, 310.0e3, 0.0),            # Zhang 2006 eq. (4)
            (Dislocation.wet_jadeit_Orzol_2006, exp(-3.3), 3.7, 0.0, 326.0e3, 0.0),             # Orzol 2006: ln A = -3.3
            (Dislocation.wet_olivine_Afonso_2004, 2.0e3, 4.0, 0.0, 471.0e3, 0.0),               # Afonso & Ranalli 2004 Table 1
            (Dislocation.rock_salt_Li_Urai_2016, 7.26e-6, 5.0, 0.0, 53.92e3, 0.0),              # Li & Urai 2016 Table 1
        )
        p = SetDislocationCreep(law)
        ε = compute_εII(p, σ * 1.0e6 / p.FT; T, P = P * 1.0e6, f = fH2O * 1.0e6) * p.FE
        @test ε ≈ A * σ^n * fH2O^r * exp(-(Q + V * P) / (NumValue(p.R) * T)) rtol = 1.0e-6
    end
    @test Dislocation.dry_dunite_Chopra_1984 in dislocation_law_list()

    # `_refit` entries: the printed law with A refitted to the paper's own data,
    # i.e. the printed law times the ratio of the prefactors
    for (refit, printed, ratio) in (
            (Dislocation.dry_olivine_Hirth_2003_refit, Dislocation.dry_olivine_Hirth_2003, 2.8e5 / 1.1e5),
            (Dislocation.wet_olivine_Karato_2003_refit, Dislocation.wet_olivine_Karato_2003, 10^-0.2),
            (Dislocation.wet_quartzite_Lu_2019_refit, Dislocation.wet_quartzite_Lu_2019, 1 / 6),
        )
        p, p0 = SetDislocationCreep(refit), SetDislocationCreep(printed)
        args = (; T, P = P * 1.0e6, f = fH2O * 1.0e6)
        @test compute_εII(p, σ * 1.0e6, args) ≈ ratio * compute_εII(p0, σ * 1.0e6, args) rtol = 1.0e-6
        @test refit in dislocation_law_list()
    end

    # Kirby (1983): A = 10^3.2 GPa^-2.4 s^-1, σ in GPa
    p = SetDislocationCreep(Dislocation.wet_quartzite_Kirby_1983)
    ε = compute_εII(p, σ * 1.0e6 / p.FT; T) * p.FE
    @test ε ≈ 10^3.2 * (σ / 1.0e3)^2.4 * exp(-160.0e3 / (NumValue(p.R) * T)) rtol = 1.0e-6

    # Jaoul et al. (1984), Table 1, rows 1 and 2: A in s^-1 kbar^-n, σ in kbar, E in kcal/mol
    for (law, A, n, E) in (
            (Dislocation.dry_quartzite_Jaoul_1984, 1.37, 2.8, 44.0),
            (Dislocation.wet_quartzite_Jaoul_1984, 0.63, 2.4, 39.0),
        )
        p = SetDislocationCreep(law)
        ε = compute_εII(p, σ * 1.0e6 / p.FT; T) * p.FE
        @test ε ≈ A * (σ / 100)^n * exp(-E * 4184 / (NumValue(p.R) * T)) rtol = 1.0e-6
    end

    # Dimanov & Dresen (2005), Tables 3a/3b: A in Pa^-n s^-1, σ in Pa
    for (law, A, n, Q) in (
            (Dislocation.dry_diopside_Dimanov_2005, 3.01e-28, 5.47, 691.0e3),
            (Dislocation.wet_diopside_Dimanov_2005, 5.16e-33, 5.52, 534.0e3),
        )
        p = SetDislocationCreep(law)
        ε = compute_εII(p, σ * 1.0e6 / p.FT; T) * p.FE
        @test ε ≈ A * (σ * 1.0e6)^n * exp(-Q / (NumValue(p.R) * T)) rtol = 1.0e-6
    end

    # HK03 Table 1, constant C_OH: f is C_OH in H/10^6 Si (dimensionless). Footnote e:
    # 90 (0.3)^3.5 (1000)^1.2 exp(-(480000 + 1e9 * 11e-6)/(8.314 * 1673)) = 2.5e-12 s^-1
    p = SetDislocationCreep(Dislocation.wet_olivine1_Hirth_2003)
    ε = compute_εII(p, 0.3e6 / p.FT; T = 1673.0, P = 1.0e9, f = 1000.0) * p.FE
    @test ε ≈ 90 * 0.3^3.5 * 1000^1.2 * exp(-(480.0e3 + 1.0e9 * 11.0e-6) / (NumValue(p.R) * 1673.0)) rtol = 1.0e-6
    @test ε ≈ 2.5e-12 rtol = 2.0e-2
end

# Each entry evaluated through GeoParams reproduces the printed law of its paper:
# σ (MPa) is the differential stress for AxialCompression; for SimpleShear the papers use
# σ = 2τ and ε̇ = γ̇/√3, so the GeoParams shear strain rate is divided by √3.
@testset "Dislocation creep database vs. printed laws (crust, evaporites, ice)" begin
    σ, T, P = 100.0, 1000.0, 1.0e8
    for (law, A, n, Q, V) in (
            (Dislocation.wet_quartzite_Gleason_1995, 1.1e-4, 4.0, 223.0e3, 0.0),
            (Dislocation.melt_quartzite_Gleason_1995, 1.8e-8, 4.0, 137.0e3, 0.0),
            (Dislocation.wet_quartzite_Gleason_1995_Holyoke_2010, 5.1e-4, 4.0, 223.0e3, 0.0),
            (Dislocation.wet_quartzite_Paterson_1990, 6.5e-8, 3.1, 135.0e3, 0.0),
            (Dislocation.dry_quartzite_Brace_1980, 5.0e-6, 3.0, 190.0e3, 0.0),
            (Dislocation.felsic_granulite_Wilks_1990, 8.0e-3, 3.1, 243.0e3, 0.0),
            (Dislocation.felsic_granulite_Wen_2021, 10^-4.4, 4.2, 260.0e3, 0.0),
            (Dislocation.weak_diabase_Mackwell_1998, 190.0, 4.7, 485.0e3, 0.0),
            (Dislocation.wet_plagioclase_An60_Rybacki_2004, 10^-1.5, 3.0, 235.0e3, 0.0),
            (Dislocation.amphibolite_Hacker_1990, 1.1e-4, 3.7, 244.0e3, 0.0),
            (Dislocation.dry_eclogite_Zhang_2007, 10^3.3, 3.5, 403.0e3, 27.2e-6),
            (Dislocation.calcite_coarse_dislocation_Gilgannon_2025, 2.16e-8 * 0.001^-0.6, 7.0, 200.0e3, 0.0),
            (Dislocation.rock_salt_Carter_1993, 1.6e-4, 5.3, 68.1e3, 0.0),
            (Dislocation.lowstress_rock_salt_Carter_1993, 8.1e-5, 3.4, 51.6e3, 0.0),
            (Dislocation.anhydrite_Dorner_2014, 4.8e-19 * 1.0e6^3.9, 3.9, 338.0e3, 0.0),
            (Dislocation.ice_dislocation_Goldsby_2001, 6.0e28, 4.0, 181.0e3, -13.0e-6),
            (Dislocation.ice_dislocation_cold_Goldsby_2001, 4.0e5, 4.0, 60.0e3, -13.0e-6),
            (Dislocation.ice_basal_slip_Goldsby_2001, 5.5e7, 2.4, 60.0e3, -13.0e-6),
            (Dislocation.ice_regime_A_Durham_2001, 10^11.8, 4.0, 91.0e3, -13.0e-6),
            (Dislocation.ice_regime_B_Durham_2001, 10^5.1, 4.0, 61.0e3, -13.0e-6),
            (Dislocation.ice_II_cold_Durham_2001, 10^1.84, 5.3, 55.0e3, 7.0e-6),
            (Dislocation.ice_II_warm_Durham_2001, 10^11.7, 5.2, 98.0e3, 4.0e-6),
            (Dislocation.ice_III_warm_Durham_2001, 10^26.4, 4.3, 151.0e3, 0.0),
            (Dislocation.ice_V_Durham_2001, 10^23.0, 6.0, 136.0e3, 29.0e-6),
            (Dislocation.ice_VI_Durham_2001, 10^6.7, 4.5, 66.0e3, 11.0e-6),
            (Dislocation.co2_ice_Durham_1999, 10^3.86, 5.6, 33.0e3, 0.0),
        )
        p = SetDislocationCreep(law)
        @test p.Apparatus == AxialCompression
        ε = compute_εII(p, σ * 1.0e6 / p.FT; T, P) * p.FE
        @test ε ≈ A * σ^n * exp(-(Q + P * V) / (NumValue(p.R) * T)) rtol = 1.0e-6
    end

    # Keefner et al. (2011): oxygen-fugacity term, f = fO2 in Pa
    p = SetDislocationCreep(Dislocation.dry_dunite_fO2_Keefner_2011)
    fO2 = 1.0e-3
    ε = compute_εII(p, σ * 1.0e6 / p.FT; T, P, f = fO2) * p.FE
    @test ε ≈ 10^2.6 * σ^3.6 * fO2^0.2 * exp(-449.0e3 / (NumValue(p.R) * T)) rtol = 1.0e-6

    # Hufford et al. (2026): σ = 2τ, ε̇ = γ̇/√3
    p = SetDislocationCreep(Dislocation.glaucophane_Hufford_2026)
    @test p.Apparatus == SimpleShear
    ε = compute_εII(p, σ * 1.0e6 / 2; T) * 2 / sqrt(3)
    @test ε ≈ 2.32e10 * σ^3 * exp(-450.0e3 / (NumValue(p.R) * T)) rtol = 1.0e-6

    # `_refit` entries are the printed law times the ratio of the prefactors
    for (refit, printed, ratio) in (
            (Dislocation.felsic_granulite_Wen_2021_refit, Dislocation.felsic_granulite_Wen_2021, 10^-0.3),
            (Dislocation.ice_dislocation_cold_Goldsby_2001_refit, Dislocation.ice_dislocation_cold_Goldsby_2001, 5.4e4 / 4.0e5),
        )
        args = (; T, P)
        @test compute_εII(SetDislocationCreep(refit), 1.0e7, args) ≈ ratio * compute_εII(SetDislocationCreep(printed), 1.0e7, args) rtol = 1.0e-6
        @test refit in dislocation_law_list()
    end
end
