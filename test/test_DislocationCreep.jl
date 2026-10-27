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
            (Dislocation.wet_quartzite_Tokle_2019, 10^-11.959, 3.0, 1.2, 115.0e3, 0.0),
            (Dislocation.lowstress_wet_quartzite_Tokle_2019, 1.75e-12, 4.0, 1.0, 125.0e3, 0.0),
            (Dislocation.dry_quartzite_Jaoul_1984, 10^-5.415, 2.8, 0.0, 184.0e3, 0.0),
            (Dislocation.wet_quartzite_Jaoul_1984, 10^-5.045, 2.8, 0.0, 163.0e3, 0.0),
        )
        p = SetDislocationCreep(law)
        ε = compute_εII(p, σ * 1.0e6 / p.FT; T, P = P * 1.0e6, f = fH2O * 1.0e6) * p.FE
        @test ε ≈ A * σ^n * fH2O^r * exp(-(Q + V * P) / (NumValue(p.R) * T)) rtol = 1.0e-6
    end

    # Kirby (1983): A = 10^3.2 GPa^-2.4 s^-1, σ in GPa
    p = SetDislocationCreep(Dislocation.wet_quartzite_Kirby_1983)
    ε = compute_εII(p, σ * 1.0e6 / p.FT; T) * p.FE
    @test ε ≈ 10^3.2 * (σ / 1.0e3)^2.4 * exp(-160.0e3 / (NumValue(p.R) * T)) rtol = 1.0e-6

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
