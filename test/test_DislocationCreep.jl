using Test
using GeoParams
import GeoParams.Dislocation

@testset "DislocationCreepLaws" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 1000km)

    # Define a linear viscous creep law ---------------------------------
    x1 = DislocationCreep()
    @test isbits(x1)
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

        # test overriding the default values
        # a = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)"; V=1e-6m^3 / mol)
        # @test Value(a.V) == 1e-6m^3 / mol
    end
end
