using Test
using GeoParams
import GeoParams.Peierls

@testset "PeierlsCreep" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1Pa * s, length = 1m)

    # Define a linear viscous creep law ---------------------------------
    x1 = PeierlsCreep()
    @test isbits(x1)
    @test Value(x1.n) == 1.0
    @test Value(x1.q) == 2.0
    @test Value(x1.A) == 5.7e11s^(-1.0)
    @test repr("text/plain", x1) isa String

    # perform a computation with the dislocation creep laws
    # Calculate EpsII, using a set of pre-defined values
    CharDim = GEO_units(;
        length = 1000km, viscosity = 1.0e19Pa * s, stress = 100MPa, temperature = 1000C
    )
    EpsII = GeoUnit(1.0s^-1.0)
    TauII = GeoUnit(1.0e6MPa)
    T = GeoUnit(600C)

    # compute a pure diffusion creep rheology
    p = SetPeierlsCreep(Peierls.dry_olivine_Goetze_1979)
    T = 600 + 273.15
    args = (; T = T)
    TauII = 1.0e9
    ε = compute_εII(p, TauII, args)
    @test ε ≈ 9.127028135349583e-20

    # same but while removing the tensor correction
    ε_notensor = compute_εII(remove_tensor_correction(p), TauII, args)
    @test ε_notensor ≈ 1.3652144989670166e-20

    # test with arrays
    τII_array = ones(10) * 1.0e9
    ε_array = similar(τII_array)
    T_array = ones(size(τII_array)) * (600.0 + 273.15)
    args_array = (; T = T_array)

    compute_εII!(ε_array, p, τII_array, args_array)
    @test ε_array[1] ≈ ε

    # compute when args are scalars
    compute_εII!(ε_array, p, τII_array, args)
    @test ε_array[1] ≈ ε
    # ===

    # test overriding the default values
    # a = SetPeierlsCreep("Dry Olivine | Goetze and Evans (1979)"; E=535.0kJ / mol)
    # @test Value(a.E) == 535.0kJ / mol

    # Do some basic checks on all creeplaws in the DB
    CharDim = GEO_units()
    creeplaw_list = peierls_law_list()       # all creeplaws in database
    for fun in creeplaw_list
        p = SetPeierlsCreep(fun)                # original creep law
        p_nd = nondimensionalize(p, CharDim)        # non-dimensionalized
        @test p_nd == SetPeierlsCreep(fun, CharDim) # check that the non-dimensionalized version is the same as the original
        p_dim = dimensionalize(p, CharDim)           # dimensionalized

        # Check that values are the same after non-dimensionalisation & dimensionalisation
        for field in fieldnames(typeof(p_dim))
            val_original = getfield(p, field)
            val_final = getfield(p_dim, field)
            if isa(val_original, GeoUnit)
                @test Value(val_original) == Value(val_final)
            end
        end

        # Perform computations with the rheology
        args = (T = 500.0, τII_old = 2.2e9)
        ε = 1.0e-7
        τ = compute_τII(p, ε, args)
        ε_test = compute_εII(p, τ, args)
        @test ε ≈ ε_test
    end
end
