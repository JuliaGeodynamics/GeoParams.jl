using Test
using GeoParams
import GeoParams.NonLinearPeierls

@testset "NonLinearPeierlsCreepLaws" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1.0Pa * s, length=1.0m)

    # Define a linear viscous creep law ---------------------------------
    x1 = NonLinearPeierlsCreep()
    @test isbits(x1)
    @test Value(x1.n) == 2.0
    @test Value(x1.q) == 1.0
    @test Value(x1.o) == 0.5
    @test Value(x1.A) == 5.7e11MPa^(-2.0) * s^(-1.0)

    # perform a computation with the peierls creep laws
    # Calculate EpsII, using a set of pre-defined values
    CharDim = GEO_units(;
        length=1000.0km, viscosity=1.0e19Pa * s, stress=100.0MPa, temperature=1000.0C
    )
    EpsII = GeoUnit(1.0s^-1.0)
    TauII = GeoUnit(1.0e3MPa)
    T     = GeoUnit(473.0C)

    # compute a pure non linear peierls creep rheology
    p     = SetNonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010)
    T     = 200.0 + 273.15
    args  = (; T=T)
    TauII = 1.0e9
    ε     = compute_εII(p, TauII, args)
    @test ε ≈ 5.638821307626487e-17

    # same but with removing the tensor correction
    ε_notensor = compute_εII(remove_tensor_correction(p), TauII, args)
    @test ε_notensor ≈ 5.491207069105064e-22

    # test with arrays
    τII_array  = ones(10) * 1.0e9
    ε_array    = similar(τII_array)
    T_array    = ones(size(τII_array)) * (200.0 + 273.15)
    args_array = (; T=T_array)
    compute_εII!(ε_array, p, τII_array, args_array)
    @test ε_array[1] ≈ ε

    # compute when args are scalars
    compute_εII!(ε_array, p, τII_array, args)
    @test ε_array[1] ≈ ε

    # wet olivine, stress-strainrate curve
    # p    = SetNonLinearPeierlsCreep("Dry Olivine | Mei et al. (2010)")
    εII  = exp10.(-22:0.5:-12)
    τII  = exp10.(9:0.05:10)               # preallocate array
    T    = 200.0 + 273.15
    args = (; T=T)
    εII  = zero(τII)
    compute_εII!(εII, p, τII, args)
    eta_array1 = @. 0.5 * τII / εII

    # test overriding the default values
    # a = SetNonLinearPeierlsCreep("Dry Olivine | Mei et al. (2010)"; E=475.0kJ / mol)
    # @test Value(a.E) == 475.0kJ / mol

    # Do some basic checks on all creeplaws in the DB
    CharDim = GEO_units()
    creeplaw_list = nonlinearpeierls_law_list() # all creeplaws in database
    for fun in creeplaw_list
        p     = SetNonLinearPeierlsCreep(fun)                # original creep law
        p_nd  = nondimensionalize(p, CharDim)                # non-dimensionalized
        @test p_nd == SetNonLinearPeierlsCreep(fun, CharDim) # check that the non-dimensionalized version is the same as the original
        p_dim = dimensionalize(p, CharDim)                   # dimensionalized

        # Check that values are the same after non-dimensionalisation & dimensionalisation
        for field in fieldnames(typeof(p_dim))
            val_original = getfield(p, field)
            val_final    = getfield(p_dim, field)
            if isa(val_original, GeoUnit)
                @test Value(val_original) == Value(val_final)
            end
        end

        # Perform computations with the rheology
        args   = (T=273 + 200.0, d=100e-6, τII_old=1e6)
        ε      = 1.0e-15
        τ_ini  = 1.75e9                                      # initial guess for stress iterations
        τ      = Peierls_stress_iterations(p, τ_ini, ε, args)
        ε_test = compute_εII(p, τ; args...)
        @test ε ≈ ε_test
    end

    # check that derivative calculation is correct
    # p        = SetNonLinearPeierlsCreep("Dry Olivine | Mei et al. (2010)")
    TauII    = 1.0e9
    args     = (; T=273 + 200.0)
    depsdtau = dεII_dτII(p, TauII; args...)
    @test depsdtau ≈ 1.339830415314779e-24
end
