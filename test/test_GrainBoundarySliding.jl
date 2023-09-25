using Test
using GeoParams

@testset "GrainBoundarySliding" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1Pa * s, length=1m)

    # Define a linear viscous creep law ---------------------------------
    x1 = GrainBoundarySliding()
    @test Value(x1.n) == 3.5
    @test Value(x1.p) == -2.0
    @test Value(x1.A) == 6500.0MPa^(-3.5) * s^(-1.0) * µm^(2.0)

    # perform a computation with the dislocation creep laws 
    # Calculate EpsII, using a set of pre-defined values
    CharDim = GEO_units(;
        length=1000km, viscosity=1e19Pa * s, stress=100MPa, temperature=1000C
    )
    EpsII = GeoUnit(1.0s^-1.0)
    TauII = GeoUnit(0.3MPa)
    P = GeoUnit(1.0e9Pa)
    T = GeoUnit(1400.0C)
    d = GeoUnit(10.0mm)
    d_nd = nondimensionalize(d, CharDim)

    # compute a pure diffusion creep rheology
    p = SetGrainBoundarySliding("Dry Olivine < 1523K | Hirth and Kohlstedt (2003)")

    T = 650.0 + 273.15

    args = (; T=T)
    TauII = 1.0e6
    ε = compute_εII(p, TauII, args)
    @test ε ≈ 8.968730687982234e-31

    # same but while removing the tensor correction
    ε_notensor = compute_εII(remove_tensor_correction(p), TauII, args)
    @test ε_notensor ≈ 1.5143914737172087e-31

    # test with arrays
    τII_array = ones(10) * 1.0e6
    ε_array = similar(τII_array)
    T_array = ones(size(τII_array)) * (650.0 + 273.15)

    args_array = (; T=T_array)

    compute_εII!(ε_array, p, τII_array, args_array)
    @test ε_array[1] ≈ ε

    # compute when args are scalars
    compute_εII!(ε_array, p, τII_array, args)
    @test ε_array[1] ≈ ε

    # ===

    # dry olivine, stress-strainrate curve
    p = SetGrainBoundarySliding("Dry Olivine < 1523K | Hirth and Kohlstedt (2003)")
    εII = exp10.(-22:0.5:-12)
    τII = zero(εII)                # preallocate array
    T = 650.0 + 273.15
    gsiz = 100.0e-6
    args = (T=T, d=gsiz)
    compute_τII!(τII, p, εII, args)

    eta_array = @. 0.5 * τII / εII

    εII = zero(τII)
    compute_εII!(εII, p, τII, args)
    eta_array1 = @. 0.5 * τII / εII

    # test overriding the default values
    a =  SetGrainBoundarySliding("Dry Olivine >= 1523K | Hirth and Kohlstedt (2003)", V=1e-6m^3/mol)
    @test Value(a.V) == 1.0e-6m^3/mol


    # Do some basic checks on all creeplaws in the DB
    CharDim = GEO_units()
    creeplaw_list = GrainBoundarySliding_info       # all creeplaws in database
    for (key, val) in creeplaw_list
        p     = SetGrainBoundarySliding(key)        # original creep law
        p_nd  = nondimensionalize(p,CharDim)    # non-dimensionalized
        p_dim = dimensionalize(p,CharDim)       # dimensionalized

        # Check that values are the same after non-dimensionalisation & dimensionalisation
        for field in fieldnames(typeof(p_dim))
            val_original = getfield(p,    field)
            val_final    = getfield(p_dim,field)
            if isa(val_original, GeoUnit)
                @test Value(val_original) == Value(val_final)        
            end
        end
        
        # Perform computations with the rheology
        args   = (T=900.0, d=100.0e-6)
        ε      = 1.0e-15
        τ      = compute_τII(p,ε,args)
        ε_test = compute_εII(p,τ,args)
        @test ε ≈ ε_test

    end

end
