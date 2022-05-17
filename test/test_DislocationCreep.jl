using Test
using GeoParams

@testset "DislocationCreepLaws" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=1000km)

    # Define a linear viscous creep law ---------------------------------
    x1 = DislocationCreep()
    @test isbits(x1)
    @test x1.n.val == 1.0
    @test x1.A.val == 1.5

    x2 = DislocationCreep(; n=3)
    @test x2.A.val == 1.5

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
    args = (; P=100., T=500., f=50e6)
    Phase = SetMaterialParams(;
        Name="Viscous Matrix",
        Phase=2,
        Density=ConstantDensity(),
        CreepLaws=DislocationCreep(; n=3NoUnits, r=1NoUnits),
        CharDim=CharDim,
    )
    TauII = 1e6
    ε = compute_εII(x1, TauII, args)



    # test with arrays
    τII_array       =   ones(10)*1e6
    ε_array         =   similar(τII_array)
    T_array         =   ones(size(τII_array))*(500)

    args_array = (;T=T_array, P=100., f=5e7 )

    compute_εII!(ε_array, x1, τII_array, args_array)
    @test ε_array[1] ≈ ε 

#    εII = computeCreepLaw_EpsII(TauII, Phase.CreepLaws[1], p)
#    @test εII ≈ 2.1263214994323903e-11 rtol = 1e-8

 #   @test εII == computeCreepLaw_EpsII(TauII, Phase.CreepLaws[1], P.val, T.val, f.val)

    # Check that once inverted, we get back the TauII that we used to calculate EpsII
 #   NewTau = computeCreepLaw_TauII(εII, Phase.CreepLaws[1], p)
 #   @test NewTau ≈ TauII.val
 #   @test NewTau == computeCreepLaw_TauII(εII, Phase.CreepLaws[1], P.val, T.val, f.val)

  

    # Given strainrate 
    #@test computeCreepLaw_EpsII(1e-13/s, x1, CreepLawParams())==1e18*2*1e-13Pa       # dimensional input       
    # -------------------------------------------------------------------

end
