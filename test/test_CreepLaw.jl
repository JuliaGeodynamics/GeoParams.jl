using Test
using GeoParams

@testset "CreepLaw" begin

    #Make sure structs are isbits
    x = LinearViscous()
    @test isbits(x)
    @test isvolumetric(x) == false

    x = PowerlawViscous()
    @test isbits(x)

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=1000km)

    # Define a linear viscous creep law ---------------------------------
    x1 = LinearViscous(; η=1e18Pa * s)
    @test x1.η.val == 1e18

    x1_ND = LinearViscous(; η=1e18Pa * s)
    @test isDimensional(x1_ND) == true
    x1_ND = nondimensionalize(x1_ND, CharUnits_GEO)                 # check that we can nondimensionalize all entries within the struct
    @test isDimensional(x1_ND) == false
    @test x1_ND.η * 1.0 == 0.1
    x1_ND = dimensionalize(x1_ND, CharUnits_GEO)                    # check that we can dimensionalize it again
    @test x1_ND.η.val == 1e18
    x1_ND = nondimensionalize(x1_ND, CharUnits_GEO)

    # Given stress
    args = (;)
    @test compute_εII(x1, 1e6Pa, args) == 5e-13 / s                # dimensional input       
    @test compute_εII(x1_ND, 1e0, args) == 5.0                   # non-dimensional

    τ = [0.0; 0.0]
    compute_εII!(τ, x1_ND, [1e0; 2.0], args)
    @test τ == [5.0; 10.0]    # vector input

    # Given strainrate 
    @test compute_τII(x1, 1e-13 / s, args) == 1e18 * 2 * 1e-13Pa       # dimensional input       
    @test compute_τII(x1_ND, 1e0, args) == 0.2                   # non-dimensional
    ε = [0.0; 0.0]
    compute_τII!(ε, x1_ND, [1e0; 2.0], args)
    @test ε == [0.2; 0.4]     # vector input
    # -------------------------------------------------------------------

    # -------------------------------------------------------------------
    # Define powerlaw viscous rheology
    x2 = PowerlawViscous()
    x2 = nondimensionalize(x2, CharUnits_GEO)
    @test NumValue(x2.ε0) == 0.001                                 # powerlaw 

    # Given strainrate 
    compute_τII(x1, 1e-13 / s, args)
    @test compute_τII(x1, 1e-13 / s, args) == 1e18 * 2 * 1e-13Pa       # dimensional input       
    @test compute_τII(x1_ND, 1e0, args) == 0.2                   # non-dimensional
    @test [compute_τII(x1_ND, EII, args) for EII in (1e0, 2.0)] == [0.2; 0.4]     # vector input
    # -------------------------------------------------------------------

    # -------------------------------------------------------------------
    # Define powerlaw viscous rheology
    x2 = PowerlawViscous()
    x2 = nondimensionalize(x2, CharUnits_GEO)
    @test x2.ε0.val == 0.001                                 # powerlaw 
    # -------------------------------------------------------------------

    # ArrheniusType rheology --------------------------------------------
    args = (;)
    x3 = ArrheniusType()
    # For a given temperature
    T = 1.0
    # and stress 
    τ1 = 1.0
    #eps     =   compute_εII(x3,τ1;T,args)
    #@test   eps ≈ 0.5

    args = (T=T,)
    eps = compute_εII(x3, τ1, args)
    @test eps ≈ 0.5

    # and strain rate
    ε1 = 0.5
    tau = compute_τII(x3, ε1, args)
    @test tau ≈ 1.0
    # using a vector input ------------------    
    T2 = [1.0; 0.9; 0.8]
    # and stress 
    ε21 = [0.0; 0.0; 0.0]
    τ21 = [1.0; 0.8; 0.9]
    args = (T=T2,)
    compute_εII!(ε21, x3, τ21, args)
    # and strain rate
    τ22 = [0.0; 0.0; 0.0]
    ε22 = [0.5; 1.0; 0.2]
    compute_τII!(τ22, x3, ε22, args)
    # derivatives
    sol1 = dεII_dτII(x3, τ22)
    @test sol1 ≈ 0.5
    sol2 = dτII_dεII(x3, ε21)
    @test sol2 ≈ 2.0
end
