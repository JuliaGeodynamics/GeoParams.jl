using Test
using GeoParams

@testset "CreepLaw" begin

    #Make sure structs are isbits
    x = LinearViscous()
    @test isbits(x)

    x = PowerlawViscous()
    @test isbits(x)

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=1000km)

    # Define a linear viscous creep law ---------------------------------
    x1 = LinearViscous(; η=1e18Pa * s)
    @test x1.η.val == 1e18

# Given stress
args = (;)
@test compute_εII(x1,1e6Pa, args)==5e-13/s                # dimensional input       
@test compute_εII(x1_ND,1e0, args)==5.0                   # non-dimensional

τ = [0.0; 0.0]
compute_εII!(τ,  x1_ND, [1e0; 2.0], args)
@test τ==[5.0; 10.0]    # vector input

# Given strainrate 
@test compute_τII(x1,1e-13/s, args)==1e18*2*1e-13Pa       # dimensional input       
@test compute_τII(x1_ND, 1e0, args)==0.2                   # non-dimensional
ε = [0.0; 0.0]
compute_τII!(ε,x1_ND, [1e0; 2.0], args)
@test ε==[0.2; 0.4]     # vector input
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Define powerlaw viscous rheology
x2  =   PowerlawViscous()
x2 = nondimensionalize(x2,CharUnits_GEO)
@test NumValue(x2.ε0) ==0.001                                 # powerlaw 

    # Given strainrate 
    @test computeCreepLaw_TauII(1e-13 / s, x1, args) == 1e18 * 2 * 1e-13Pa       # dimensional input       
    @test computeCreepLaw_TauII(1e0, x1_ND, args) == 0.2                   # non-dimensional
    @test computeCreepLaw_TauII([1e0; 2.0], x1_ND, args) == [0.2; 0.4]     # vector input
    # -------------------------------------------------------------------

    # -------------------------------------------------------------------
    # Define powerlaw viscous rheology
    x2 = PowerlawViscous()
    x2 = nondimensionalize(x2, CharUnits_GEO)
    @test x2.ε0.val == 0.001                                 # powerlaw 

    # perform computations

    # -------------------------------------------------------------------

    #include("test_DislocationCreep.jl")

end
