using Test
using GeoParams

@testset "CreepLaw" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
                
# Define a linear viscous creep law ---------------------------------
x1      =   LinearViscous(η=1e18Pa*s)
@test x1.η.val == 1e18Pa*s

x1_ND   =   LinearViscous(η=1e18Pa*s)
@test  isDimensional(x1_ND)==true 
Nondimensionalize!(x1_ND,CharUnits_GEO)                 # check that we can nondimensionalize all entries within the struct
@test  isDimensional(x1_ND)==false 
@test x1_ND.η*1.0==0.1                            
Dimensionalize!(x1_ND,CharUnits_GEO)                    # check that we can dimensionalize it again
@test x1_ND.η.val==1e18Pa*s
Nondimensionalize!(x1_ND,CharUnits_GEO)       

# perform a computation with the viscous creep laws 

# Given stress
@test ComputeCreepLaw_EpsII(1e6Pa, x1, CreepLawVariables())==5e-13/s                # dimensional input       
@test ComputeCreepLaw_EpsII(1e0, x1_ND, CreepLawVariables())==5.0                   # non-dimensional
@test ComputeCreepLaw_EpsII([1e0; 2.0], x1_ND, CreepLawVariables())==[5.0; 10.0]    # vector input

# Given strainrate 
@test ComputeCreepLaw_TauII(1e-13/s, x1, CreepLawVariables())==1e18*2*1e-13Pa       # dimensional input       
@test ComputeCreepLaw_TauII(1e0, x1_ND, CreepLawVariables())==0.2                   # non-dimensional
@test ComputeCreepLaw_TauII([1e0; 2.0], x1_ND, CreepLawVariables())==[0.2; 0.4]     # vector input
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Define powerlaw viscous rheology
x2  =   PowerlawViscous()
Nondimensionalize!(x2,CharUnits_GEO)
@test x2.ε0.val ==0.001                                 # powerlaw 

# perform computations


# -------------------------------------------------------------------


#include("test_DislocationCreep.jl")

end