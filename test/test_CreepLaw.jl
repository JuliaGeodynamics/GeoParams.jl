using Test
using GeoParams

@testset "CreepLaw" begin

#Make sure structs are isbits
x = LinearViscous()
@test isbits(x)
 
x = PowerlawViscous()
@test isbits(x)

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
                
# Define a linear viscous creep law ---------------------------------
x1      =   LinearViscous(η=1e18Pa*s)
@test x1.η.val == 1e18

x1_ND   =   LinearViscous(η=1e18Pa*s)
@test  isDimensional(x1_ND)==true 
x1_ND = nondimensionalize(x1_ND,CharUnits_GEO)                 # check that we can nondimensionalize all entries within the struct
@test  isDimensional(x1_ND)==false 
@test x1_ND.η*1.0==0.1                            
x1_ND = dimensionalize(x1_ND,CharUnits_GEO)                    # check that we can dimensionalize it again
@test x1_ND.η.val==1e18
x1_ND = nondimensionalize(x1_ND,CharUnits_GEO)       

# perform a computation with the viscous creep laws 

# Given stress
vars = CreepLawVariables() 
args = (P=vars.P, T=vars.T, f=vars.f, d=vars.d)
@test computeCreepLaw_EpsII(1e6Pa, x1, args) ==5e-13/s                # dimensional input       
@test computeCreepLaw_EpsII(1e0, x1_ND, args)==5.0                   # non-dimensional
@test computeCreepLaw_EpsII([1e0; 2.0], x1_ND, args)==[5.0; 10.0]    # vector input

# Given strainrate 
@test computeCreepLaw_TauII(1e-13/s, x1, args)==1e18*2*1e-13Pa       # dimensional input       
@test computeCreepLaw_TauII(1e0, x1_ND, args)==0.2                   # non-dimensional
@test computeCreepLaw_TauII([1e0; 2.0], x1_ND, args)==[0.2; 0.4]     # vector input
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Define powerlaw viscous rheology
x2  =   PowerlawViscous()
x2 = nondimensionalize(x2,CharUnits_GEO)
@test x2.ε0.val ==0.001                                 # powerlaw 

# perform computations


# -------------------------------------------------------------------


#include("test_DislocationCreep.jl")

end