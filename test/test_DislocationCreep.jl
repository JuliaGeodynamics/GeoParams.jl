using Test
using GeoParams

@testset "DislocationCreepLaws" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
                
# Define a linear viscous creep law ---------------------------------
x1      =   DislocationCreep()
@test x1.n.val == 1.0
@test x1.A.val == 1.5MPa^-1*s^-1

x2      =   DislocationCreep(n=3)
@test x2.A.val == 1.5MPa^-3*s^-1

# perform a computation with the viscous creep laws 

# Given stress
#@test ComputeCreepLaw_EpsII(1e6Pa, x1, CreepLawParams())==5e-13/s                # dimensional input       

# Given strainrate 
#@test ComputeCreepLaw_TauII(1e-13/s, x1, CreepLawParams())==1e18*2*1e-13Pa       # dimensional input       
# -------------------------------------------------------------------


end