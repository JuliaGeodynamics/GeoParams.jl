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

# perform a computation with the dislocation creep laws 
    # Calculate EpsII, using a set of pre-defined values
CharDim = GEO_units()
EpsII  = GeoUnit(0s^-1)
        Nondimensionalize!(EpsII,CharDim)
TauII   = GeoUnit(100MPa)
        Nondimensionalize!(TauII,CharDim)
P       = GeoUnit(100MPa)
        Nondimensionalize!(P,CharDim)
T       = GeoUnit(500C)
        Nondimensionalize!(T,CharDim)
f       = GeoUnit(50MPa)
        Nondimensionalize!(f,CharDim)
p       = CreepLawVariables(P=P,T=T,f=f)
Phase   = SetMaterialParams(Name="Viscous Matrix", Phase=2,
                                     Density   = ConstantDensity(),
                                     CreepLaws = DislocationCreep(n=3NoUnits, r=1NoUnits), CharDim = CharDim)
EpsII.val = ComputeCreepLaw_EpsII(TauII,Phase.CreepLaws[1],p)
@test round(EpsII.val - 2.1263214994323903e-11, digits=14) == 0     
    # Check that once inverted, we get back the TauII that we used to calculate EpsII
NewTau = ComputeCreepLaw_TauII(EpsII,Phase.CreepLaws[1],p)
@test round(NewTau - TauII.val,digits=14) == 0


# Given stress
#@test ComputeCreepLaw_EpsII(1e6Pa, x1, CreepLawParams())==5e-13/s                # dimensional input       

# Given strainrate 
#@test ComputeCreepLaw_TauII(1e-13/s, x1, CreepLawParams())==1e18*2*1e-13Pa       # dimensional input       
# -------------------------------------------------------------------


end