using Test
using GeoParams

@testset "DiffusionCreepLaws" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1Pa*s, length=1m);
                
# Define a linear viscous creep law ---------------------------------
x1      =   DiffusionCreep()
@test x1.n.val == 1.0
@test x1.p.val == -3.0
@test x1.A.val == 1.5MPa^-1.0*s^-1*m^(3.0)

x2      =   DiffusionCreep(n=3, p=-3.0)
@test x2.A.val == 1.5MPa^-3.0*s^-1*m^(3.0)

# perform a computation with the dislocation creep laws 
    # Calculate EpsII, using a set of pre-defined values
CharDim = GEO_units(length=1000km, viscosity=1e19Pa*s, stress=100MPa, temperature=1000C)
EpsII   = GeoUnit(1.0s^-1.0)
        Nondimensionalize!(EpsII,CharDim)
TauII   = GeoUnit(0.3MPa)
        Nondimensionalize!(TauII,CharDim)
P       = GeoUnit(1.0e9Pa)
        Nondimensionalize!(P,CharDim)
T       = GeoUnit(1400C)
        Nondimensionalize!(T,CharDim)
f       = GeoUnit(1000NoUnits)
        Nondimensionalize!(f,CharDim)
d       = GeoUnit(10mm)
        Nondimensionalize!(d,CharDim)
c       = CreepLawVariables(P=P,T=T,f=f,d=d)
Phase   = SetMaterialParams(Name="Viscous Matrix", Phase=1,
                                     Density   = ConstantDensity(),
                                     CreepLaws = DiffusionCreep(n=1.0NoUnits, r=1.0NoUnits, p=-3.0NoUnits, A=1.0e6MPa^(-1.0)*s^(-1.0)*µm^(3.0), E=335000.0J/mol, V=4.0e-6m^(3.0)/mol, Apparatus="Invariant"), CharDim = CharDim)
EpsII.val = ComputeCreepLaw_EpsII(TauII,Phase.CreepLaws[1],c)
# Tested value is Dimensional so need to non-dimensionalize it back to compare
@test EpsII.val  ≈ Nondimensionalize(7.823165072337786e-15s^(-1.0), CharDim) rtol = 1e-12
    # Check that once inverted, we get back the TauII that we used to calculate EpsII
NewTau = ComputeCreepLaw_TauII(EpsII,Phase.CreepLaws[1],c)
#@show NewTau
@test NewTau ≈ TauII.val


# Given stress
#@test ComputeCreepLaw_EpsII(1e6Pa, x1, CreepLawParams())==5e-13/s                # dimensional input       

# Given strainrate 
#@test ComputeCreepLaw_TauII(1e-13/s, x1, CreepLawParams())==1e18*2*1e-13Pa       # dimensional input       
# -------------------------------------------------------------------


end
