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
# Test data and example from Hirth & Kohlstedt, (2000), table 1, footnote e 
CharDim = GEO_units(;length=1000km, stress=100MPa, viscosity=1Pas, temperature=1000C)
EpsII   = GeoUnit(0.0s^-1.0)
TauII   = GeoUnit(0.3MPa)
P       = GeoUnit(1000MPa)
T       = GeoUnit(1673K)
f       = GeoUnit(1000NoUnits)
z       = CreepLawVariables(P=nondimensionalize(P,CharDim),T=nondimensionalize(T,CharDim),f=nondimensionalize(f,CharDim))
Phase   = SetMaterialParams(Name="Viscous Matrix", Phase=2,
                                     Density   = ConstantDensity(),
                                     CreepLaws = DislocationCreep(n=3.5NoUnits, r=1.2NoUnits, A=90MPa^(-3.5)/s,E=480000.0J/mol,V=11e-6m^(3.0)/mol,Apparatus=3), CharDim = CharDim)

# calculate strainrate with non-dimensional values, create GeoUnit() with the calculated non-dim. value and also dimensionalize it back
EpsII = computeCreepLaw_EpsII(nondimensionalize(TauII,CharDim),Phase.CreepLaws[1],z)
EpsII_geo = GeoUnit(EpsII, s^(-1.0), false)
EpsII_dim = dimensionalize(EpsII_geo, CharDim)

# check whether those calculated non-dim. and dim. values agree with example calculation from Hirth&Kohlstedt, 2000 for disl. creep
@test EpsII  ≈ nondimensionalize(2.4748138964055293e-12s^(-1.0), CharDim) rtol = 1e-12
@test EpsII_dim.val ≈ 2.4748138964055293e-12 rtol = 1e-12

# Check that once inverted, we get back the TauII that we used to calculate EpsII, for both non-dim. and dim. NewTau
NewTau = computeCreepLaw_TauII(EpsII_geo,Phase.CreepLaws[1],z)
NewTau_geo = GeoUnit(NewTau, MPa, false)
NewTau_dim = dimensionalize(NewTau_geo, CharDim)
@test NewTau_geo ≈ nondimensionalize(0.3MPa,CharDim) rtol = 1e-12
@test NewTau_dim ≈ TauII.val


# Given strainrate 
#@test computeCreepLaw_EpsII(1e-13/s, x1, CreepLawParams())==1e18*2*1e-13Pa       # dimensional input       
# -------------------------------------------------------------------

end
