using Test
using GeoParams


# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);

# Define a struct for a first phase
Phase = SetMaterialParams(Name="test1", Phase=22,
                CreepLaws= (PowerlawViscous(), LinearViscous(η=1e21Pa*s)),
                Gravity  = ConstantGravity(g=11m/s^2),
                Density  = ConstantDensity(),
                CharDim  = CharUnits_GEO);

@test   Phase.Density[1].ρ.val==2.9e-16
@test   Phase.Gravity[1].g*1.0==1.1000000000000002e19
@test   Phase.CreepLaws[1].η0*1.0==0.1
@test   Phase.CreepLaws[2].η*1.0==100.0

Phase1 = SetMaterialParams(Name="test1",
                CreepLaws= (LinearViscous(), LinearViscous(η=1e21Pa*s)),
                Density  = ConstantDensity());

#GeoParams.MaterialParameters.Dimensionalize!(Phase, CharUnits_GEO)
