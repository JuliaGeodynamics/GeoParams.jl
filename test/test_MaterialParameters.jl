using Test
using GeoParams


# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);

# Define a struct for a first phase
Phase = SetMaterialParams(Name="test1",
                CreepLaws= (LinearViscous(), LinearViscous(eta=1e21Pa*s)),
                Density  =  1000kg/m^3,
                CharDim=CharUnits_GEO);

Phase1 = SetMaterialParams(Name="test1",
                CreepLaws= (LinearViscous(), LinearViscous(eta=1e21Pa*s)),
                Density  =  1000kg/m^3);

#GeoParams.MaterialParameters.Nondimensionalize!(Phase, CharUnits_GEO)
