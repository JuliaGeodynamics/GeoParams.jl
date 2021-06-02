using Test
using GeoParams


# This tests the MaterialParameters structure

# Define a linear viscous creep law
x1  =   GeoParams.MaterialParameters.LinearViscous()

x2  =   GeoParams.MaterialParameters.LinearViscous(eta=1e3)

# Define a struct for a first phase
Phase = MaterialParams(Name="test1",
                CreepLaws=(GeoParams.MaterialParameters.LinearViscous, GeoParams.MaterialParameters.LinearViscous(eta=1e20)),
                Density = 1000kg/m^3);

# Evaluate viscous strainrate, given a stress

