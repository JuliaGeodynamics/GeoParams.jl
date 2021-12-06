using Test
using GeoParams
@testset "Density.jl" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
                
# Define a linear viscous creep law
x1      =   ConstantDensity(ρ=2900kg/m^3)
@test x1.ρ.val == 2900kg/m^3

Nondimensionalize!(x1,CharUnits_GEO)
@test x1.ρ.val ≈ 2.9e-16


x2      =   PT_Density()
@test x2.α.val==3e-5/K
@test x2.ρ0.val==2900kg/m^3

Nondimensionalize!(x2,CharUnits_GEO)
@test x2.T0.val≈0.21454659702313156


# Compute with density
@test ComputeDensity(1.0,1.0, x2) ≈ 2.8419999999999996e-16
@test ComputeDensity(1.0,1.0, x1) ≈ 2.9e-16

# Read Phase diagram interpolation object
fname   =   "./test_data/Peridotite.in"
PD_data =   Read_LaMEM_Perple_X_Diagram(fname);
@test PD_data.meltFrac(1500,1e7) ≈ 0.24368492372485706
@test PD_data.Rho(1500,1e7) ≈ 3042.836820256982
@test PD_data.meltRho(1500,1e7) ≈ 2662.227167592414
@test PD_data.rockRho(1500,1e7) ≈ 3165.467673917775

end