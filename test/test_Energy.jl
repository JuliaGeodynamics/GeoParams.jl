using Test
using GeoParams
@testset "EnergyParameters.jl" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1km);
           
# Heat capacity ---------

# Constant heat capacity
cp1      =   ConstantHeatCapacity()
@test cp1.cp.val == 1050J/mol/K

Nondimensionalize!(cp1,CharUnits_GEO)
@test cp1.cp.val ≈ 1.3368075000000003e-10

@test ComputeHeatCapacity(100.0, cp1) ≈ 1.3368075000000003e-10  # compute


# Temperature-dependent heat capacity
# dimensional
T        =   (250:100:1250)*K;
cp2      =   T_HeatCapacity_Whittacker()
Cp       =   ComputeHeatCapacity(T,cp2)
@test sum(Cp) ≈ 2587.515181409116J/K/mol

# nondimensional
cp2_nd   =   T_HeatCapacity_Whittacker()
Nondimensionalize!(cp2_nd,CharUnits_GEO)
T_nd     =   Float64.(T/CharUnits_GEO.Temperature)
Cp_nd    =   ComputeHeatCapacity(T_nd,cp2_nd)
@test sum(Cp_nd) ≈ 3.294294953211017e-10

# Dimensionalize again and double-check the results
@test ustrip(sum(abs.(Cp_nd*CharUnits_GEO.heatcapacity - Cp))) < 1e-12


end