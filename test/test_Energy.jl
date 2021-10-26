using Test
using GeoParams
@testset "EnergyParameters.jl" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=10km);
           
# Heat capacity ---------

# Constant heat capacity
cp1      =   ConstantHeatCapacity()
@test cp1.cp.val == 1050J/kg/K

Nondimensionalize!(cp1,CharUnits_GEO)
@test cp1.cp.val ≈ 1.3368075000000002e22

@test ComputeHeatCapacity(100.0, cp1) ≈ 1.3368075000000002e22  # compute


# Temperature-dependent heat capacity
# dimensional
T        =   (250:100:1250)*K;
cp2      =   T_HeatCapacity_Whittacker()
Cp       =   ComputeHeatCapacity(T,cp2)
@test sum(Cp) ≈ 11667.035717418683J/kg/K

# nondimensional
cp2_nd   =   T_HeatCapacity_Whittacker()
Nondimensionalize!(cp2_nd,CharUnits_GEO)
T_nd     =   Float64.(T/CharUnits_GEO.Temperature)
Cp_nd    =   ComputeHeatCapacity(T_nd,cp2_nd)
@test sum(Cp_nd) ≈ 1.4853886523631602e23

# Dimensionalize again and double-check the results
@test ustrip(sum(abs.(Cp_nd*CharUnits_GEO.heatcapacity - Cp))) < 1e-11

# -----------------------

# Conductivity ----------

# Constant 

# Constant conductivity
cond      =   ConstantConductivity()
@test Value(cond.k) == 3.0Watt/K/m

Nondimensionalize!(cond,CharUnits_GEO)
@test  Value(cond.k) ≈ 3.8194500000000007

@test ComputeConductivity(100.0, cond) ≈ 3.8194500000000007 # compute

# Temperature-dependent conductivity
# dimensional
T        =   (250:100:1250)*K;
cond2    =   T_Conductivity_Whittacker()
k        =   ComputeConductivity(T,cond2)
@test sum(k) ≈ 27.503366436682285Watt/K/m

# nondimensional
cond2_nd =   T_Conductivity_Whittacker()
Nondimensionalize!(cond2_nd,CharUnits_GEO)
T_nd     =   Float64.(T/CharUnits_GEO.Temperature)
k_nd     =   ComputeConductivity(T_nd,cond2_nd)
@test sum(k_nd) ≈ 35.01591097886205

# Dimensionalize again and double-check the results
@test ustrip(sum(abs.(k_nd*CharUnits_GEO.conductivity - k))) < 1e-11


# -----------------------

end