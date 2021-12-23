using Test
using GeoParams
@testset "EnergyParameters.jl" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=10km);
           
# Heat capacity ---------

# Constant heat capacity
cp1      =   ConstantHeatCapacity()
@test cp1.cp.val == 1050.0

Nondimensionalize!(cp1,CharUnits_GEO)
@test cp1.cp.val ≈ 1.3368075000000002e22

@test ComputeHeatCapacity(100.0, cp1) ≈ 1.3368075000000002e22  # compute


# Temperature-dependent heat capacity
# dimensional
T        =   250:100:1250;
cp2      =   T_HeatCapacity_Whittacker()
Cp       =   ComputeHeatCapacity(T,cp2)
@test sum(Cp) ≈ 11667.035717418683

# nondimensional
cp2_nd   =   T_HeatCapacity_Whittacker()
Nondimensionalize!(cp2_nd,CharUnits_GEO)
T_nd     =   Float64.(T*K/CharUnits_GEO.Temperature)
Cp_nd    =   ComputeHeatCapacity(T_nd,cp2_nd)
@test sum(Cp_nd) ≈ 1.4853886523631602e23

# Dimensionalize again and double-check the results
@test sum(abs.(ustrip.(Cp_nd*CharUnits_GEO.heatcapacity) - Cp)) < 1e-11

# Test with arrays
T_array     =  ustrip.(T)*ones(10)'
Cp_array    =  zeros(size(T_array))
P_array     =  zeros(size(T_array))
ComputeHeatCapacity!(Cp_array, P_array, T_array,cp1)
@test Cp_array[1] ≈ 1.3368075000000002e22

T_array     =  ustrip.(T)*ones(10)'
Cp_array    =  zeros(size(T_array))
P_array     =  zeros(size(T_array))
ComputeHeatCapacity!(Cp_array, P_array, T_array,cp2)
@test sum(Cp_array[:,1]) ≈ 11667.035717418683

T_array     =  ustrip.(T)*ones(10)'
Cp_array    =  zeros(size(T_array))
P_array     =  zeros(size(T_array))
ComputeHeatCapacity!(Cp_array, T_array,cp2)
@test sum(Cp_array[:,1]) ≈ 11667.035717418683


# Check that it works if we give a phase array
MatParam    =   Array{MaterialParams, 1}(undef, 2);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                    HeatCapacity  = ConstantHeatCapacity());

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                    HeatCapacity  = T_HeatCapacity_Whittacker());

# test computing material properties
n = 100;
Phases              = ones(Int64,n,n,n);
Phases[:,:,20:end] .= 2

Cp = zeros(size(Phases))
T  =  ones(size(Phases))*1500
P  =  zeros(size(Phases))

ComputeHeatCapacity!(Cp, Phases, P, T, MatParam)       # computation routine w/out P (not used in most heat capacity formulations)
@test sum(Cp[1,1,:]) ≈ 121399.0486067196

ComputeHeatCapacity!(Cp, Phases, T, MatParam)       # computation routine w/out P (not used in most heat capacity formulations)
@test sum(Cp[1,1,:]) ≈ 121399.0486067196

# test if we provide phase ratios
PhaseRatio  = zeros(n,n,n,3);
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz)
    PhaseRatio[I] = 1.0  
end
ComputeHeatCapacity!(Cp, PhaseRatio, P, T, MatParam)
@test sum(Cp[1,1,:]) ≈ 121399.0486067196


# -----------------------

# Conductivity ----------

# Constant 

# Constant conductivity
cond      =   ConstantConductivity()
@test Value(cond.k) == 3.0
@test cond.k.unit==u"W"/K/m

Nondimensionalize!(cond,CharUnits_GEO)
@test  Value(cond.k) ≈ 3.8194500000000007

@test ComputeConductivity(100.0, cond) ≈ 3.8194500000000007 # compute

# Temperature-dependent conductivity
# dimensional
T        =   (250:100:1250);
cond2    =   T_Conductivity_Whittacker()
k        =   ComputeConductivity(T,cond2)
@test sum(k) ≈ 27.503366436682285

# nondimensional
cond2_nd =   T_Conductivity_Whittacker()
Nondimensionalize!(cond2_nd,CharUnits_GEO)
T_nd     =   Float64.(ustrip(T/CharUnits_GEO.Temperature))
k_nd     =   ComputeConductivity(T_nd,cond2_nd)
@test sum(k_nd) ≈ 35.01591097886205

# Dimensionalize again and double-check the results
@test sum(abs.(ustrip.(k_nd*CharUnits_GEO.conductivity) - k)) < 1e-11

# Check if we use arrays
T_array     =  ustrip.(T)*ones(100)'
k_array     =  copy(T_array)
P_array     =  copy(T_array)

ComputeConductivity!(k_array,P_array,T_array, cond)
@test k_array[1] ≈ 3.8194500000000007

ComputeConductivity!(k_array,P_array,T_array, cond2)
@test sum(k_array) ≈ 2750.3366436682285

k_TP    =   Set_TP_Conductivity["LowerCrust"]
ComputeConductivity!(k_array, P_array, T_array, k_TP)
@test sum(k_array) ≈ 2055.7129327367625

# Check that it works if we give a phase array
MatParam    =   Array{MaterialParams, 1}(undef, 3);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                    Conductivity  = ConstantConductivity());

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                    Conductivity  = T_Conductivity_Whittacker());

MatParam[3] =   SetMaterialParams(Name="MantleLithosphere", Phase=3,
                    Conductivity  = Set_TP_Conductivity["Mantle"]);

# test computing material properties
n = 100;
Phases              = ones(Int64,n,n,n);
Phases[:,:,20:end] .= 2
Phases[:,:,60:end] .= 3

PhaseRatio  = zeros(n,n,n,3);
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz)
    PhaseRatio[I] = 1.0  
end

k   =  zeros(size(Phases))
T   =  ones(size(Phases))*1500
P   =  zeros(size(Phases))

ComputeConductivity!(k, Phases, P, T, MatParam) 
@test sum(k) ≈ 1.9216938849389635e6

ComputeConductivity!(k, PhaseRatio, P, T, MatParam) 
@test sum(k) ≈ 1.9216938849389635e6

#


######

# TP-dependent conductivity for different predefines cases
T       =   (250:100:1250);
P       = 1e6*ones(size(T))/ustrip(uconvert(Pa,1MPa))  # must be in MPa!
List    = ["LowerCrust"   "Mantle"        "OceanicCrust"  "UpperCrust"]
Sol_kT  = [20.55712932736763 28.700405819019323 20.55712932736763 19.940302462417037]
for i=1:length(List)
    k_TP    =   Set_TP_Conductivity[List[i]]
    k       =   ComputeConductivity(P,T,k_TP)           # note that P must be in MPa
    @test sum(k) ≈ Sol_kT[i]

    k_TP_nd  =  deepcopy(k_TP)
    Nondimensionalize!(k_TP_nd,CharUnits_GEO)
    T_nd     =   Float64.(ustrip.(T/CharUnits_GEO.Temperature))
    P_nd     =   Float64.(ustrip(P/CharUnits_GEO.stress))
    k_nd     =   ComputeConductivity(P_nd,T_nd,k_TP_nd)

    @test ustrip(sum(abs.(ustrip.(k_nd*CharUnits_GEO.conductivity) - k))) < 1e-11

end


T = [200 300; 400 500]
k1        =   ComputeConductivity(T,cond2)


# -----------------------


# Latent heat -----------
a = ConstantLatentHeat()
Q_L = ComputeLatentHeat(a)
@test Q_L == 400

Nondimensionalize!(a,CharUnits_GEO)
Q_L = ComputeLatentHeat(a)
@test Q_L ≈ 4e21
# -----------------------

# Radioactive heat ------
a = ConstantRadioactiveHeat()
H_r = ComputeRadioactiveHeat(a)
@test H_r ≈ 1.0e-6

Nondimensionalize!(a,CharUnits_GEO)
H_r = ComputeRadioactiveHeat(a)
@test H_r == 0.1
# -----------------------


# Shear heating -------
Χ       = ConstantShearheating(1.0)
 
# Define parameters as vectors
τ       = [1 2 3 4]*1e6    
ε       = [1 0.1 0.1 1]   
ε_el    = [0.01 0.01 0.01 0.01]

τ_2D       = [1 2; 3 4]*1e6     
ε_2D       = [1 0.1; 0.1 1]   
ε_el_2D    = [0.01 0.01; 0.01 0.01]   

# With elasticity
H_s1 = ComputeShearheating(τ,   ε,    ε_el,     Χ)
H_s2 = ComputeShearheating(τ_2D,ε_2D, ε_el_2D,  Χ)
@test H_s1 ≈ 5.4e6
@test H_s2 ≈ 5.4e6

# No elasticity
H_s3 = ComputeShearheating(τ,   ε,     Χ)
H_s4 = ComputeShearheating(τ_2D,ε_2D,  Χ)
@test H_s3 ≈ 5.5e6
@test H_s4 ≈ 5.5e6

# Now in non-dimensional units
τ       = [1 2 3 4]     
ε       = [1 0.1 0.1 1]  
ε_el    = [0.01 0.01 0.01 0.01]   

τ_2D       = [1 2; 3 4]     
ε_2D       = [1 0.1; 0.1 1]   
ε_el_2D    = [0.01 0.01; 0.01 0.01]  
Nondimensionalize!(Χ,CharUnits_GEO)

H_s1 = ComputeShearheating(τ,   ε,    ε_el,     Χ)
H_s2 = ComputeShearheating(τ_2D,ε_2D, ε_el_2D,  Χ)
H_s3 = ComputeShearheating(τ,   ε,     Χ)
H_s4 = ComputeShearheating(τ_2D,ε_2D,  Χ)
@test H_s1 ≈ 5.4
@test H_s2 ≈ 5.4
@test H_s3 ≈ 5.5
@test H_s4 ≈ 5.5
# -----------------------

end

