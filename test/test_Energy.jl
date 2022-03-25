using Test
using GeoParams
@testset "EnergyParameters.jl" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=10km);
           
# Heat capacity ---------

# Constant heat capacity
cp1      =  ConstantHeatCapacity()
info     =  param_info(cp1)
@test isbits(cp1)
@test cp1.cp.val == 1050.0

cp1     = nondimensionalize(cp1,CharUnits_GEO)
@test cp1.cp.val ≈ 1.3368075000000002e22
@test compute_heatcapacity(cp1, 10.0, 100.0) ≈ 1.3368075000000002e22  # compute

# Temperature-dependent heat capacity
# dimensional
T        =   250.0:100:1250;
cp2      =   T_HeatCapacity_Whittington()
Cp       =   similar(T)
@test isbits(cp2)
compute_heatcapacity!(Cp, cp2, T)
@test sum(Cp) ≈ 11667.035717418683

# nondimensional
cp2_nd   =   T_HeatCapacity_Whittington()
cp2_nd   =   nondimensionalize(cp2_nd,CharUnits_GEO)
T_nd     =   Float64.(T*K/CharUnits_GEO.Temperature)
Cp_nd    =   similar(T)
compute_heatcapacity!(Cp_nd, cp2_nd, T_nd)
@test sum(Cp_nd) ≈ 1.4853886523631602e23

# Dimensionalize again and double-check the results
@test sum(abs.(ustrip.(Cp_nd*CharUnits_GEO.heatcapacity) - Cp)) < 1e-11

# Test with arrays
T_array     =  T*ones(10)'
Cp_array    =  similar(T_array)
P_array     =  similar(T_array)
compute_heatcapacity!(Cp_array, cp1, P_array, T_array)
@test Cp_array[1] ≈ 1.3368075000000002e22

Cp_array    =  similar(T_array)
P_array     =  similar(T_array)
compute_heatcapacity!(Cp_array, cp2, P_array, T_array)
@test sum(Cp_array[:,1]) ≈ 11667.035717418683

T_array     =  T*ones(10)'
Cp_array    =  zeros(size(T_array))
P_array     =  zeros(size(T_array))
compute_heatcapacity!(Cp_array, cp2, T_array)
@test sum(Cp_array[:,1]) ≈ 11667.035717418683


# Check that it works if we give a phase array
MatParam    =   Array{MaterialParams, 1}(undef, 2);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                    HeatCapacity  = ConstantHeatCapacity());

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                    HeatCapacity  = T_HeatCapacity_Whittington());

Mat_tup = Tuple(MatParam)

Mat_tup = Tuple(MatParam)

# test computing material properties
n = 100;
Phases              = ones(Int64,n,n,n);
Phases[:,:,20:end] .= 2

Cp = zeros(size(Phases))
T  =  ones(size(Phases))*1500
P  =  zeros(size(Phases))

compute_heatcapacity!(Cp, Mat_tup, Phases, P, T)    # computation routine w/out P (not used in most heat capacity formulations)     
num_alloc = @allocated compute_heatcapacity!(Cp, Mat_tup, Phases, P, T)
@test sum(Cp[1,1,:]) ≈ 121399.0486067196
@test num_alloc == 0

compute_heatcapacity!(Cp, Mat_tup, Phases, T)       # computation routine w/out P (not used in most heat capacity formulations)
@test sum(Cp[1,1,:]) ≈ 121399.0486067196

# test if we provide phase ratios
PhaseRatio  = zeros(n,n,n,3);
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz)
    PhaseRatio[I] = 1.0  
end
compute_heatcapacity!(Cp, Mat_tup, PhaseRatio, P, T)
num_alloc = @allocated compute_heatcapacity!(Cp, Mat_tup, PhaseRatio, P, T)
@test sum(Cp[1,1,:]) ≈ 121399.0486067196
@test num_alloc == 0


# -----------------------

# Conductivity ----------

# Constant 

# Constant conductivity
cond      =   ConstantConductivity()
@test isbits(cond)
@test NumValue(cond.k) == 3.0
@test cond.k.unit==u"W"/K/m

cond = nondimensionalize(cond,CharUnits_GEO)
@test  NumValue(cond.k) ≈ 3.8194500000000007

@test compute_conductivity(cond, 100.0) ≈ 3.8194500000000007 # compute

# Temperature-dependent conductivity
# dimensional
T        =   Vector{Float64}(250:100:1250);
cond2    =   T_Conductivity_Whittington()
k        =   compute_conductivity(cond2, T)
@test isbits(cond2)
@test sum(k) ≈ 27.503366436682285

# nondimensional
cond2_nd =   T_Conductivity_Whittington()
cond2_nd =   nondimensionalize(cond2_nd,CharUnits_GEO)
T_nd     =   Float64.(ustrip.(T/CharUnits_GEO.Temperature))
k_nd     =   compute_conductivity(cond2_nd, T_nd)
@test sum(k_nd) ≈ 35.01591097886205

# Dimensionalize again and double-check the results
@test sum(abs.(ustrip.(k_nd*CharUnits_GEO.conductivity) - k)) < 1e-11



# Temperature-dependent parameterised conductivity
# dimensional
T        =   Vector{Float64}(250:100:1250);
cond2    =   T_Conductivity_Whittington_parameterised()
k        =   compute_conductivity(cond2, T)
@test isbits(cond2)
@test sum(k) ≈ 27.553653387829254

# nondimensional
cond2_nd =   T_Conductivity_Whittington_parameterised()
cond2_nd =   nondimensionalize(cond2_nd,CharUnits_GEO)
T_nd     =   Float64.(ustrip.(T/CharUnits_GEO.Temperature))
k_nd     =   compute_conductivity(cond2_nd, T_nd)
@test sum(k_nd) ≈ 35.079933810714806

# Dimensionalize again and double-check the results
@test sum(abs.(ustrip.(k_nd*CharUnits_GEO.conductivity) - k)) < 1e-11


# Check if we use arrays
T_array     =  ustrip.(T)*ones(100)'
k_array     =  copy(T_array)
P_array     =  copy(T_array)

compute_conductivity!(k_array,cond,P_array,T_array)
@test k_array[1] ≈ 3.8194500000000007

compute_conductivity!(k_array,cond2,P_array,T_array) 
@test sum(k_array) ≈ 2755.3653387829254

k_TP    =   Set_TP_Conductivity("LowerCrust")
compute_conductivity!(k_array, k_TP, P_array, T_array) 
@test sum(k_array) ≈ 2055.7129327367625

# Check that it works if we give a phase array
MatParam    =   Array{MaterialParams, 1}(undef, 3);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                    Conductivity  = ConstantConductivity());

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                    Conductivity  = T_Conductivity_Whittington());

MatParam[3] =   SetMaterialParams(Name="MantleLithosphere", Phase=3,
                    Conductivity = Set_TP_Conductivity("Mantle"));

Mat_tup = Tuple(MatParam)
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

compute_conductivity!(k, Mat_tup, Phases, P, T) 
num_alloc = @allocated compute_conductivity!(k, Mat_tup, Phases, P, T) 
@test sum(k) ≈ 1.9216938849389635e6
@test num_alloc == 0

compute_conductivity!(k, Mat_tup, PhaseRatio, P, T) 
num_alloc = @allocated compute_conductivity!(k, Mat_tup, PhaseRatio, P, T) 
@test num_alloc == 0
@test sum(k) ≈ 1.9216938849389635e6


#


######

# TP-dependent conductivity for different predefines cases
T       =   Vector{Float64}(250:100:1250);
P       = 1e6*ones(size(T))/ustrip(uconvert(Pa,1MPa))  # must be in MPa!
List    = ["LowerCrust"   "Mantle"        "OceanicCrust"  "UpperCrust"]
Sol_kT  = [20.55712932736763 28.700405819019323 20.55712932736763 19.940302462417037]
for i=1:length(List)
    k_TP    =   Set_TP_Conductivity(List[i])
    k       =   compute_conductivity(k_TP,P,T)           # note that P must be in MPa
    @test sum(k) ≈ Sol_kT[i]

    k_TP_nd  =   deepcopy(k_TP)
    k_TP_nd  =   nondimensionalize(k_TP_nd,CharUnits_GEO)
    T_nd     =   Float64.(ustrip.(T/CharUnits_GEO.Temperature))
    P_nd     =   Float64.(ustrip(P/CharUnits_GEO.stress))
    k_nd     =   compute_conductivity(k_TP_nd,P_nd,T_nd)

    @test ustrip(sum(abs.(ustrip.(k_nd*CharUnits_GEO.conductivity) - k))) < 1e-11

end


T = [200. 300.; 400. 500.]
k1        =   compute_conductivity(cond2,T)


# -----------------------


# Latent heat -----------
a = ConstantLatentHeat()
Q_L = compute_latent_heat(a)
@test isbits(a)
@test Q_L == 400

a   = nondimensionalize(a,CharUnits_GEO)
Q_L = compute_latent_heat(a)
@test Q_L ≈ 4e21
# -----------------------

# Radioactive heat ------
a = ConstantRadioactiveHeat()
H_r = compute_radioactive_heat(a)
@test isbits(a)
@test H_r ≈ 1.0e-6

a = nondimensionalize(a,CharUnits_GEO)
H_r = compute_radioactive_heat(a)
@test H_r == 0.1
# -----------------------


# Shear heating -------
Χ       = ConstantShearheating(1.0)
@test isbits(Χ)
 
# Define parameters as vectors
τ       = [1 2 3 4]*1e6    
ε       = [1 0.1 0.1 1]   
ε_el    = [0.01 0.01 0.01 0.01]

τ_2D       = [1 2; 3 4]*1e6     
ε_2D       = [1 0.1; 0.1 1]   
ε_el_2D    = [0.01 0.01; 0.01 0.01]   

# With elasticity
H_s1 = compute_shearheating(Χ, τ,   ε,    ε_el   )
H_s2 = compute_shearheating(Χ, τ_2D,ε_2D, ε_el_2D)
@test H_s1 ≈ 5.4e6
@test H_s2 ≈ 5.4e6

# No elasticity
H_s3 = compute_shearheating(Χ, τ,   ε   )
H_s4 = compute_shearheating(Χ, τ_2D,ε_2D)
@test H_s3 ≈ 5.5e6
@test H_s4 ≈ 5.5e6

# Now in non-dimensional units
τ       = [1 2 3 4]     
ε       = [1 0.1 0.1 1]  
ε_el    = [0.01 0.01 0.01 0.01]   

τ_2D       = [1 2; 3 4]     
ε_2D       = [1 0.1; 0.1 1]   
ε_el_2D    = [0.01 0.01; 0.01 0.01]  
Χ       = nondimensionalize(Χ,CharUnits_GEO)
info    = param_info(Χ)

H_s1 = compute_shearheating(Χ, τ,   ε,    ε_el,  )
H_s2 = compute_shearheating(Χ, τ_2D,ε_2D, ε_el_2D)
H_s3 = compute_shearheating(Χ, τ,   ε)
H_s4 = compute_shearheating(Χ, τ_2D,ε_2D)
@test H_s1 ≈ 5.4
@test H_s2 ≈ 5.4
@test H_s3 ≈ 5.5
@test H_s4 ≈ 5.5
# -----------------------

end

