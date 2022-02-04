using Test
using GeoParams
@testset "MeltingParam.jl" begin

#Make sure structure is isbits
x = MeltingParam_Caricchi()
@test isbits(x)

# This tests the various melting parameterizations
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=10km);
       

T       =   Vector(250:100:1250)*K .+ 273.15K
T_nd    =   Float64.(T/CharUnits_GEO.Temperature)

# Caricchi parameterization [in ND numbers, which is anyways the typical use case]
p        =  MeltingParam_Caricchi()
phi_dim  =  zeros(size(T))
compute_meltfraction!(phi_dim, p, zeros(size(T)), ustrip.(T))

phi_dim1  = zeros(size(phi_dim))
compute_meltfraction!(phi_dim1, p, zeros(size(T)), ustrip.(T)) # in-place routine

p_nd     =  p
p_nd     =  nondimensionalize(p_nd, CharUnits_GEO)
phi_nd   = zeros(size(T))
compute_meltfraction!(phi_nd, p_nd, zeros(size(T_nd)),T_nd)

# Do this computation manually, using the actual expression of Caricchi
T_C         =   Vector(250:100:1250)    # in celcius
Phi_solid   =   1.0 .- 1.0./(1.0 .+ exp.((800.0 .- T_C)./23.0)); 
Phi_anal    =   1.0 .- Phi_solid

@test sum(phi_dim  - Phi_anal)   < 1e-12
@test sum(phi_dim1 - Phi_anal)   < 1e-12
@test sum(phi_nd - Phi_anal)    < 1e-12


# Test computation of density for the whole computational domain, using arrays 
MatParam    =   Array{MaterialParams, 1}(undef, 3);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                        Melting  = PerpleX_LaMEM_Diagram("test_data/Peridotite.in"));

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                        Melting   = MeltingParam_Caricchi());

# No melting parameterization for this phase
MatParam[3] =   SetMaterialParams(Name="UpperCrust", Phase=3,
                        Density   = PT_Density());

# test computing material properties
n = 100;
Phases              = ones(Int64,n,n,n);
Phases[:,:,20:end] .= 2
Phases[:,:,80:end] .= 3

ϕ = zeros(size(Phases))
T =  ones(size(Phases))*1500
P =  ones(size(Phases))*10

compute_meltfraction!(ϕ, MatParam, Phases, P,T)
@test sum(ϕ)/n^3 ≈ 0.6463001302812086


# test computing material properties when we have PhaseRatios, instead of Phase numbers
PhaseRatio  = zeros(n,n,n,3);
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz)
    PhaseRatio[I] = 1.0  
end

compute_meltfraction!(ϕ, MatParam, PhaseRatio, P,T)
@test sum(ϕ)/n^3 ≈ 0.6463001302812086


end

