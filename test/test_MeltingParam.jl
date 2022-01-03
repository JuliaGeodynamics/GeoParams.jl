using Test
using GeoParams
@testset "MeltingParam.jl" begin

# This tests the various melting parameterizations
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=10km);
       

T   =   (250:100:1250)*K .+ 273.15K
T_nd= Float64.(T/CharUnits_GEO.Temperature)

# Caricchi parameterization [in ND numbers, which is anyways the typical use case]
p        =  MeltingParam_Caricchi()
phi_dim  =  ComputeMeltingParam(0,T, p)
phi_dim1  = zeros(size(phi_dim))

ComputeMeltingParam!(phi_dim1, 0,T, p) # in-place routine


p_nd     =  p
p_nd     =  Nondimensionalize(p_nd, CharUnits_GEO)
phi_nd   =  ComputeMeltingParam(0,T_nd, p_nd)

# Do this computation manually, using the actual expression of Caricchi
T_C         =   250:100:1250    # in celcius
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

ComputeMeltingParam!(ϕ, Phases, P,T, MatParam)
@test sum(ϕ)/n^3 ≈ 0.6463001302812086


# test computing material properties when we have PhaseRatios, instead of Phase numbers
PhaseRatio  = zeros(n,n,n,3);
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz)
    PhaseRatio[I] = 1.0  
end

ComputeMeltingParam!(ϕ, PhaseRatio, P,T, MatParam)
@test sum(ϕ)/n^3 ≈ 0.6463001302812086


end

