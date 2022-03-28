using Test
using GeoParams, LinearAlgebra
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

#------------------------------
# 5th order polynomial
p        =  MeltingParam_5thOrder();
compute_meltfraction!(phi_dim, p, zeros(size(T)), ustrip.(T))
@test sum(phi_dim) ≈ 4.708427909521561

# experimental data to create the fit
data = [
1115 100;
1050 90;
969 85;
932 80;
907 70;
880 60;
850 54;
825 51.7;
800 52.9;
775 46.3;
750 44.9;
725 29.9;
700 14.9;
690 0
] ;

data[:,2] = data[:,2]/100;
Tdata = data[:,1] .+ 273.15;
phi = zeros(size(Tdata))
compute_meltfraction!(phi, p, zeros(size(Tdata)), ustrip.(Tdata))

@test norm(data[:,2]-phi) ≈ 0.07151515017819135
#------------------------------



#------------------------------
# 4th order polynomial
p        =  MeltingParam_4thOrder();
compute_meltfraction!(phi_dim, p, zeros(size(T)), ustrip.(T))
@test sum(phi_dim) ≈ 4.853749635538406

# experimental data to create the fit
data = [
1000 100;
990 100;
975 93;
950 89.2;
925 76.3;
900 69.6;
875 59;
850 54;
825 51.7;
800 52.9;
775 46.3;
750 44.9;
725 29.9;
700 14.9;
690 0] ;

data[:,2] = data[:,2]/100;
Tdata = data[:,1] .+ 273.15;
phi = zeros(size(Tdata))
compute_meltfraction!(phi, p, zeros(size(Tdata)), ustrip.(Tdata))

@test norm(data[:,2]-phi) ≈ 0.0678052542705406
#------------------------------


#------------------------------
# Quadratic parameterisation
p        =  MeltingParam_Quadratic();
compute_meltfraction!(phi_dim, p, zeros(size(T)), ustrip.(T))
@test sum(phi_dim) ≈ 5.0894901144641
#------------------------------


# Test computation of melt parameterization for the whole computational domain, using arrays 
MatParam    =   Array{MaterialParams, 1}(undef, 4);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                        Melting  = PerpleX_LaMEM_Diagram("test_data/Peridotite.in"));

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                        Melting   = MeltingParam_Caricchi());

# No melting parameterization for this phase
MatParam[3] =   SetMaterialParams(Name="UpperCrust", Phase=3,
                        Melting   = MeltingParam_5thOrder(),
                        Density   = PT_Density());

# No melting parameterization for this phase
MatParam[4] =   SetMaterialParams(Name="LowerCrust", Phase=4,
                        Density   = PT_Density());

Mat_tup = Tuple(MatParam)

# test computing material properties
n = 100;
Phases              = ones(Int64,n,n,n);
Phases[:,:,20:end] .= 2
Phases[:,:,80:end] .= 3
Phases[:,:,90:end] .= 4

ϕ = zeros(size(Phases))
T =  ones(size(Phases))*1500
P =  ones(size(Phases))*10

compute_meltfraction!(ϕ, Mat_tup, Phases, P,T) #allocations coming from computing meltfraction using PhaseDiagram_LookupTable

@test sum(ϕ)/n^3 ≈ 0.7463001302812086


# test computing material properties when we have PhaseRatios, instead of Phase numbers
PhaseRatio  = zeros(n,n,n,4);
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz)
    PhaseRatio[I] = 1.0  
end

compute_meltfraction!(ϕ, Mat_tup, PhaseRatio, P,T)
@test sum(ϕ)/n^3 ≈ 0.7463001302812086


end

