using Test
using GeoParams
@testset "SeismicVelocity.jl" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
                
# Define constant velocities 
x1      =   ConstantSeismicVelocity(Vp=8.05km/s, Vs=3.5km/s)
@test  Value(x1.Vp)==8.05km/s
@test  Value(x1.Vs)==3.5km/s

x1 = Nondimensionalize(x1,CharUnits_GEO)
@test NumValue(x1.Vp) ≈ 8.050000000000001e9
@test NumValue(x1.Vs) ≈ 3.5e9

# Compute
@test ComputePwaveVelocity(1.0,1.0, x1) ≈ 8.050000000000001e9
@test ComputeSwaveVelocity(1.0,1.0, x1) ≈ 3.5e9

# Read Phase diagram interpolation object
fname   =   "./test_data/Peridotite.in"
PD_data =   PerpleX_LaMEM_Diagram(fname);
@test PD_data.Vp(1500,1e7) ≈ 6.5290725233303935
@test PD_data.Vs(1500,1e7) ≈ 2.4874400647487658

@test ComputePwaveVelocity(1e7, 1500, PD_data) ≈  6.5290725233303935
@test ComputeSwaveVelocity(1e7, 1500, PD_data) ≈  2.4874400647487658

# Do the same but non-dimensionalize the result
CharDim  =  GEO_units();
PD_data1 =  PerpleX_LaMEM_Diagram(fname, CharDim=CharDim);

rho_ND   =  PD_data1.Rho(Nondimensionalize(1500K,CharDim),Nondimensionalize(1e8*Pa,CharDim)) 
Vp_ND    =  PD_data1.Vp(Nondimensionalize(1500K,CharDim),Nondimensionalize(1e8*Pa,CharDim)) 
Vs_ND    =  PD_data1.Vs(Nondimensionalize(1500K,CharDim),Nondimensionalize(1e8*Pa,CharDim)) 

# redimensionalize and check with value from original structure that did not use non-dimensionalization 
@test   ustrip(Dimensionalize(rho_ND,kg/m^3,CharDim)) ≈ PD_data.Rho(1500,1e8) 
@test   ustrip(Dimensionalize(Vp_ND, km/s,  CharDim)) ≈ PD_data.Vp(1500,1e8) 
@test   ustrip(Dimensionalize(Vs_ND, km/s,  CharDim)) ≈ PD_data.Vs(1500,1e8) 



# Test computation of velocity for the whole computational domain, using arrays 
MatParam    =   Array{MaterialParams, 1}(undef, 3);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=0,
                        SeismicVelocity   = PerpleX_LaMEM_Diagram("test_data/Peridotite.in"));

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=1,
                        SeismicVelocity   = ConstantSeismicVelocity());

MatParam[3] =   SetMaterialParams(Name="UpperCrust", Phase=2,
                        SeismicVelocity   = ConstantSeismicVelocity(Vp=10km/s, Vs=3km/s));

# test computing material properties
Phases              = ones(Int64,400,400)*0;
Phases[:,20:end] .= 1
Phases[:,300:end] .= 2

Vp      = zeros(size(Phases))
Vs      = zeros(size(Phases))
T       =  ones(size(Phases))
P       =  ones(size(Phases))*10

ComputePwaveVelocity!(Vp, Phases, P,T, MatParam)
@test sum(Vp)/400^2 ≈ 8.541562850000005

# test computing material properties when we have PhaseRatios, instead of Phase numbers
PhaseRatio  = zeros(size(Phases)...,length(MatParam));
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz+1)
    PhaseRatio[I] = 1.0  
end

ComputeSwaveVelocity!(Vs, PhaseRatio, P,T, MatParam)
@test sum(Vs)/400^2 ≈ 4.0739837

end