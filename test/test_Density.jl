using Test
using GeoParams
@testset "Density.jl" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
                
# Define a linear viscous creep law
x1      =   ConstantDensity(ρ=2900kg/m^3)
@test x1.ρ.val == 2900

x1 = nondimensionalize(x1,CharUnits_GEO)
@test x1.ρ.val ≈ 2.9e-16

x2      =   PT_Density()
@test x2.α.val==3e-5
@test x2.ρ0.val==2900.0

x2 = nondimensionalize(x2,CharUnits_GEO)
@test x2.T0.val≈0.21454659702313156

# Compute with density while specifying P & T (not used in case of )
@test compute_density(x2,1.0,1.0) ≈ 2.8419999999999996e-16
@test compute_density(x1,1.0,1.0) ≈ 2.9e-16

# test to allocations
rho = [0.0]
P   = 1.0;
T   = 1.0

x   = ConstantDensity()
num_alloc = @allocated compute_density!(rho, x, P, T)
num_alloc = @allocated compute_density!(rho, x, P, T)
@show num_alloc
@test num_alloc == 0

# This does NOT allocate if I test this with @btime;
#   yet it does while running the test here
x   = PT_Density()
num_alloc = @allocated compute_density!(rho, x, P, T)
num_alloc = @allocated compute_density!(rho, x, P, T)
@show num_alloc
@test num_alloc <= 32

# This does NOT allocate if I test this with @btime;
#   yet it does while running the test here
x   = Compressible_Density()
num_alloc = @allocated compute_density!(rho, x, P, T)
num_alloc = @allocated compute_density!(rho, x, P, T)
@show num_alloc
@test num_alloc <=32

# Read Phase diagram interpolation object
fname   =   "./test_data/Peridotite.in"
PD_data =   PerpleX_LaMEM_Diagram(fname);
@test PD_data.meltFrac(1500,1e7) ≈ 0.24368492372485706
@test PD_data.Rho(1500,1e7) ≈ 3042.836820256982
@test PD_data.meltRho(1500,1e7) ≈ 2662.227167592414
@test PD_data.rockRho(1500,1e7) ≈ 3165.467673917775

@test compute_density(PD_data, 1e7, 1500.0) ≈ 3042.836820256982

# Do the same but non-dimensionalize the result
CharDim  =  GEO_units();
PD_data1 =  PerpleX_LaMEM_Diagram(fname, CharDim=CharDim);

rho_ND   =  PD_data1.Rho(nondimensionalize(1500.0K,CharDim), nondimensionalize(1e8*Pa,CharDim)) 
Vp_ND    =  PD_data1.Vp( nondimensionalize(1500.0K,CharDim), nondimensionalize(1e8*Pa,CharDim)) 
Vs_ND    =  PD_data1.Vs( nondimensionalize(1500.0K,CharDim), nondimensionalize(1e8*Pa,CharDim)) 

# redimensionalize and check with value from original structure that did not use non-dimensionalization 
@test   ustrip(dimensionalize(rho_ND,kg/m^3,CharDim)) ≈ PD_data.Rho(1500.0,1e8) 
@test   ustrip(dimensionalize(Vp_ND, km/s,  CharDim)) ≈ PD_data.Vp( 1500.0,1e8) 
@test   ustrip(dimensionalize(Vs_ND, km/s,  CharDim)) ≈ PD_data.Vs( 1500.0,1e8) 


# Test computation of density for the whole computational domain, using arrays 
MatParam    =   Vector{AbstractMaterialParamsStruct}(undef, 4)

MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=0,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PerpleX_LaMEM_Diagram("test_data/Peridotite.in"));

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=1,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = ConstantDensity(ρ=2900kg/m^3));

MatParam[3] =   SetMaterialParams(Name="UpperCrust", Phase=2,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PT_Density());
MatParam[4] =   SetMaterialParams(Name="UpperCrust", Phase=3,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = Compressible_Density());

Mat_tup = Tuple(MatParam);  # create a tuple to avoid allocations


MatParam1 = Vector{AbstractMaterialParamsStruct}(undef, 4)
MatParam1[1] = SetMaterialParams(Name="Crust", Phase=0,
                            CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pas)),
                            Density   = ConstantDensity(ρ=2900kg/m^3))
MatParam1[2] = SetMaterialParams(Name="Lower Crust", Phase=1,
                            CreepLaws= (PowerlawViscous(n=5.), LinearViscous(η=1e21Pas)),
                            Density  = Compressible_Density(ρ0=3000kg/m^3))
MatParam1[3] = SetMaterialParams(Name="Lower Crust", Phase=2,
                            CreepLaws= LinearViscous(η=1e21Pas),
                            Density  = ConstantDensity())
MatParam1[4] = SetMaterialParams(Name="Lower Crust", Phase=3,
                            CreepLaws= LinearViscous(η=1e21Pas),
                            Density  = ConstantDensity())
Mat_tup1 = Tuple(MatParam1)

# test computing material properties
Phases              = ones(Int64,400,400)*0;
Phases[:,20:end]   .= 1
Phases[:,200:end]  .= 2
Phases[:,300:end]  .= 3

#Phases .= 2;
rho     = zeros(size(Phases))
T       =  ones(size(Phases))
P       =  ones(size(Phases))*10

# Test computing density when Mat_tup1 is provided as a tuple
compute_density!(rho, Mat_tup1, Phases, P=P,T=T) 
num_alloc = @allocated compute_density!(rho, Mat_tup1, Phases, P=P,T=T)   #      287.416 μs (0 allocations: 0 bytes)
@test sum(rho)/400^2 ≈ 2945.000013499999
@test num_alloc == 32

# If we employ a phase diagram many allocations occur:
compute_density!(rho, Mat_tup, Phases, P=P,T=T)   #        37.189 ms (1439489 allocations: 26.85 MiB)     - the allocations are from the phase diagram
@test sum(rho)/400^2 ≈ 2920.6151145725003

# test computing material properties when we have PhaseRatios, instead of Phase numbers
PhaseRatio  = zeros(size(Phases)...,length(Mat_tup1));
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz+1)
    PhaseRatio[I] = 1.0  
end

compute_density!(rho, Mat_tup1, PhaseRatio, P=P,T=T)
num_alloc = @allocated compute_density!(rho, Mat_tup1, PhaseRatio, P=P,T=T) #   136.776 μs (0 allocations: 0 bytes)
# using BenchmarkTools
# @btime compute_density!($rho, $Mat_tup1, $PhaseRatio, P=$P,T=$T)
@test sum(rho)/400^2 ≈ 2945.000013499999
@test num_alloc == 32           # for some reason this does indicate allocations but @btime does not

# Test calling the routine with only pressure as input. 
# This is ok for Mat_tup1, as it only has constant & P-dependent densities.
# Note, however, that if you have P & T dependent densities and do this it will use 0 as defualt value for T 
compute_density!(rho, Mat_tup1, PhaseRatio, P=P)
@test sum(rho)/400^2 ≈ 2945.000013499999


end