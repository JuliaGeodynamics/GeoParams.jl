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
PD_data =   PerpleX_LaMEM_Diagram(fname);
@test PD_data.meltFrac(1500,1e7) ≈ 0.24368492372485706
@test PD_data.Rho(1500,1e7) ≈ 3042.836820256982
@test PD_data.meltRho(1500,1e7) ≈ 2662.227167592414
@test PD_data.rockRho(1500,1e7) ≈ 3165.467673917775


@test ComputeDensity(1e7, 1500, PD_data) ≈ 3042.836820256982

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



# Test computation of density for the whole computational domain, using arrays 
MatParam    =   Array{MaterialParams, 1}(undef, 3);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PerpleX_LaMEM_Diagram("test_data/Peridotite.in"));

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = ConstantDensity(ρ=2900kg/m^3));

MatParam[3] =   SetMaterialParams(Name="UpperCrust", Phase=3,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PT_Density());

# test computing material properties
Phases              = ones(Int64,400,400);
Phases[:,20:end] .= 2
Phases[:,300:end] .= 3

rho     = zeros(size(Phases))
T       =  ones(size(Phases))
P       =  ones(size(Phases))*10

ComputeDensity!(rho, Phases, P,T, MatParam)
@test sum(rho)/400^2 ≈ 2920.6148898225


end