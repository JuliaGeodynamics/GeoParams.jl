using Test
using GeoParams, LinearAlgebra
@testset "ZirconSaturation.jl" begin

#Make sure structure is isbits
x = Tierney()
@test isbits(x)

# This tests the various melting parameterizations
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=10km);
       
T       =   Vector(800:20:1100.)*K 
T_nd    =   Float64.(T/CharUnits_GEO.Temperature)

# Tierney et al. parameterization [in ND numbers, which is anyways the typical use case]
p        =  Tierney()
Fzrs_dim =  zeros(size(T))
compute_zirconsaturation!(Fzrs_dim, p,  zeros(size(T)), ustrip.(T))

Fzrs_dim1  = zeros(size(Fzrs_dim))
compute_zirconsaturation!(Fzrs_dim1, p,  zeros(size(T)), ustrip.(T)) # in-place routine

p_nd     =  p
p_nd     =  nondimensionalize(p_nd, CharUnits_GEO)
Fzrs_nd  =  zeros(size(T))
compute_zirconsaturation!(Fzrs_nd, p_nd,  zeros(size(T)), T_nd)

# Do this computation manually
Fzrs_anal   =   1.62 .- 1.8e4*exp.(-1e4./ ustrip.(T)); 
Fzrs_anal[T.<Value(p.T_s)] .= 1.
Fzrs_anal[T.>Value(p.T_zrs)] .= 0.

@test sum(Fzrs_dim  - Fzrs_anal)   < 1e-12
@test sum(Fzrs_dim1 - Fzrs_anal)   < 1e-12
@test sum(Fzrs_nd - Fzrs_anal)    < 1e-12


# Test computation of Zircon saturation for the whole computational domain, using arrays 
MatParam    =   Array{MaterialParams, 1}(undef, 3);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                            Melting   = MeltingParam_Caricchi(),
                    ZirconSaturation  = Tierney());

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                                Melting   = MeltingParam_Caricchi(),
                        ZirconSaturation  = Tierney());

# No melting parameterization for this phase
MatParam[3] =   SetMaterialParams(Name="UpperCrust", Phase=3,
                    ZirconSaturation  = Tierney(),
                        Density   = PT_Density());

Mat_tup = Tuple(MatParam)

# test computing material properties
n = 100;
Phases              = ones(Int64,n,n);
Phases[:,:,20:end] .= 2;
Phases[:,:,80:end] .= 3;

Fzrs =  zeros(size(Phases));
ϕ    = zeros(size(Phases));
T    =  ones(size(Phases))*1000.;
P    =  ones(size(Phases));

compute_zirconsaturation!(Fzrs, Mat_tup, Phases, P, T) #allocations coming from computing meltfraction using PhaseDiagram_LookupTable
#compute_meltfraction!(ϕ, Mat_tup, Phases, P,T) 

@test sum(Fzrs)/n^2 ≈ 0.8028012642752728


# test computing material properties when we have PhaseRatios, instead of Phase numbers
PhaseRatio  = zeros(n,n,3);
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz)
    PhaseRatio[I] = 1.0  
end

compute_zirconsaturation!(Fzrs, Mat_tup, PhaseRatio, P, T)
@test sum(Fzrs)/n^2 ≈  0.8028012642752728


end

