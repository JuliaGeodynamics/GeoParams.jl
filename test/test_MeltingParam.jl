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
Nondimensionalize!(p_nd, CharUnits_GEO)
phi_nd   =  ComputeMeltingParam(0,T_nd, p_nd)

# Do this computation manually, using the actual expression of Caricchi
T_C         =   250:100:1250    # in celcius
Phi_solid   =   1.0 .- 1.0./(1.0 .+ exp.((800.0 .- T_C)./23.0)); 
Phi_anal    =   1.0 .- Phi_solid

@test sum(phi_dim  - Phi_anal)   < 1e-12
@test sum(phi_dim1 - Phi_anal)   < 1e-12
@test sum(phi_nd - Phi_anal)    < 1e-12

end

