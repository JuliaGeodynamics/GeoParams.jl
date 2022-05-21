using Test
using GeoParams

@testset "DiffusionCreepLaws" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1Pa*s, length=1m);
                
# Define a linear viscous creep law ---------------------------------
x1      =   DiffusionCreep()
@test Value(x1.n) == 1.0
@test Value(x1.p) == -3.0
@test Value(x1.A) == 1.5MPa^-1.0*s^-1*m^(3.0)

x2      =   DiffusionCreep(n=3, p=-3.0)
@test Value(x2.A) == 1.5MPa^-3.0*s^-1*m^(3.0)

# perform a computation with the dislocation creep laws 
    # Calculate EpsII, using a set of pre-defined values
CharDim = GEO_units(length=1000km, viscosity=1e19Pa*s, stress=100MPa, temperature=1000C)
EpsII   = GeoUnit(1.0s^-1.0)
EpsII_nd= nondimensionalize(EpsII,CharDim)
TauII   = GeoUnit(0.3MPa)
TauII_nd= nondimensionalize(TauII,CharDim)
P       = GeoUnit(1.0e9Pa)
P_nd    = nondimensionalize(P,CharDim)
T       = GeoUnit(1400C)
T_nd    = nondimensionalize(T,CharDim)
f       = GeoUnit(1000NoUnits)
f_nd    = nondimensionalize(f,CharDim)
d       = GeoUnit(10mm)
d_nd    = nondimensionalize(d,CharDim)


# compute a pure diffusion creep rheology
p = SetDiffusionCreep("Dry Plagioclase | Bürgmann & Dresen (2008)")

T = 650+273.15;

args = (;T=T )
TauII = 1e6
ε = compute_εII(p, TauII, args)


# test with arrays
τII_array       =   ones(10)*1e6
ε_array         =   similar(τII_array)
T_array         =   ones(size(τII_array))*(650. + 273.15)

args_array = (;T=T_array )

compute_εII!(ε_array, p, τII_array, args_array)
@test ε_array[1] ≈ ε 

# compute when args are scalars
compute_εII!(ε_array, p, τII_array, args)
@test ε_array[1] ≈ ε 


# ===



# dry anorthtite, stress-strainrate curve
p = SetDiffusionCreep("Dry Anorthite | Bürgmann & Dresen (2008)")
EpsII = exp10.(-22:.5:-12);
T     = 650 + 273.15;
gsiz  = 100e-6;
P     = 0.
args  = (T=T, d=gsiz, P=P)
τII_array = zero(EpsII)
compute_τII!(τII_array, p, EpsII, args)
eta_array = τII_array./(2*EpsII)

εII_array = zero(τII_array)
compute_εII!(εII_array, p, τII_array, args)
eta_array1 = τII_array./(2*εII_array)


# matlab script
eII = 1e-22; PPa = 0.0
gsiz        = 100;
TK          = 650+273.15;

logA   = [12.1  12.7]; #Logarithm of pre-exponential factor
npow   = [   1     3]; #Power law exponent
Qact   = [ 460   641]; #Activation energy (KJ)
m_gr0  = [   3     0]; #Grain size Exponent (will convert to negative)
r_fug  = [   0     0]; #Exponent of Fugacity
Vact   = [  24    24]; #Activation Volume cm-3
fugH   = [  1      1]; #Fugacity of water MPa 


# Conversion Factors and constants ---------------------------------------------------
R      = 8.314; #Gas Constant
MPa2Pa = 1e6;   #MPa  -> Pa
cm32m3 = 1e-6;  #cm3  -> m3
J2kJ   = 1e-3;  #Joul -> kJoule

A0     = 10.0.^(logA);
m_gr   = -m_gr0;
PMPa   =  PPa/MPa2Pa;

FG_e   = 1/(2^((npow[i_flow]-1)./npow[i_flow])*3^((npow[i_flow]+1)./(2*npow[i_flow])))
FG_s   = 1/(3^((npow[i_flow]+1)./2));    


iflow = 1;
mu1    =         FG_e.*eII.^(1/npow[i_flow]-1)*A0[i_flow]^(-1.0/npow[i_flow])*gsiz^(-m_gr[i_flow]/npow[i_flow])*fugH[i_flow]^(-r_fug[i_flow]/npow[i_flow])*exp((Qact[i_flow]+PMPa*MPa2Pa.*Vact[i_flow]*cm32m3*J2kJ)/(R*J2kJ*TK*npow[i_flow]));
mu     =        mu1.*MPa2Pa; #In Pa.s




# How the same is implemented in GeoParams
n = 1.0NoUnits                         # power-law exponent
r = 0.0NoUnits                         # exponent of water-fugacity
p = -3.0NoUnits                        # grain size exponent
#A = (10^12.1)MPa^(-1)*μm^3.0*s^(-1)    # material specific rheological parameter
A = 1.258925411794166e-12*m^3/Pa/s

E = 460.0e3J/mol                       # activation energy
V = 24e-6m^3/mol                       # activation Volume
τ = 1e6Pa
P = 0.0Pa
d = 100μm
R = 8.3145J/K/mol
T = 923.15K


ε =  A * τ^n * f^r * d^p * exp( -(E + P*V)/(R*T) )

η = τ/(2ε)  # viscosity

#p = SetDislocationCreep("Dry Anorthite | Rybecki and Dresen (2000)")
#compute_τII!(τII_array, p, EpsII, args)




#Phase   = SetMaterialParams(Name="Viscous Matrix", Phase=1,
#                                     Density   = ConstantDensity(),
#                                     CreepLaws = DiffusionCreep(n=1.0NoUnits, r=1.0NoUnits, p=-3.0NoUnits, A=1.0e6MPa^(-1.0)*s^(-1.0)*µm^(3.0), E=335000.0J/mol, V=4.0e-6m^(3.0)/mol, Apparatus="Invariant"), CharDim = CharDim)
#EpsII.val = ComputeCreepLaw_EpsII(TauII,Phase.CreepLaws[1],c)
# Tested value is Dimensional so need to non-dimensionalize it back to compare
#@test EpsII.val  ≈ nondimensionalize(7.823165072337786e-15s^(-1.0), CharDim) rtol = 1e-12
    # Check that once inverted, we get back the TauII that we used to calculate EpsII
#NewTau = ComputeCreepLaw_TauII(EpsII,Phase.CreepLaws[1],c)
#@show NewTau
#@test NewTau ≈ TauII.val


# Given stress
#@test ComputeCreepLaw_EpsII(1e6Pa, x1, CreepLawParams())==5e-13/s                # dimensional input       

# Given strainrate 
#@test ComputeCreepLaw_TauII(1e-13/s, x1, CreepLawParams())==1e18*2*1e-13Pa       # dimensional input       
# -------------------------------------------------------------------


end
