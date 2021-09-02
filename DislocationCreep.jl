module CreepLaw

# This implements viscous creep laws and routines to compute with them
#
# Note that various simple creep laws are defined in this file; 
# more complex ones (such as DislocationCreep) are in separate files 
# in the same directory
#
# In case you want to add new creep laws, have a look at how the ones
# here are implemented. Please add tests as well!

using Base: Float64, String
using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import Base.show

export  DislocationCreep,
        SetDislocationCreep

# Dislocation Creep ------------------------------------------------
"""
    DislocationCreep(n = 1.0NoUnits, r = 0.00.0NoUnits, A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = "AxialCompression" )
    
Defines the flow law parameter of a dislocation creep law 

The (isotropic) dislocation creep law, as used by experimtalists, is given by  
```math  
     \\dot{\\gamma} = A \\sigma_\\mathrm{d}^n f_\\mathrm{H2O}^r \\exp(-\\frac{E+PV}{RT})
```

where ``n`` is the power law exponent, ``r`` is the exponent of fugacity dependence, ``A`` is a pre-exponential factor [MPa^(n+r)] (if manually defined, n and r must be either pre-defined or substituted),  
``E`` is the activation energy [kJ/mol], ``V`` is the activation volume [m^3/mol]. ``\\dot{\\gamma}`` is the ordinary strain rate [1/s], 
and ``\\sigma_\\mathrm{d}`` is the differential stress which are converted into second invariants using the apparatus type that can be
either "AxialCompression", "SimpleShear" or "Unknown".
If the flow law paramters are given as a function of second invariants, choose apparatus = "Unknown"
"""
@with_kw_noshow mutable struct DislocationCreep <: AbstractCreepLaw
    equation::LaTeXString   =   L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}" 
    n::GeoUnit        = 1.0NoUnits         # power-law exponent
    r::GeoUnit        = 0.0NoUnits         # exponent of water-fugacity dependence
    A::GeoUnit        = 1.5MPa^(-n-r)/s    # pre-exponential factor
    E::GeoUnit        = 476.0kJ/mol        # activation energy
    V::GeoUnit        = 6e-6m^3/mol        # activation volume
    R::GeoUnit        = 8.314J/mol/K       # Universal gas constant
    Apparatus::string = "AxialCompression" # type of experimental apparatus, either AxialCompression, SimpleShear or Unknown
end

# Calculation routines for dislocation creep rheologies
function ComputeCreepLaw_EpsII(EpsII, a::DislocationCreep, p::CreepLawVariables, CharDim)
    @unpack n           = a
    @unpack r           = a
    @unpack A           = a
    @unpack E           = a
    @unpack V           = a
    @unpack Apparatus   = a
    if typeof(CharDim) <: GeoUnits
        n =   Nondimensionalize!(n, CharDim);
        r =   Nondimensionalize!(r, CharDim);
        A =   Nondimensionalize!(A, CharDim);
        E =   Nondimensionalize!(E, CharDim);
        V =   Nondimensionalize!(V, CharDim);
    else
        error("CharDim should be of type GeoUnits")
    end
    if Apparatus == "AxialCompression"
        FT = sqrt(3.0)               # relation between differential stress recorded by apparatus and TauII
        FE = 2.0/sqrt(3.0)           # relation between gamma recorded by apparatus and EpsII
    elseif Apparatus == "SimpleShear"
        FT = 2.0                     # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = 2.0
    elseif Apparatus == "Unknown"
        FT = 1.0
        FE = 1.0
    end
    return EpsII = A*(TauII*FT)^n*f^r*exp(-(E+P*V)/(R*T))/FE;
end

# EpsII .= A.*(TauII.*FT).^n.*f.^r.*exp.(-(E.+P.*V)./(R.*T))./FE; Once we have a 

function ComputeCreepLaw_TauII(TauII, a::DislocationCreep, p::CreepLawVariables, CharDim)
    @unpack n           = a
    @unpack r           = a
    @unpack A           = a
    @unpack E           = a
    @unpack V           = a
    @unpack Apparatus   = a
    if typeof(CharDim) <: GeoUnits
        n =   Nondimensionalize!(n, CharDim);
        r =   Nondimensionalize!(r, CharDim);
        A =   Nondimensionalize!(A, CharDim);
        E =   Nondimensionalize!(E, CharDim);
        V =   Nondimensionalize!(V, CharDim);
    else
        error("CharDim should be of type GeoUnits")
    end
    @unpack Apparatus   = a
    if Apparatus == "AxialCompression"
        FT = sqrt(3.0)               # relation between differential stress recorded by apparatus and TauII
        FE = 2.0/sqrt(3.0)           # relation between gamma recorded by apparatus and EpsII
    elseif Apparatus == "SimpleShear"
        FT = 2.0                     # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = 2.0
    elseif Apparatus == "Unknown"
        FT = 1.0
        FE = 1.0
    end
    return TauII = (A*EpsII*FE)^(-1/n)*f^(-r/n)*exp((E+P*V)/(n*R*T))/FT;
end


# Print info 
function show(io::IO, g::DislocationCreep)  
    print(io, "DislocationCreep: n=$(g.n.val), r=$(g.r.val), A=$(g.A.val), E=$(g.E.val), V=$(g.V.val), apparatus=$(g.Apparatus.val)" )  
end
#-------------------------------------------------------------------------

# add routine SetDislocationCreep, to look up flow law params for any hard coded flow laws
function SetDislocationCreep(CreepLawName::String)

AvailableCreepLaws = [  "Dry Olivine Hirth and Kohlstedt 2003",
                        "Wet Olivine Hirth and Kohlstedt 2003"
                     ]


if      CreepLawName == "Dry Olivine Hirth and Kohlstedt 2003"
   A = 1e8MPa 
   n = 3.04

elseif  CreepLawName == "Wet Olivine Hirth and Kohlstedt 2003"
   


else
    error("Unknown dislocation creep law")
end



DislocationCreep()


end