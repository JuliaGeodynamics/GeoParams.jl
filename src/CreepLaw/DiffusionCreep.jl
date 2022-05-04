using ..MaterialParameters: MaterialParamsInfo
import GeoParams.param_info


export  DiffusionCreep,
        SetDiffusionCreep

const AxialCompression, SimpleShear, Invariant = 1,2,3

#=----Diffusion Creep---
Defines diffusion creep law parameters

n is the power-law exponent
r is the water-fugacity exponent
p is the negative defined grain size exponent (value is set negative in the calculation of σ and ε)
A is a   material specific rheological parameter
E is the activation energy
V is the activation volume
R is the universal gas constant
Apparatus defines the appartus type that shall be recreated (Axial Compression, Simple Shear, Invariant)
=#

@with_kw_noshow mutable struct DiffusionCreep{T,N,U1,U2,U3,U4,U5} <: AbstractCreepLaw{T}
    Name::NTuple{N,Char}        = ""
    n::GeoUnit{T,U1}            = 1.0NoUnits         # power-law exponent
    r::GeoUnit{T,U1}            = 0.0NoUnits         # exponent of water-fugacity
    p::GeoUnit{T,U1}            = -3.0NoUnits         # grain size exponent
    A::GeoUnit{T,U2}            = 1.5MPa^(-n-r)*s^(-1)*m^(-p)    # material specific rheological parameter
    E::GeoUnit{T,U3}            = 500kJ/mol          # activation energy
    V::GeoUnit{T,U4}            = 6e-6m^3/mol        # activation volume
    R::GeoUnit{T,U5}            = 8.314J/mol/K       # universal gas constant
    Apparatus::Int32            = AxialCompression # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
end


DiffusionCreep(args...) = DiffusionCreep(NTuple{length(args[1]), Char}(collect.(args[1])), convert.(GeoUnit,args[2:end-1])..., args[end])

function param_info(s::DiffusionCreep)
    name = String(collect(s.Name))
    eq = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"
    if name == "" 
        return MaterialParamsInfo(Equation=eq)
    end
    inf = DiffusionCreep_info[name][2]
    return MaterialParamsInfo(Equation=eq, Comment=inf.Comment, BibTex_Reference=inf.BibTex_Reference)
end

# Calculation routines for linear viscous rheologies
# All inputs must be non-dimensionalized (or converted to consitent units) GeoUnits
function computeCreepLaw_EpsII(TauII, a::DiffusionCreep; P, T, f, d, kwargs...)
    @unpack_val n,r,p,A,E,V,R = a
    
    FT, FE = CorrectionFactor(a);    
   
    return A*(TauII*FT)^n*f^r*d^p*exp(-(E + P*V)/(R*T))/FE; 
end

# EpsII .= A.*(TauII.*FT).^n.*f.^r.*d.^p.*exp.(-(E.+P.*V)./(R.*T))./FE; Once we have a 
# All inputs must be non-dimensionalized (or converted to consistent units) GeoUnits
function computeCreepLaw_TauII(EpsII, a::DiffusionCreep; P, T, f, d, kwargs...)
    @unpack_val n,r,p,A,E,V,R = a

    FT, FE = CorrectionFactor(a);    

    return A^(-1/n)*(EpsII*FE)^(1/n)*f^(-r/n)*d^(-p/n)*exp((E + P*V)/(n * R*T))/FT;
end



# Print info 
function show(io::IO, g::DiffusionCreep)  
    print(io, "DiffusionCreep: Name = $(String(collect(g.Name))), n=$(g.n.val), r=$(g.r.val), p=$(g.p.val), A=$(g.A.val), E=$(g.E.val), V=$(g.V.val), Apparatus=$(g.Apparatus)" )  
end


# This computes correction factors to go from experimental data to tensor format
# A nice discussion 
function CorrectionFactor(a::DiffusionCreep{_T}) where {_T}
    if a.Apparatus == AxialCompression
        FT = sqrt(one(_T)*3)               # relation between differential stress recorded by apparatus and TauII
        FE = one(_T)*2/sqrt(one(_T)*3)     # relation between gamma recorded by apparatus and EpsII
    elseif a.Apparatus == SimpleShear
        FT = one(_T)*2                     # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = one(_T)*2 
    end
    return FT,FE
end

# Add pre-defined creep laws 
"""
    SetDiffusionCreep["Name of Diffusion Creep"]
This is a dictionary with pre-defined creep laws    
"""
SetDiffusionCreep(name::String) = DiffusionCreep_info[name][1]

# predefined diffusion creep laws are to be added in the dictionary as it is done for dislocation creep laws (see 'DislocationCreep.jl')!

DiffusionCreep_info = Dict([])
