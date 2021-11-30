export  DiffusionCreep,
        SetDiffusionCreep

#=----Diffusion Creep---
Defines diffusion creep law parameters

n is the power-law exponent
r is the water-fugacity exponent
p is the positiv defined grain size exponent (value is set negative in the calculation of σ and ε)
A is a material specific rheological parameter
E is the activation energy
V is the activation volume
R is the universal gas constant
Apparatus defines the appartus type that shall be recreated (Axial Compression, Simple Shear, Invariant)
=#

@with_kw_noshow mutable struct DiffusionCreep <: AbstractCreepLaw
    equation::LaTeXString = "to be updated"
    n::GeoUnit            = 1.0NoUnits         # power-law exponent
    r::GeoUnit            = 0.0NoUnits         # exponent of water-fugacity
    p::GeoUnit            = 3.0NoUnits         # grain size exponent
    A::GeoUnit            = 1.5MPa^(-n-r)/s    # material specific rheological parameter
    E::GeoUnit            = 500kJ/mol          # activation energy
    V::GeoUnit            = 6e-6m^3/mol        # activation volume
    R::GeoUnit            = 8.314J/mol/K       # universal gas constant
    Apparatus             = "AxialCompression" # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    Comment::String       = ""                 # comment when implementing new creep laws
    BibTex_Reference      = ""                 # BibTex reference
end



function ComputeDiffCreepLaw_EpsII(TauII, a::DiffusionCreep, p::CreepLawVariables)
    @unpack n         = a
    @unpack r         = a
    @unpack p         = a
    @unpack A         = a
    @unpack E         = a
    @unpack V         = a
    @unpack R         = a
    @unpack Apparatus = a
    @unpack P         = p
    @unpack T         = p
    @unpack d         = p
    @unpack f         = p
    if Apparatus == "AxialCompression"
        FT = sqrt(3.0)NoUnits
        FE = 2.0/sqrt(3.0)NoUnits
    elseif Apparatus == "SimpleShear"
        FT = 2.0NoUnits
        FE = 2.0NoUnits
    elseif Apparatus == "Invariant"
        FT = 1.0NoUnits
        FE = 1.0NoUnits
    end
    return A.val*(TauII.val*FT)^n.val*d.val^p.val*f.val^r.val*exp(-(E.val+P.val*V.val)/(R.val*T.val))/FE
end

function ComputeDiffCreepLaw_TauII(EpsII, a::DiffusionCreep, p::CreepLawVariables)
    @unpack n = a
    @unpack r = a
    @unpack p = a
    @unpack A = a
    @unpack E = a
    @unpack V = a
    @unpack R = a
    @unpack Apparatus = a
    @unpack P = a
    @unpack T = a
    @unpack d = a
    @unpack f = a
    if Apparatus == "AxialCompression"
        FT = sqrt(3.0)NoUnits               # relation between differential stress recorded by apparatus and TauII
        FE = 2.0/sqrt(3.0)NoUnits           # relation between gamma recorded by apparatus and EpsII
    elseif Apparatus == "SimpleShear"
        FT = 2.0NoUnits                     # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = 2.0NoUnits
    elseif Apparatus == "Invariant"
        FT = 1.0NoUnits
        FE = 1.0NoUnits
    end
    return A.val^(-1/n.val)*(EpsII.val*FE)^(1/n.val)*d.val^(-p.val/n.val)*f.val^(-r.val/n.val)*exp((E.val+P.val*V.val)/(n.val*R.val*T.val))/FT;
end


# Print info 
function show(io::IO, g::DislocationCreep)  
    print(io, "DiffusionCreep: n=$(g.n.val), r=$(g.r.val), p=$(g.p.val), A=$(g.A.val), E=$(g.E.val), V=$(g.V.val), Apparatus=$(g.Apparatus)" )  
end

# predefined diffusion creep laws are to be added in the dictionary as it is done for dislocation creep laws (see 'DislocationCreep.jl')!

SetDiffusionCreep = Dict([])
