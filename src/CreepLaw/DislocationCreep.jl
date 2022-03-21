# This implements viscous creep laws and routines to compute with them
#
# Note that various simple creep laws are defined in this file; 
# more complex ones (such as DislocationCreep) are in separate files 
# in the same directory
#
# In case you want to add new creep laws, have a look at how the ones
# here are implemented. Please add tests as well!

using ..MaterialParameters: MaterialParamsInfo
import GeoParams.param_info

export  DislocationCreep,
        SetDislocationCreep

const AxialCompression, SimpleShear, Invariant = 1,2,3




# Dislocation Creep ------------------------------------------------
"""
    DislocationCreep(n = 1.0NoUnits, r = 0.00.0NoUnits, A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = "AxialCompression" )
    
Defines the flow law parameter of a dislocation creep law 

The (isotropic) dislocation creep law, as used by experimtalists, is given by  
```math  
     \\dot{\\gamma} = A \\sigma_\\mathrm{d}^n f_\\mathrm{H2O}^r \\exp(-\\frac{E+PV}{RT})
```

where ``n`` is the power law exponent,  
``r`` is the exponent of fugacity dependence, 
``A`` is a pre-exponential factor [MPa^(n+r)] (if manually defined, n and r must be either pre-defined or substituted),  
``E`` is the activation energy [kJ/mol], ``V`` is the activation volume [m^3/mol]. ``\\dot{\\gamma}`` is the ordinary strain rate [1/s], 
and ``\\sigma_\\mathrm{d}`` is the differential stress which are converted into second invariants using the apparatus type that can be
either "AxialCompression", "SimpleShear" or "Invariant".
If the flow law paramters are already given as a function of second invariants, choose apparatus = "Invariant"

```julia-repl 
julia> x2      =   DislocationCreep(n=3)
DislocationCreep: n=3, r=0.0, A=1.5 MPa^-3 s^-1, E=476.0 kJ mol^-1, V=6.0e-6 m^3 mol^-1, Apparatus=AxialCompression
```
"""
@with_kw_noshow struct DislocationCreep{T,N,U1,U2,U3,U4,U5} <: AbstractCreepLaw{T}
    Name::NTuple{N,Char}    =   ""               # The name is encoded as a NTuple{Char} to make it isbits    
    n::GeoUnit{T,U1}        = 1.0NoUnits         # power-law exponent
    r::GeoUnit{T,U1}        = 0.0NoUnits         # exponent of water-fugacity dependence
    A::GeoUnit{T,U2}        = 1.5MPa^(-n-r)/s    # pre-exponential factor
    E::GeoUnit{T,U3}        = 476.0kJ/mol        # activation energy
    V::GeoUnit{T,U4}        = 6e-6m^3/mol        # activation volume
    R::GeoUnit{T,U5}        = 8.314J/mol/K       # Universal gas constant
    Apparatus::Int32        = AxialCompression   # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
end
DislocationCreep(args...) = DislocationCreep(NTuple{length(args[1]), Char}(collect.(args[1])), convert.(GeoUnit,args[2:end-1])..., args[end])

function param_info(s::DislocationCreep)
    name = String(collect(s.Name))
    eq = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"
    if name == "" 
        return MaterialParamsInfo(Equation=eq)
    end
    inf = DislocationCreep_info[name][2]
    return MaterialParamsInfo(Equation=eq, Comment=inf.Comment, BibTex_Reference=inf.BibTex_Reference)
end

# Calculation routines for linear viscous rheologies
# All inputs must be non-dimensionalized (or converted to consitent units) GeoUnits
function computeCreepLaw_EpsII(TauII, a::DislocationCreep, p::CreepLawVariables)
    @unpack_val n,r,A,E,V,R = a
    @unpack_val P,T,f       = p
    
    FT, FE = CorrectionFactor(a)
   
    return A*(TauII*FT)^n*f^r*exp(-(E + P*V)/(R*T))/FE
end

function computeCreepLaw_EpsII(TauII, a::DislocationCreep, P::_R, T::_R, f::_R) where _R<:Real
    @unpack_val n,r,A,E,V,R = a
    
    FT, FE = CorrectionFactor(a);    
   
    return A*(TauII*FT)^n*f^r*exp(-(E + P*V)/(R*T))/FE; 
end

# EpsII .= A.*(TauII.*FT).^n.*f.^r.*exp.(-(E.+P.*V)./(R.*T))./FE; Once we have a 
# All inputs must be non-dimensionalized (or converted to consistent units) GeoUnits
function computeCreepLaw_TauII(EpsII, a::DislocationCreep, p::CreepLawVariables)
    @unpack_val n,r,A,E,V,R = a
    @unpack_val P,T,f       = p

    FT, FE = CorrectionFactor(a)    

    return A^(-1/n)*(EpsII*FE)^(1/n)*f^(-r/n)*exp((E + P*V)/(n * R*T))/FT;
end


# EpsII .= A.*(TauII.*FT).^n.*f.^r.*exp.(-(E.+P.*V)./(R.*T))./FE; Once we have a 
# All inputs must be non-dimensionalized (or converted to consistent units) GeoUnits
function computeCreepLaw_TauII(EpsII, a::DislocationCreep, P::_R, T::_R, f::_R) where _R<:Real
    @unpack_val n,r,A,E,V,R = a

    FT, FE = CorrectionFactor(a);    

    return A^(-1/n)*(EpsII*FE)^(1/n)*f^(-r/n)*exp((E + P*V)/(n * R*T))/FT
end


# Print info 
function show(io::IO, g::DislocationCreep)  
    print(io, "DislocationCreep: Name = $(String(collect(g.Name))), n=$(g.n.val), r=$(g.r.val), A=$(g.A.val), E=$(g.E.val), V=$(g.V.val), Apparatus=$(g.Apparatus)" )  
end
#-------------------------------------------------------------------------

# This computes correction factors to go from experimental data to tensor format
# A nice discussion 
function CorrectionFactor(a::DislocationCreep{_T}) where {_T}

    FT = one(_T) 
    FE = one(_T)
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
    SetDislocationCreep["Name of Dislocation Creep"]

This is a dictionary with pre-defined creep laws    
"""
SetDislocationCreep(name::String) = DislocationCreep_info[name][1]

DislocationCreep_info = Dict([

# Olivine rheology 
("Dry Olivine | Hirth & Kohlstedt (2003)", 
# after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
# Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters
    (DislocationCreep(
        Name = "Dry Olivine | Hirth & Kohlstedt (2003)",
        n = 3.5NoUnits,
        A = 1.1e5MPa^(-3.5)/s, 
        E = 530.0kJ/mol,
        V = 15e-6m^3/mol,
        Apparatus = AxialCompression,
        r = 0.0NoUnits),
    MaterialParamsInfo(Comment = "Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
        BibTex_Reference = parse_bibtex("""
            @incollection{eiler_rheology_2003,
            address = {Washington, D. C.},
            title = {Rheology of the upper mantle and the mantle wedge: {A} view from the experimentalists},
            volume = {138},
            isbn = {978-0-87590-997-4},
            shorttitle = {Rheology of the upper mantle and the mantle wedge},
            url = {http://www.agu.org/books/gm/v138/138GM06/138GM06.shtml},
            language = {en},
            urldate = {2019-10-09},
            booktitle = {Geophysical {Monograph} {Series}},
            publisher = {American Geophysical Union},
            author = {Hirth, Greg and Kohlstedt, David},
            editor = {Eiler, John},
            year = {2003},
            doi = {10.1029/138GM06},
            pages = {83--105},
            }
        """))
    )
)

# Olivine rheology 
("Wet Olivine | Hirth & Kohlstedt (2003)", 
    # After Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    #   Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
    #  Note that this assumes C_OH=1000
    (DislocationCreep(
        Name = "Wet Olivine | Hirth & Kohlstedt (2003)",
        n = 3.5NoUnits,
        A = 90MPa^(-3.5)/s, 
        E = 480kJ/mol,
        V = 11e-6m^3/mol,
        r   = 1.2NoUnits,
        Apparatus = AxialCompression),
    MaterialParamsInfo(Comment = "Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
        BibTex_Reference = parse_bibtex("""
            @incollection{HirthKohlstedt_OlivineRheology_2003,
            address = {Washington, D. C.},
            title = {Rheology of the upper mantle and the mantle wedge: {A} view from the experimentalists},
            volume = {138},
            isbn = {978-0-87590-997-4},
            shorttitle = {Rheology of the upper mantle and the mantle wedge},
            url = {http://www.agu.org/books/gm/v138/138GM06/138GM06.shtml},
            language = {en},
            urldate = {2019-10-09},
            booktitle = {Geophysical {Monograph} {Series}},
            publisher = {American Geophysical Union},
            author = {Hirth, Greg and Kohlstedt, David},
            editor = {Eiler, John},
            year = {2003},
            doi = {10.1029/138GM06},
            pages = {83--105},
            }
        """))
    )
)



]); # end of setting pre-defined creep laws
