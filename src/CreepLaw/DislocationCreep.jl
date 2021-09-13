# This implements viscous creep laws and routines to compute with them
#
# Note that various simple creep laws are defined in this file; 
# more complex ones (such as DislocationCreep) are in separate files 
# in the same directory
#
# In case you want to add new creep laws, have a look at how the ones
# here are implemented. Please add tests as well!

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
@with_kw_noshow mutable struct DislocationCreep <: AbstractCreepLaw
    equation::LaTeXString   =   L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}" # TO BE UPDATED
    n::GeoUnit        = 1.0NoUnits         # power-law exponent
    r::GeoUnit        = 0.0NoUnits         # exponent of water-fugacity dependence
    A::GeoUnit        = 1.5MPa^(-n-r)/s    # pre-exponential factor
    E::GeoUnit        = 476.0kJ/mol        # activation energy
    V::GeoUnit        = 6e-6m^3/mol        # activation volume
    R::GeoUnit        = 8.314J/mol/K       # Universal gas constant
    Apparatus         = "AxialCompression" # type of experimental apparatus, either AxialCompression, SimpleShear or Invariant
    Comment::String   = ""                 # Some remarks you want to add about this creep law implementation
    BibTex_Reference  = ""                 # BibTeX reference
end

# Calculation routines for linear viscous rheologies
# All inputs must be non-dimensionalized (or converted to consitent units) GeoUnits
function ComputeCreepLaw_EpsII(TauII, a::DislocationCreep, p::CreepLawVariables)
    @unpack n           = a
    @unpack r           = a
    @unpack A           = a
    @unpack E           = a
    @unpack V           = a
    @unpack R           = a
    @unpack Apparatus   = a
    @unpack P           = p
    @unpack T           = p
    @unpack f           = p
    if Apparatus == "AxialCompression"
        FT = sqrt(3.0)NoUnits               # relation between differential stress recorded by apparatus and TauII
        FE = 2.0/sqrt(3.0)NoUnits            # relation between gamma recorded by apparatus and EpsII
    elseif Apparatus == "SimpleShear"
        FT = 2.0NoUnits                      # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = 2.0NoUnits 
    elseif Apparatus == "Invariant"
        FT = 1.0NoUnits 
        FE = 1.0NoUnits 
    end
    return A.val*(TauII.val*FT)^n.val*f.val^r.val*exp(-(E.val+P.val*V.val)/(R.val*T.val))/FE; 
end

# EpsII .= A.*(TauII.*FT).^n.*f.^r.*exp.(-(E.+P.*V)./(R.*T))./FE; Once we have a 
# All inputs must be non-dimensionalized (or converted to consitent units) GeoUnits
function ComputeCreepLaw_TauII(EpsII, a::DislocationCreep, p::CreepLawVariables)
    @unpack n           = a
    @unpack r           = a
    @unpack A           = a
    @unpack E           = a
    @unpack V           = a
    @unpack R           = a
    @unpack Apparatus   = a
    @unpack P           = p
    @unpack T           = p
    @unpack f           = p
    if Apparatus == "AxialCompression"
        FT = sqrt(3.0)               # relation between differential stress recorded by apparatus and TauII
        FE = 2.0/sqrt(3.0)           # relation between gamma recorded by apparatus and EpsII
    elseif Apparatus == "SimpleShear"
        FT = 2.0                     # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = 2.0
    elseif Apparatus == "Invariant"
        FT = 1.0
        FE = 1.0
    end
    return A.val^(-1/n.val)*(EpsII.val*FE)^(1/n.val)*f.val^(-r.val/n.val)*exp((E.val+P.val*V.val)/(n.val*R.val*T.val))/FT;
end


# Print info 
function show(io::IO, g::DislocationCreep)  
    print(io, "DislocationCreep: n=$(g.n.val), r=$(g.r.val), A=$(g.A.val), E=$(g.E.val), V=$(g.V.val), Apparatus=$(g.Apparatus)" )  
end
#-------------------------------------------------------------------------

# Add pre-defined creep laws 
"""
    SetDislocationCreep["Name of Dislocation Creep"]

This is a dictionary with pre-defined creep laws    
"""
SetDislocationCreep = Dict([

# Olivine rheology 
("Dry Olivine | Hirth & Kohlstedt (2003)", 
# after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
# Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters
    DislocationCreep(
        n = 3.5NoUnits,
        A = 1.1e5MPa^(-3.5)/s, 
        E = 530kJ/mol,
        V = 15e-6m^3/mol,
        Apparatus =   "AxialCompression",
        r = 0NoUnits,
        Comment = "Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
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
        """)
    )
)

# Olivine rheology 
("Wet Olivine | Hirth & Kohlstedt (2003)", 
    # After Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    #   Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
    #  Note that this assumes C_OH=1000
    DislocationCreep(
        n = 3.5NoUnits,
        A = 90MPa^(-3.5)/s, 
        E = 480kJ/mol,
        V = 11e-6m^3/mol,
        r   = 1.2NoUnits,
        Apparatus =   "AxialCompression",
        Comment = "Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
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
        """);
    )
)



]); # end of setting pre-defined creep laws
