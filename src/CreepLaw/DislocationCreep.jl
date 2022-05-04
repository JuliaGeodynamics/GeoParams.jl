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
    DislocationCreep(n = 1.0NoUnits, r = 0.00.0NoUnits, A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = AxialCompression )
    
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
either AxialCompression, SimpleShear or Invariant.
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

function computeCreepLaw_EpsII(TauII, a::DislocationCreep; P::_R, T::_R, f::_R, args...) where _R<:Real
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
function computeCreepLaw_TauII(EpsII, a::DislocationCreep; P::_R, T::_R, f::_R, args...) where _R<:Real
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

# Quartz Diorite rheology
("Quartz Diorite | Hansen & Carter (1982)", 
    #  After Hansen (1982), 'Semibrittle creep of selected crustal rocks at 1000 MPa.' and, Hansen & Carter (1982),
    #  'Flow properties of continental lithosphere.'
    #  Hansen (1982), Fig. 53, page 184 in PDF viewer and table 18, page 224
    #  Carter & Tsenn (1986), table 4, page 18 in PDF viewer
    (DislocationCreep(
        Name = "Quartz Diorite | Hansen & Carter (1982)",
        n = 2.25NoUnits,
        A = 3.5e-2MPa^(-2.25)/s, 
        E = 212kJ/mol,
        V = 0m^3/mol,
        r   = 0NoUnits,
        Apparatus =   AxialCompression),
        MaterialParamsInfo(Comment = "Verified with the original publication Hansen(1982)(NM). Values checked, plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
             @article{carter1982stress,
             title={Stress magnitudes in natural rock salt},
             author={Carter, Neville L and Hansen, Francis D and Senseny, Paul E},
             journal={Journal of Geophysical Research: Solid Earth},
             volume={87},
             number={B11},
             pages={9289--9300},
             year={1982},
             publisher={Wiley Online Library}
             }
        """))
    )
)
                
# Diabase rheology
("Diabase | Caristan (1982)", 
    #  After Caristan (1982), 'The transition from high temperature creep to fracture in Maryland diabase.'
    #  and, Bremond (1999),
    #  'Hydrothermalism and diapirism in the Archean: gravitational instability constraints'
    #  Bremond (1999), page 5 in text
    (DislocationCreep(
        Name = "Diabase | Caristan (1982)",
        n = 3.05NoUnits,
        A = 6.0e-2MPa^(-3.05)/s, 
        E = 276kJ/mol,
        V = 0m^3/mol,
        r   = 0NoUnits,
        Apparatus =   AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Bremond (1999)), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
             @article{caristan1982transition,
             title={The transition from high temperature creep to fracture in Maryland diabase},
             author={Caristan, Y},
             journal={Journal of Geophysical Research: Solid Earth},
             volume={87},
             number={B8},
             pages={6781--6790},
             year={1982},
             publisher={Wiley Online Library}
             }
        """))
    )
)

# Tumut Pond Serpentinite rheology
("Tumut Pond Serpentinite | Raleigh and Paterson (1965)", 
    #  After Raleigh and Paterson (1965), 'Experimental deformation of serpentinite and its tectonic implications'
    #  and, Bremond (1999),
    #  'Hydrothermalism and diapirism in the Archean: gravitational instability constraints'
    #  Bremond (1999), page 5 in text
    (DislocationCreep(
        Name = "Tumut Pond Serpentinite | Raleigh and Paterson (1965)",
        n = 2.8NoUnits,
        A = 6.3e-7MPa^(-2.8)/s, 
        E = 66kJ/mol,
        V = 0m^3/mol,
        r   = 0NoUnits,
        Apparatus =   AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Bremond (1999)), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
             @article{raleigh1965experimental,
             title={Experimental deformation of serpentinite and its tectonic implications},
             author={Raleigh, CB and Paterson, MS},
             journal={Journal of Geophysical Research},
             volume={70},
             number={16},
             pages={3965--3985},
             year={1965},
             publisher={Wiley Online Library}
             }
        """))
    )
)

# Maryland strong diabase rheology
("Maryland strong diabse | Mackwell et al. (1998)", 
    #  After Mackwell et al. (1998), 'High-temperatur deformation of dry diabase with application to tectonics on Venus'
    #  Mackwell et al. (1998), page 980, equation in text
    (DislocationCreep(
        Name = "Maryland strong diabse | Mackwell et al. (1998)",
        n = 4.7NoUnits,
        A = 8MPa^(-4.7)/s, 
        E = 485kJ/mol,
        V = 0m^3/mol,
        r = 0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Mackwell et al. (1998))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{mackwell1998high,
            title={High-temperature deformation of dry diabase with application to tectonics on Venus},
            author={Mackwell, SJ and Zimmerman, ME and Kohlstedt, DL},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={103},
            number={B1},
            pages={975--984},
            year={1998},
            publisher={Wiley Online Library}
            }
        """))
    )
)       

# Mafic Granulite rheology
("Mafic Granulite | Wilks and Carter (1990)", 
    #  After Li, Gerya and Burg (2010), table 2
    #  referring to Ranalli (1995), 'Rheology of the Earth' (Book), page 334, table 10.3
    #  referring to Wilks and Carter (1990), 'Rheology of some continental lower crustal rocks', Fig. 6, Pikwitonei Granulite
    (DislocationCreep(
        Name = "Mafic Granulite | Wilks and Carter (1990)",
        n = 4.2NoUnits,
        A = 1.4e4MPa^(-4.2)/s, 
        E = 445kJ/mol,
        V = 0m^3/mol,
        r = 0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Wilks and Carter (1990))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{wilks1990rheology,
            title={Rheology of some continental lower crustal rocks},
            author={Wilks, Kenneth R and Carter, Neville L},
            journal={Tectonophysics},
            volume={182},
            number={1-2},
            pages={57--77},
            year={1990},
            publisher={Elsevier}
            }
        """))
    )
)

# Wet Quartzite rheology
("Wet Quartzite | Ueda et al. (2008)", 
    #  Ueda et al. (2008), table 1
    (DislocationCreep(
        Name = "Wet Quartzite | Ueda et al. (2008)",
        n = 2.3NoUnits,
        A = 1*exp10(-3.5)MPa^(-2.3)/s, 
        E = 154kJ/mol,
        V = 0m^3/mol,
        r = 0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Ueda et al. (2008))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{ueda2008subduction,
            title={Subduction initiation by thermal--chemical plumes: numerical studies},
            author={Ueda, Kosuke and Gerya, Taras and Sobolev, Stephan V},
            journal={Physics of the Earth and Planetary Interiors},
            volume={171},
            number={1-4},
            pages={296--312},
            year={2008},
            publisher={Elsevier}
            }
        """))
    )
)

# Granite rheology
("Granite | Carter and Tsenn (1987)", 
    #  Huismans et al. (2001), table 2
    #  referring to Carter and Tsenn (1987), 'Flow properties of continental lithosphere', table 4, Westerly Granite (dry)
    #  referring to Hansen and Carter (1983), 'Semibrittle Creep Of Dry And Wet Westerly Granite At 1000 MPa', not accessable
    (DislocationCreep(
        Name = "Granite | Carter and Tsenn (1987)",
        n = 3.3NoUnits,
        A = 1.0*exp10(-5.7)MPa^(-3.3)/s,
        E = 186.5kJ/mol,
        V = 0m^3/mol,
        r = 0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Carter and Tsenn (1987))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{carter1987flow,
            title={Flow properties of continental lithosphere},
            author={Carter, Neville L and Tsenn, Michael C},
            journal={Tectonophysics},
            volume={136},
            number={1-2},
            pages={27--63},
            year={1987},
            publisher={Elsevier}
            }
        """))
    )
)

# Plagioclase An75 rheology
("Plagioclase An75 | Ji and Zhao (1993)", 
    #  Ranalli (1995), page 334, table 10.3
    #  referring to Ji and Zhao (1993), 'Flow laws of multiphase rocks calculated from experimental data on the constituent phases', table 2 , plagioclase (Ab25An75)
    #  referring to Shelton and Tullis (1981), 'Experimental flow laws for crustal rocks', not accessable
   (DislocationCreep(
        Name = "Plagioclase An75 | Ji and Zhao (1993)",
        n = 3.2NoUnits,
        A = 3.27e-4MPa^(-3.2)/s, 
        E = 238kJ/mol,
        V = 0m^3/mol,
        r = 0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Ji and Zhao (1993))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{ji1993flow,
            title={Flow laws of multiphase rocks calculated from experimental data on the constituent phases},
            author={Ji, Shaocheng and Zhao, Pinglao},
            journal={Earth and Planetary Science Letters},
            volume={117},
            number={1-2},
            pages={181--187},
            year={1993},
            publisher={Elsevier}
            }
        """))
    )
)

# Dry Anorthite rheology
("Dry Anorthite | Rybecki and Dresen (2000)", 
    #  Rybecki and Dresen (2000), table 2 + table 3
    (DislocationCreep(
        Name = "Dry Anorthite | Rybecki and Dresen (2000)",
        n = 3.0NoUnits,
        A = 1.0*exp10(-12.7)MPa^(-3.0)/s, 
        E = 648kJ/mol,
        V = 0m^3/mol,
        r = 0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Rybecki and Dresen (2000))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{rybacki2000dislocation,
            title={Dislocation and diffusion creep of synthetic anorthite aggregates},
            author={Rybacki, Erik and Dresen, Georg},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={105},
            number={B11},
            pages={26017--26036},
            year={2000},
            publisher={Wiley Online Library}
            }
        """))
    )
)

# Wet Anorthite rheology
("Wet Anorthite | Rybecki and Dresen (2000)", 
    #  Rybecki and Dresen (2000), table 2 + table 3
    (DislocationCreep(
        Name = "Wet Anorthite | Rybecki and Dresen (2000)",
        n = 3.0NoUnits,
        A = 1.0*exp10(-2.6)MPa^(-3.0)/s, 
        E = 356kJ/mol,
        V = 0m^3/mol,
        r = 0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Rybecki and Dresen (2000))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{rybacki2000dislocation,
            title={Dislocation and diffusion creep of synthetic anorthite aggregates},
            author={Rybacki, Erik and Dresen, Georg},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={105},
            number={B11},
            pages={26017--26036},
            year={2000},
            publisher={Wiley Online Library}
            }
        """))
    )
)

# Wet Olivine rheology
("Wet Olivine | Hirth and Kohlstedt (2003)", 
    #  Hirth and Kohlstedt (2003), table 1, no constant C_OH
    (DislocationCreep(
        Name = "Wet Olivine | Hirth and Kohlstedt (2003)",
        n = 3.5NoUnits,
        A = 1600.0MPa^(-3.5)/s, 
        E = 520kJ/mol,
        V = 11.0e-6m^3/mol,
        r = 1.2NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Hirth and Kohlstedt (2003))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{hirth2003rheology,
            title={Rheology of the upper mantle and the mantle wedge: A view from the experimentalists},
            author={Hirth, Greg and Kohlstedf, D},
            journal={Geophysical monograph-american geophysical union},
            volume={138},
            pages={83--106},
            year={2003},
            publisher={AGU AMERICAN GEOPHYSICAL UNION}
            }
        """))
    )
)
                
# Wet Quarzite rheology
("Wet Quarzite | Kirby (1983)", 
    #  Li, Gerya and Burg (2010), table 2
    #  Ranalli (1995), 'Rheology of the Earth' (Book)
    #  referring to Kirby (1983), table 2, first quartzite (wet)
    #  referring to Koch et al. (1981), unpublished manuscript...
    (DislocationCreep(
        Name = "Wet Quarzite | Kirby (1983)",
        n = 2.3NoUnits,
        A = 3.2e-4MPa^(-2.3)/s, 
        E = 154kJ/mol,
        V = 0m^3/mol,
        r = 0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Kirby (1983))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{kirby1983rheology,
            title={Rheology of the lithosphere},
            author={Kirby, Stephen H},
            journal={Reviews of Geophysics},
            volume={21},
            number={6},
            pages={1458--1487},
            year={1983},
            publisher={Wiley Online Library}
            }
        """))
    )
)

# Wet Upper Mantle Olivine rheology
("Wet Upper Mantle Olivine | Afonso and Ranalli (2004)", 
    #  Schmalholz, Kaus, Burg (2009), table 1
    #  referring to Afonso and Ranalli (2004), table 1, wet peridotite
    #  referring to Chopra and Paterson papers, but values dont fit the Afonso and Ranalli (2004) ones
    (DislocationCreep(
        Name = "Wet Upper Mantle Olivine | Afonso and Ranalli (2004)",
        n = 4.0NoUnits,
        A = 2.0e3MPa^(-4.0)/s, 
        E = 471kJ/mol,
        V = 0.0m^3/mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Afonso and Ranalli (2004))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{afonso2004crustal,
            title={Crustal and mantle strengths in continental lithosphere: is the jelly sandwich model obsolete?},
            author={Afonso, Juan Carlos and Ranalli, Giorgio},
            journal={Tectonophysics},
            volume={394},
            number={3-4},
            pages={221--232},
            year={2004},
            publisher={Elsevier}
            }
        """))
    )
)

# Granite rheology
("Granite | Tirel et al. (2008)", 
    #  Tirel et al. (2008), table 1
    #  referring to Kirby and Kronenberg (1987), table 3
    #  different values for n and A in Kirby and Kronenberg (1987) compared with Tirel et al. (2008)
    (DislocationCreep(
        Name = "Granite | Tirel et al. (2008)",
        n = 3.2NoUnits,
        A = 1.25e-9MPa^(-3.2)/s, 
        E = 123kJ/mol,
        V = 0.0m^3/mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression),
        MaterialParamsInfo(Comment = "Values checked (Tirel et al. (2008))(NM), plots are not reproduced (NM).",
        BibTex_Reference = parse_bibtex("""
            @article{tirel2008dynamics,
            title={Dynamics and structural development of metamorphic core complexes},
            author={Tirel, C{\'e}line and Brun, Jean-Pierre and Burov, Evgueni},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={113},
            number={B4},
            year={2008},
            publisher={Wiley Online Library}
            }
        """))
    )
)

]); # end of setting pre-defined creep laws
