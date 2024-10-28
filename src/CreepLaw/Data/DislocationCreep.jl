module Dislocation

# This contains all the dislocation creep data
using GeoParams
export dislocation_database, dislocation_database_info

# function dislocation_law_list()
#     m = @__MODULE__
#     @show m
#     out = string.(names(m; all=true, imported=true))
#     filter!(x -> !startswith(x, "#"), out)
#     return [getfield(m, Symbol(x)) for x in out if !isnothing(tryparse(Int, string(x[end]))) || endswith(x, "a") || endswith(x, "b")]
# end

"""
Dislocation creep data for dry olivine after Hirth, G. & Kohlstedt (2003)
"""
function dry_olivine_Hirth_2003()
    # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    # Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters

    data = DislocationCreep(;
        Name="Dry Olivine | Hirth & Kohlstedt (2003)",
        n=3.5NoUnits,
        r=0.0NoUnits,
        A=1.1e5MPa^(-7//2) / s,
        E=530.0kJ / mol,
        V=14e-6m^3 / mol,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

"""
Dislocation creep data for wet olivine (1), constant water fugacity after Hirth, G. & Kohlstedt (2003)
"""
function wet_olivine1_Hirth_2003()
    # After Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    #   Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
    #  Note that this assumes C_OH=1000

    data = DislocationCreep(;
        Name="Wet Olivine 1 | Hirth & Kohlstedt (2003)",
        n=3.5NoUnits,
        A=90MPa^(-7//2) / s,
        E=480kJ / mol,
        V=11e-6m^3 / mol,
        r=1.2NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

"""
Dislocation creep data for wet olivine (2), 5th creep law in table 1 of Hirth & Kohlstedt (2003)
"""
function wet_olivine2_Hirth_2003()
    # After Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    #   Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
    #  Note that this assumes C_OH=1000

    data = DislocationCreep(;
        Name="Wet Olivine 2 | Hirth & Kohlstedt (2003)",
        n=3.5NoUnits,
        A=1600MPa^(-3) / s,
        E=520.0kJ / mol,
        V=22.0e-6m^3 / mol,
        r=1.2NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

"""
Dislocation creep data for Quartz Diorite after Hansen (1982), 'Semibrittle creep of selected crustal rocks at 1000 MPa.' and, Hansen & Carter (1982),  'Flow properties of continental lithosphere.'
"""
function quartz_diorite_HansenCarter_1982()
    #  After Hansen (1982), 'Semibrittle creep of selected crustal rocks at 1000 MPa.' and, Hansen & Carter (1982),
    #  'Flow properties of continental lithosphere.'
    #  Hansen (1982), Fig. 53, page 184 in PDF viewer and table 18, page 224
    #  Carter & Tsenn (1986), table 4, page 18 in PDF viewer

    data = DislocationCreep(;
        Name="Quartz Diorite | Hansen & Carter (1982)",
        n=2.25NoUnits,
        A=3.5e-2MPa^(-9//4) / s,
        E=212kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Verified with the original publication Hansen(1982)(NM). Values checked, plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

"""
Dislocation creep data for Diabase after Caristan (1982), 'The transition from high temperature creep to fracture in Maryland diabase.' and, Bremond (1999), 'Hydrothermalism and diapirism in the Archean: gravitational instability constraints' Bremond (1999), page 5 in text
"""
function diabase_Caristan_1982()
    #  After Caristan (1982), 'The transition from high temperature creep to fracture in Maryland diabase.'
    #  and, Bremond (1999),
    #  'Hydrothermalism and diapirism in the Archean: gravitational instability constraints'
    #  Bremond (1999), page 5 in text

    data = DislocationCreep(;
        Name="Diabase | Caristan (1982)",
        n=3.05NoUnits,
        A=6.12e-2MPa^(-61//20) / s,
        E=276kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Bremond (1999)), plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function serpentinite_Raleigh_1965()
    #  After Raleigh and Paterson (1965), 'Experimental deformation of serpentinite and its tectonic implications'
    #  and, Bremond (1999),
    #  'Hydrothermalism and diapirism in the Archean: gravitational instability constraints'
    #  Bremond (1999), page 5 in text

    data = DislocationCreep(;
        Name="Tumut Pond Serpentinite | Raleigh and Paterson (1965)",
        n=2.8NoUnits,
        A=6.3e-7MPa^(-14//5) / s,
        E=66kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Bremond (1999)), plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function strong_diabase_Mackwell_1998()
    #  After Mackwell et al. (1998), 'High-temperature deformation of dry diabase with application to tectonics on Venus'
    #  Mackwell et al. (1998), page 980, equation in text
    data = DislocationCreep(;
        Name="Maryland strong diabase | Mackwell et al. (1998)",
        n=4.7NoUnits,
        A=8MPa^(-47//10) / s,
        E=485kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Mackwell et al. (1998))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function mafic_granulite_Wilks_1990()
    #  After Li, Gerya and Burg (2010), table 2
    #  referring to Ranalli (1995), 'Rheology of the Earth' (Book), page 334, table 10.3
    #  referring to Wilks and Carter (1990), 'Rheology of some continental lower crustal rocks', Fig. 6, Pikwitonei Granulite

    data = DislocationCreep(;
        Name="Mafic Granulite | Wilks and Carter (1990)",
        n=4.2NoUnits,
        A=1.4e4MPa^(-21//5) / s,
        E=445kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Wilks and Carter (1990))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function wet_quartzite_Ueda_2008()
    #  Ueda et al. (2008), table 1
    data = DislocationCreep(;
        Name="Wet Quartzite | Ueda et al. (2008)",
        n=2.3NoUnits,
        A=1 * exp10(-3.5)MPa^(-23//10) / s,
        E=154kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Ueda et al. (2008))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function granite_Carter_1987()
    #  Huismans et al. (2001), table 2
    #  referring to Carter and Tsenn (1987), 'Flow properties of continental lithosphere', table 4, Westerly Granite (dry)
    #  referring to Hansen and Carter (1983), 'Semibrittle Creep Of Dry And Wet Westerly Granite At 1000 MPa', not accessible
    data = DislocationCreep(;
        Name="Granite | Carter and Tsenn (1987)",
        n=3.3NoUnits,
        A=1.0 * exp10(-5.7)MPa^(-33//10) / s,
        E=186.5kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Carter and Tsenn (1987))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function plagioclase_An75_Ji_1993()
    #  Ranalli (1995), page 334, table 10.3
    #  referring to Ji and Zhao (1993), 'Flow laws of multiphase rocks calculated from experimental data on the constituent phases', table 2 , plagioclase (Ab25An75)
    #  referring to Shelton and Tullis (1981), 'Experimental flow laws for crustal rocks', not accessible
    data = DislocationCreep(;
        Name="Plagioclase An75 | Ji and Zhao (1993)",
        n=3.2NoUnits,
        A=3.27e-4MPa^(-16//5) / s,
        E=238kJ / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Ji and Zhao (1993))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function dry_anorthite_Rybacki_2006()
    # Rybacki, Gottschalk, Wirth and Dresen (2006), table 5
    data = DislocationCreep(;
        Name="Dry Anorthite | Rybacki et al. (2006)",
        n=3.0NoUnits,
        A=exp10(12.7)MPa^(-3) / s,
        E=641kJ / mol,
        V=24e-6m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values taken from Bürgmann & Rybacki (2008) Supplementary table 1. (BK), plots are not reproduced (NM).",
        BibTex_Reference="
            @article{rybacki2006influence,
            title={Influence of water fugacity and activation volume on the flow properties of fine-grained anorthite aggregates},
            author={Rybacki, Erik and Gottschalk, Matthias and Wirth, Richard and Dresen, Georg},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={111},
            number={B3},
            year={2006},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end

function wet_anorthite_Rybacki_2006()
    # Rybacki, Gottschalk, Wirth and Dresen (2006), table 5
    data = DislocationCreep(;
        Name="Wet Anorthite | Rybacki et al. (2006)",
        n=3.0NoUnits,
        A=(10^0.2)MPa^(-3) / s,
        E=345kJ / mol,
        V=38e-6m^3 / mol,
        r=1NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Rybacki, Gottschalk, Wirth and Dresen (2006))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
            @article{rybacki2006influence,
            title={Influence of water fugacity and activation volume on the flow properties of fine-grained anorthite aggregates},
            author={Rybacki, Erik and Gottschalk, Matthias and Wirth, Richard and Dresen, Georg},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={111},
            number={B3},
            year={2006},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end

function wet_olivine_Hirth_2003()
    #  Hirth and Kohlstedt (2003), table 1, no constant C_OH
    data = DislocationCreep(;
        Name="Wet Olivine | Hirth and Kohlstedt (2003)",
        n=3.5NoUnits,
        A=1600.0MPa^(-7//2) / s,
        E=520kJ / mol,
        V=22e-6m^3 / mol,
        r=1.2NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Hirth and Kohlstedt (2003))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
            @article{hirth2003rheology,
            title={Rheology of the upper mantle and the mantle wedge: A view from the experimentalists},
            author={Hirth, Greg and Kohlstedf, D},
            journal={Geophysical monograph-american geophysical union},
            volume={138},
            pages={83--106},
            year={2003},
            publisher={AGU AMERICAN GEOPHYSICAL UNION}
            }
    ",
    )
    return data, info
end

function wet_quartzite_Kirby_1983()
    #  Li, Gerya and Burg (2010), table 2
    #  Ranalli (1995), 'Rheology of the Earth' (Book)
    #  referring to Kirby (1983), table 2, first quartzite (wet)
    #  referring to Koch et al. (1981), unpublished manuscript...
    data = DislocationCreep(;
        Name="Wet Quarzite | Kirby (1983)",
        n=2.4NoUnits,
        A=1.25e-2MPa^(-2.4), #A is given in 10^(3.2)^ GPa^(2.4)^ in paper
        E=160e3J / mol,
        V=0m^3 / mol,
        r=0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Kirby (1983))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function wet_olivine_Afonso_2004()
    #  Schmalholz, Kaus, Burg (2009), table 1
    #  referring to Afonso and Ranalli (2004), table 1, wet peridotite
    #  referring to Chopra and Paterson papers, but values dont fit the Afonso and Ranalli (2004) ones
    data = DislocationCreep(;
        Name="Wet Upper Mantle Olivine | Afonso and Ranalli (2004)",
        n=4.0NoUnits,
        A=2.0e3MPa^(-4) / s,
        E=471kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Afonso and Ranalli (2004))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function granite_Tirel_2008()
    #  Tirel et al. (2008), table 1
    #  referring to Kirby and Kronenberg (1987), table 3
    #  different values for n and A in Kirby and Kronenberg (1987) compared with Tirel et al. (2008)
    data = DislocationCreep(;
        Name="Granite | Tirel et al. (2008)",
        n=3.2NoUnits,
        A=1.25e-9MPa^(-3.2) / s,
        E=123kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (Tirel et al. (2008))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
            @article{tirel2008dynamics,
            title={Dynamics and structural development of metamorphic core complexes},
            author={Tirel, C{\'e}line and Brun, Jean-Pierre and Burov, Evgueni},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={113},
            number={B4},
            year={2008},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end

function dry_olivine_Gerya_2019()
    # after Exercise 6.1 of Numerical Geodynamics
    data = DislocationCreep(;
        Name="Dry Olivine | Gerya (2019)",
        n=3.5NoUnits,
        #A = 2.5e-17Pa^(-7//2)/s,
        A=uconvert(MPa^(-7//2) / s, 2.5e-17Pa^(-7//2) / s),
        E=532.0kJ / mol,
        V=0.0m^3 / mol,
        Apparatus=AxialCompression,   # used in book (according to matlab script)
        r=0.0NoUnits,
    )
    info = MaterialParamsInfo(;
        Comment="This is from Exercise 6.1, indicated to be valid for olivine under upper mantle conditions.",
        BibTex_Reference="
            @book{Gerya_2019,
            title={Introduction to Numerical Geodynamic Modelling},
            ISBN={978-1-107-14314-2},
            note={Google-Books-ID: 8XGSDwAAQBAJ},
            publisher={Cambridge University Press},
            author={Gerya, Taras},
            year={2019}, month={May},
            language={en} }
    ",
    )
    return data, info
end

function rock_salt_Li_Urai_2016()
    data = DislocationCreep(;
        Name="Rock salt | Li & Urai (2016)",
        n=5.0NoUnits,
        A=7.26e-6MPa^(-5) / s,
        E=53.92kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked in (Li & Urai, (2016)) are different from given source (Wawersik & Zeuch, (1986))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
        @article{li2016rheology,
        title={Rheology of rock salt for salt tectonics modeling},
        author={Li, Shi-Yuan and Urai, Janos L},
        journal={Petroleum science},
        volume={13},
        number={4},
        pages={712--724},
        year={2016},
        publisher={Springer}
        }
    ",
    )
    return data, info
end

function salado_rock_salt_Li_2016()
    data = DislocationCreep(;
        Name="Salado Rock salt | Li & Urai (2016)",
        n=5.0NoUnits,
        A=7.26e-6MPa^(-5) / s,
        E=53.92kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked in (Li & Urai, (2016)) are different from given source (Wawersik & Zeuch, (1986))(NM), plots are not reproduced (NM).",
        BibTex_Reference="
        @article{li2016rheology,
        title={Rheology of rock salt for salt tectonics modeling},
        author={Li, Shi-Yuan and Urai, Janos L},
        journal={Petroleum science},
        volume={13},
        number={4},
        pages={712--724},
        year={2016},
        publisher={Springer}
        }
    ",
    )
    return data, info
end

function wet_olivine_Mei_2000b()
    #  Mei & Kohlstedt (2000b), table 1
    data = DislocationCreep(;
        Name="Wet Olivine | Mei & Kohlstedt (2000b)",
        n=3.02NoUnits,
        A=(1.5e3)MPa^(-4) / s,
        E=470.0kJ / mol,
        V=20.0e-6m^3 / mol,
        r=0.98NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="The A value is not exactly the same as in Mei & Kohlstedt (2000b) but is approximated (NM).",
        BibTex_Reference="
        @article{mei2000influence,
        title={Influence of water on plastic deformation of olivine aggregates: 2. Dislocation creep regime},
        author={Mei, S and Kohlstedt, DL},
        journal={Journal of Geophysical Research: Solid Earth},
        volume={105},
        number={B9},
        pages={21471--21481},
        year={2000},
        publisher={Wiley Online Library}
        }
    ",
    )
    return data, info
end

function dry_olivine_Karato_2003()
    #  Mei & Kohlstedt (2000b), abstract
    data = DislocationCreep(;
        Name="Dry Olivine | Karato & Jung (2003)",
        n=3.0NoUnits,
        A=(10^6.1)MPa^(-3) / s,
        E=510.0kJ / mol,
        V=14.0e-6m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @article{karato2003effects,
        title={Effects of pressure on high-temperature dislocation creep in olivine},
        author={Karato, Shun-Ichiro and Jung, Haemyeong},
        journal={Philosophical Magazine},
        volume={83},
        number={3},
        pages={401--414},
        year={2003},
        publisher={Taylor \\& Francis}
        }
    ",
    )
    return data, info
end

function wet_olivine_Karato_2003()
    #  Karato & Jung (2003), abstract
    data = DislocationCreep(;
        Name="Wet Olivine | Karato & Jung (2003)",
        n=3.0NoUnits,
        A=(10^2.9)MPa^(-4.2) / s,
        E=470.0kJ / mol,
        V=24.0e-6m^3 / mol,
        r=1.2NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM, PA)",
        BibTex_Reference="
        @article{karato2003effects,
        title={Effects of pressure on high-temperature dislocation creep in olivine},
        author={Karato, Shun-Ichiro and Jung, Haemyeong},
        journal={Philosophical Magazine},
        volume={83},
        number={3},
        pages={401--414},
        year={2003},
        publisher={Taylor \\& Francis}
        }
    ",
    )
    return data, info
end

function wet_clinopyroxene_Chen_2006()
    #  Chen et al. (2006), section 4. Discussion, equation (3)
    data = DislocationCreep(;
        Name="Wet Clinopyroxene | Chen et al. (2006)",
        n=2.7NoUnits,
        A=(10^6.7)MPa^(-27//10) / s,
        E=670.0kJ / mol,
        V=0.0m^3 / mol,
        r=3.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @article{chen2006water,
        title={Water weakening of clinopyroxene in the dislocation creep regime},
        author={Chen, S and Hiraga, T and Kohlstedt, David L},
        journal={Journal of Geophysical Research: Solid Earth},
        volume={111},
        number={B8},
        year={2006},
        publisher={Wiley Online Library}
        }
    ",
    )
    return data, info
end

function dry_clinopyroxene_Bystricky_Mackwell_2001()
    #  Bystricky & Mackwell (2001), section 4. Discussion, equation (3)
    data = DislocationCreep(;
        Name="Dry Clinopyroxene | Bystricky & Mackwell (2001)",
        n=4.7NoUnits,
        A=(10^9.8)MPa^(-47//10) / s,
        E=760.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @article{bystricky2001creep,
        title={Creep of dry clinopyroxene aggregates},
        author={Bystricky, Misha and Mackwell, Stephen},
        journal={Journal of Geophysical Research: Solid Earth},
        volume={106},
        number={B7},
        pages={13443--13454},
        year={2001},
        publisher={Wiley Online Library}
        }
    ",
    )
    return data, info
end

function dry_diopside_Dimanov_2005()
    #  Dimanov & Dresen (2005), table 3b
    data = DislocationCreep(;
        Name="Dry Diopside | Dimanov & Dresen (2005)",
        n=5.5NoUnits,
        A=uconvert(MPa^(-55//10) / s, 3.01e-28Pa^(-55//10) / s),
        E=691.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @article{dimanov2005rheology,
        title={Rheology of synthetic anorthite-diopside aggregates: Implications for ductile shear zones},
        author={Dimanov, Alexandre and Dresen, Georg},
        journal={Journal of Geophysical Research: Solid Earth},
        volume={110},
        number={B7},
        year={2005},
        publisher={Wiley Online Library}
        }
    ",
    )
    return data, info
end

function wet_diopside_Dimanov_2005()
    #  Dimanov & Dresen (2005), table 3b
    data = DislocationCreep(;
        Name="Wet Diopside | Dimanov & Dresen (2005)",
        n=5.5NoUnits,
        A=uconvert(MPa^(-55//10) / s, 5.16e-33Pa^(-55//10) / s),
        E=534.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @article{dimanov2005rheology,
        title={Rheology of synthetic anorthite-diopside aggregates: Implications for ductile shear zones},
        author={Dimanov, Alexandre and Dresen, Georg},
        journal={Journal of Geophysical Research: Solid Earth},
        volume={110},
        number={B7},
        year={2005},
        publisher={Wiley Online Library}
        }
    ",
    )
    return data, info
end

function wet_omphacite_Zhang_2006()
    #  Zhang et al. (2006), equation (4)
    data = DislocationCreep(;
        Name="Wet Omphacite | Zhang et al. (2006)",
        n=3.5NoUnits,
        A=(10^-2)MPa^(-7//2) / s,
        E=310.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @article{zhang2006rheology,
        title={Rheology of omphacite at high temperature and pressure and significance of its lattice preferred orientations},
        author={Zhang, Junfeng and Green II, Harry W and Bozhilov, Krassimir N},
        journal={Earth and Planetary Science Letters},
        volume={246},
        number={3-4},
        pages={432--443},
        year={2006},
        publisher={Elsevier}
        }
    ",
    )
    return data, info
end

function wet_jadeit_Orzol_2006()
    #  Orzol et al. (2006), page 11
    data = DislocationCreep(;
        Name="Wet Jadeit | Orzol et al. (2006)",
        n=3.7NoUnits,
        A=(10^-3.3)MPa^(-37//10) / s,
        E=326.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @article{orzol2006experimental,
        title={Experimental deformation of synthetic wet jadeite aggregates},
        author={Orzol, J and St{\"o}ckhert, B and Trepmann, CA and Rummel, F},
        journal={Journal of Geophysical Research: Solid Earth},
        volume={111},
        number={B6},
        year={2006},
        publisher={Wiley Online Library}
        }
    ",
    )
    return data, info
end

function dry_anorthite_Rybacki_2000()
    #  Rybacki & Dresen (2000), table 5
    data = DislocationCreep(;
        Name="Dry Anorthite | Rybacki & Dresen (2000)",
        n=3.0NoUnits,
        A=(10^12.7)MPa^(-3) / s,
        E=648.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function wet_anorthite_Rybacki_2000()
    #  Rybacki & Dresen (2000), table 5
    data = DislocationCreep(;
        Name="Wet Anorthite | Rybacki & Dresen (2000)",
        n=3.0NoUnits,
        A=(10^2.6)MPa^(-3) / s,
        E=356.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
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
    ",
    )
    return data, info
end

function wet_quartzite_Rutter_2004()
    #  Rutter & Brodie (2004), table 5
    data = DislocationCreep(;
        Name="Wet Quartzite | Rutter & Brodie (2004)",
        n=2.97NoUnits,
        A=(10^-4.93)MPa^(-2.97) / s,
        E=242.0kJ / mol,
        V=0.0m^3 / mol,
        r=1.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @article{rutter2004experimental,
        title={Experimental intracrystalline plastic flow in hot-pressed synthetic quartzite prepared from Brazilian quartz crystals},
        author={Rutter, EH and Brodie, KH},
        journal={Journal of Structural Geology},
        volume={26},
        number={2},
        pages={259--270},
        year={2004},
        publisher={Elsevier}
        }
    ",
    )
    return data, info
end

function wet_quartzite_Hirth_2001()
    #  Hirth et al. (2001), table 5
    data = DislocationCreep(;
        Name="Wet Quartzite | Hirth et al. (2001)",
        n=4.0NoUnits,
        A=(10^-11.2)MPa^(-4) / s,
        E=135.0kJ / mol,
        V=0.0m^3 / mol,
        r=1.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
        @article{hirth2001evaluation,
        title={An evaluation of quartzite flow laws based on comparisons between experimentally and naturally deformed rocks},
        author={Hirth, Greg and Teyssier, Christian and Dunlap, James W},
        journal={International Journal of Earth Sciences},
        volume={90},
        number={1},
        pages={77--87},
        year={2001},
        publisher={Springer}
        }
    ",
    )
    return data, info
end

function dry_quartzite_Jaoul_1984()
    #  Jaoul et al. (1984), table 1, first entry
    data = DislocationCreep(;
        Name="Dry Quartzite | Jaoul et al. (1984)",
        n=2.8NoUnits,
        A=(2.8899e-3)MPa^(-2.8) / s,
        E=184.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
            @article{jaoul1984effect,
            title={The effect of varying water contents on the creep behavior of Heavitree quartzite},
            author={Jaoul, Olivier and Tullis, Jan and Kronenberg, Andreas},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={89},
            number={B6},
            pages={4298--4312},
            year={1984},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end

function wet_quartzite_Jaoul_1984()
    #  Jaoul et al. (1984), table 1, second entry
    data = DislocationCreep(;
        Name="Wet Quartzite | Jaoul et al. (1984)",
        n=1.4NoUnits,
        A=(104.957e-3)MPa^(-14//5) / s,
        E=146.44kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
            @article{jaoul1984effect,
            title={The effect of varying water contents on the creep behavior of Heavitree quartzite},
            author={Jaoul, Olivier and Tullis, Jan and Kronenberg, Andreas},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={89},
            number={B6},
            pages={4298--4312},
            year={1984},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end

function wet_quartzite_Tokle_2019()
    #  Tokle et al. (2019), table 1, 2nd extrapolated fit
    data = DislocationCreep(;
        Name="Wet Quartzite | Tokle et al. (2019)",
        n=3.0NoUnits,
        A=(10^-11.959)MPa^(-3) / s,
        E=115.0kJ / mol,
        V=0.0m^3 / mol,
        r=1.2NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
            @article{tokle2019flow,
            title={Flow laws and fabric transitions in wet quartzite},
            author={Tokle, Leif and Hirth, Greg and Behr, Whitney M},
            journal={Earth and Planetary Science Letters},
            volume={505},
            pages={152--161},
            year={2019},
            publisher={Elsevier}
            }
    ",
    )
    return data, info
end

function wet_quartzite_Lu_2019()
    #  Lu and Jiang (2019), section 3, equation (6)
    data = DislocationCreep(;
        Name="Wet Quartzite | Lu and Jiang (2019)",
        n=4.0NoUnits,
        A=(6.0e-15)MPa^(-6.7) / s,
        E=132.0kJ / mol,
        V=35.3e-6m^3 / mol,
        r=2.7NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
            @article{lu2019quartz,
            title={Quartz flow law revisited: the significance of pressure dependence of the activation enthalpy},
            author={Lu, Lucy X and Jiang, Dazhi},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={124},
            number={1},
            pages={241--256},
            year={2019},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end

function lowP_wet_quartzite_Lusk_2021()
    #  Lusk et al. (2021), abstract, 1st law
    data = DislocationCreep(;
        Name="low pressure wet Quartzite | Lusk et al. (2021)",
        n=-3.99NoUnits,
        A=(10^-9.3)MPa^(-7//2) / s,
        E=118.0kJ / mol,
        V=2.59e-6m^3 / mol,
        r=0.49NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
            @article{lusk2021natural,
            title={Natural and Experimental Constraints on a Flow Law for Dislocation-Dominated Creep in Wet Quartz},
            author={Lusk, Alexander DJ and Platt, John P and Platt, Jason A},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={126},
            number={5},
            pages={e2020JB021302},
            year={2021},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end

function wet_quartzite_Lusk_2021()
    #  Lusk et al. (2021), table 2, full data set
    data = DislocationCreep(;
        Name="Wet Quartzite | Lusk et al. (2021)",
        n=2.1NoUnits,
        A=(10^-7.9)MPa^(-1.51) / s,
        E=94.0kJ / mol,
        V=1.44e-6m^3 / mol,
        r=0.2NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
            @article{lusk2021natural,
            title={Natural and Experimental Constraints on a Flow Law for Dislocation-Dominated Creep in Wet Quartz},
            author={Lusk, Alexander DJ and Platt, John P and Platt, Jason A},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={126},
            number={5},
            pages={e2020JB021302},
            year={2021},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end


function highP_wet_quartzite_Lusk_2021()
    #  Lusk et al. (2021), abstract, 2nd law
    data = DislocationCreep(;
        Name="high pressure wet Quartzite | Lusk et al. (2021)",
        n=2.0NoUnits,
        A=(10^-6.36)MPa^(-2.3) / s,
        E=77.0kJ / mol,
        V=2.59e-6m^3 / mol,
        r=0.49NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (NM).",
        BibTex_Reference="
            @article{lusk2021natural,
            title={Natural and Experimental Constraints on a Flow Law for Dislocation-Dominated Creep in Wet Quartz},
            author={Lusk, Alexander DJ and Platt, John P and Platt, Jason A},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={126},
            number={5},
            pages={e2020JB021302},
            year={2021},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end

function high_stress_wet_dunite_Chopra_1981()
    # Chopra & Paterson (1981) Table IV, Anita Bay dunite ($$\sigma$$ > ~ 100MPa)
    data = DislocationCreep(;
        Name="High Stress Wet Dunite | Chopra & Paterson (1981)",
        n=3.35NoUnits,
        A=(10^3.98)MPa^(-3.35) / s,
        E=444.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (PA).",
        BibTex_Reference="
            @article{CHOPRA1981453,
            title = {The experimental deformation of dunite},
            journal = {Tectonophysics},
            volume = {78},
            number = {1},
            pages = {453-473},
            year = {1981},
            note = {The Effect of Deformation on Rocks},
            issn = {0040-1951},
            doi = {https://doi.org/10.1016/0040-1951(81)90024-X},
            url = {https://www.sciencedirect.com/science/article/pii/004019518190024X},
            author = {P.N. Chopra and M.S. Paterson},
            abstract = {Deformation experiments have been carried out on two dunites (Anita Bay, of 100 μm grain size, and Åheim, of 900 μm grain size) at strain rates from 10−3 to 10−6 s−1 and temperatures from 1000°C to 1300°C in a gas-medium deformation apparatus at 300 MPa confining pressure. Most of the tests were under “wet” conditions defined by the presence of small amounts of water from hydrous minerals initially present. Constant strain rate and relaxation experiments, covering ranges of flow stress down to about 70 MPa and 7 MPa, respectively, show that there is a change in flow law in going below about 100 MPa differential stress, and that the coarser-grained rock is stronger than the finer-grained one. Power law parameters above the transition are n = 4.48 ± 0.31 and Q = 498 ± 38 kJ mol−1 for Åheim dunite and n = 3.35 ± 0.17 and Q = 444 ± 24 kJ mol−1 for Anita Bay dunite, while below the transition relaxation tests on Anita Bay dunite give n = 2.44 ±0.18 and Q = 386 ± 27 kJ mol−1. It is concluded that there is a weakening effect of water, that this effect is mainly in the grain boundaries and that grain boundary sliding is probably a significant deformation mechanism at the lower stresses under wet conditions.}
            }
    ",
    )
    return data, info
end



function low_stress_wet_dunite_Chopra_1981()
    # Chopra & Paterson (1981) Table IV, Anita Bay dunite ($$\sigma$$ < ~ 100MPa)
    data = DislocationCreep(;
        Name="Low Stress Wet Dunite | Chopra & Paterson (1981)",
        n=2.44NoUnits,
        A=(10^3.88)MPa^(-2.44) / s,
        E=386.0kJ / mol,
        V=0.0m^3 / mol,
        r=0.0NoUnits,
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Values checked (PA).",
        BibTex_Reference="
            @article{CHOPRA1981453,
            title = {The experimental deformation of dunite},
            journal = {Tectonophysics},
            volume = {78},
            number = {1},
            pages = {453-473},
            year = {1981},
            note = {The Effect of Deformation on Rocks},
            issn = {0040-1951},
            doi = {https://doi.org/10.1016/0040-1951(81)90024-X},
            url = {https://www.sciencedirect.com/science/article/pii/004019518190024X},
            author = {P.N. Chopra and M.S. Paterson},
            abstract = {Deformation experiments have been carried out on two dunites (Anita Bay, of 100 μm grain size, and Åheim, of 900 μm grain size) at strain rates from 10−3 to 10−6 s−1 and temperatures from 1000°C to 1300°C in a gas-medium deformation apparatus at 300 MPa confining pressure. Most of the tests were under “wet” conditions defined by the presence of small amounts of water from hydrous minerals initially present. Constant strain rate and relaxation experiments, covering ranges of flow stress down to about 70 MPa and 7 MPa, respectively, show that there is a change in flow law in going below about 100 MPa differential stress, and that the coarser-grained rock is stronger than the finer-grained one. Power law parameters above the transition are n = 4.48 ± 0.31 and Q = 498 ± 38
            }
    ",
    )
    return data, info
end


@inline dislocation_database(f::F) where {F} = first(f())
@inline dislocation_database_info(f::F) where {F} = last(f())

end
