# This contains predefined dislocation creep values - Feel free to expand
export DislocationCreep_info
"""
    SetDislocationCreep["Name of Dislocation Creep"]

Sets predefined dislocation creep data from a dictionary

"""
const DislocationCreep_info = Dict(
    [

        # Olivine rheology
        (
            "Dry Olivine | Hirth & Kohlstedt (2003)",
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            # Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters
            (
                DislocationCreep(;
                    Name = "Dry Olivine | Hirth & Kohlstedt (2003)",
                    n = 3.5NoUnits,
                    r = 0.0NoUnits,
                    A = 1.1e5MPa^(-7 // 2) / s,
                    E = 530.0kJ / mol,
                    V = 14.0e-6m^3 / mol,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Test rheology
        (
            "Test Disl",
            # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            # Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters
            (
                DislocationCreep(;
                    Name = "Test Disl",
                    n = 3.5NoUnits,
                    r = 0.0NoUnits,
                    A = 7.440446357925734e7MPa^(-7 // 2) / s,
                    E = 530.0kJ / mol,
                    V = 14.0e-6m^3 / mol,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Law to check Thorsten Beckers book rheology and its conversion of the A factor",
                    BibTex_Reference = "",
                ),
            ),
        )

        # Olivine rheology, constant water fugacity
        (
            "1. Wet Olivine | Hirth & Kohlstedt (2003)",
            # After Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
            #   Inside the subduction Factory 83?105. Table 1, "wet dislocation" parameters
            #  Note that this assumes C_OH=1000
            (
                DislocationCreep(;
                    Name = "1. Wet Olivine | Hirth & Kohlstedt (2003)",
                    n = 3.5NoUnits,
                    A = 90MPa^(-7 // 2) / s,
                    E = 480kJ / mol,
                    V = 11.0e-6m^3 / mol,
                    r = 1.2NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Still to be verified with the original publication (BK). Values checked, plots are not reproduced (DK).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Wet olivine rheology, 5th creep law in table 1 of Hirth & Kohlstedt (2003)
        (
            "2. Wet Olivine | Hirth & Kohlstedt (2003)",
            #  Hirth & Kohlstedt (2003), table 1
            (
                DislocationCreep(;
                    Name = "2. Wet Olivine | Hirth & Kohlstedt (2003)",
                    n = 3.0NoUnits,
                    A = 1600MPa^(-3) / s,
                    E = 520.0kJ / mol,
                    V = 22.0m^3 / mol,
                    r = 1.2NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Quartz Diorite rheology
        (
            "Quartz Diorite | Hansen & Carter (1982)",
            #  After Hansen (1982), 'Semibrittle creep of selected crustal rocks at 1000 MPa.' and, Hansen & Carter (1982),
            #  'Flow properties of continental lithosphere.'
            #  Hansen (1982), Fig. 53, page 184 in PDF viewer and table 18, page 224
            #  Carter & Tsenn (1986), table 4, page 18 in PDF viewer
            (
                DislocationCreep(;
                    Name = "Quartz Diorite | Hansen & Carter (1982)",
                    n = 2.25NoUnits,
                    A = 3.5e-2MPa^(-9 // 4) / s,
                    E = 212kJ / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Verified with the original publication Hansen(1982)(NM). Values checked, plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Diabase rheology
        (
            "Diabase | Caristan (1982)",
            #  After Caristan (1982), 'The transition from high temperature creep to fracture in Maryland diabase.'
            #  and, Bremond (1999),
            #  'Hydrothermalism and diapirism in the Archean: gravitational instability constraints'
            #  Bremond (1999), page 5 in text
            (
                DislocationCreep(;
                    Name = "Diabase | Caristan (1982)",
                    n = 3.05NoUnits,
                    A = 6.0e-2MPa^(-61 // 20) / s,
                    E = 276kJ / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Bremond (1999)), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Tumut Pond Serpentinite rheology
        (
            "Tumut Pond Serpentinite | Raleigh and Paterson (1965)",
            #  After Raleigh and Paterson (1965), 'Experimental deformation of serpentinite and its tectonic implications'
            #  and, Bremond (1999),
            #  'Hydrothermalism and diapirism in the Archean: gravitational instability constraints'
            #  Bremond (1999), page 5 in text
            (
                DislocationCreep(;
                    Name = "Tumut Pond Serpentinite | Raleigh and Paterson (1965)",
                    n = 2.8NoUnits,
                    A = 6.3e-7MPa^(-14 // 5) / s,
                    E = 66kJ / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Bremond (1999)), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Maryland strong diabase rheology
        (
            "Maryland strong diabase | Mackwell et al. (1998)",
            #  After Mackwell et al. (1998), 'High-temperature deformation of dry diabase with application to tectonics on Venus'
            #  Mackwell et al. (1998), page 980, equation in text
            (
                DislocationCreep(;
                    Name = "Maryland strong diabase | Mackwell et al. (1998)",
                    n = 4.7NoUnits,
                    A = 8MPa^(-47 // 10) / s,
                    E = 485kJ / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Mackwell et al. (1998))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Mafic Granulite rheology
        (
            "Mafic Granulite | Wilks and Carter (1990)",
            #  After Li, Gerya and Burg (2010), table 2
            #  referring to Ranalli (1995), 'Rheology of the Earth' (Book), page 334, table 10.3
            #  referring to Wilks and Carter (1990), 'Rheology of some continental lower crustal rocks', Fig. 6, Pikwitonei Granulite
            (
                DislocationCreep(;
                    Name = "Mafic Granulite | Wilks and Carter (1990)",
                    n = 4.2NoUnits,
                    A = 1.4e4MPa^(-21 // 5) / s,
                    E = 445kJ / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Wilks and Carter (1990))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Wet Quartzite rheology
        (
            "Wet Quartzite | Ueda et al. (2008)",
            #  Ueda et al. (2008), table 1
            (
                DislocationCreep(;
                    Name = "Wet Quartzite | Ueda et al. (2008)",
                    n = 2.3NoUnits,
                    A = 1 * exp10(-3.5)MPa^(-23 // 10) / s,
                    E = 154kJ / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Ueda et al. (2008))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Granite rheology
        (
            "Granite | Carter and Tsenn (1987)",
            #  Huismans et al. (2001), table 2
            #  referring to Carter and Tsenn (1987), 'Flow properties of continental lithosphere', table 4, Westerly Granite (dry)
            #  referring to Hansen and Carter (1983), 'Semibrittle Creep Of Dry And Wet Westerly Granite At 1000 MPa', not accessible
            (
                DislocationCreep(;
                    Name = "Granite | Carter and Tsenn (1987)",
                    n = 3.3NoUnits,
                    A = 1.0 * exp10(-5.7)MPa^(-33 // 10) / s,
                    E = 186.5kJ / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Carter and Tsenn (1987))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Plagioclase An75 rheology
        (
            "Plagioclase An75 | Ji and Zhao (1993)",
            #  Ranalli (1995), page 334, table 10.3
            #  referring to Ji and Zhao (1993), 'Flow laws of multiphase rocks calculated from experimental data on the constituent phases', table 2 , plagioclase (Ab25An75)
            #  referring to Shelton and Tullis (1981), 'Experimental flow laws for crustal rocks', not accessible
            (
                DislocationCreep(;
                    Name = "Plagioclase An75 | Ji and Zhao (1993)",
                    n = 3.2NoUnits,
                    A = 3.27e-4MPa^(-16 // 5) / s,
                    E = 238kJ / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Ji and Zhao (1993))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Dry Anorthite rheology
        (
            "Dry Anorthite | Rybacki et al. (2006)",
            # Rybacki, Gottschalk, Wirth and Dresen (2006), table 5
            (
                DislocationCreep(;
                    Name = "Dry Anorthite | Rybacki et al. (2006)",
                    n = 3.0NoUnits,
                    A = exp10(12.7)MPa^(-3) / s,
                    E = 641kJ / mol,
                    V = 24.0e-6m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values taken from Bürgmann & Rybacki (2008) Supplementary table 1. (BK), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Wet Anorthite rheology
        (
            "Wet Anorthite | Rybacki et al. (2006)",
            # Rybacki, Gottschalk, Wirth and Dresen (2006), table 5
            (
                DislocationCreep(;
                    Name = "Wet Anorthite | Rybacki et al. (2006)",
                    n = 3.0NoUnits,
                    A = (10^0.2)MPa^(-3) / s,
                    E = 345kJ / mol,
                    V = 38m^3 / mol,
                    r = 1NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Rybacki, Gottschalk, Wirth and Dresen (2006))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Wet Olivine rheology
        (
            "Wet Olivine | Hirth and Kohlstedt (2003)",
            #  Hirth and Kohlstedt (2003), table 1, no constant C_OH
            (
                DislocationCreep(;
                    Name = "Wet Olivine | Hirth and Kohlstedt (2003)",
                    n = 3.5NoUnits,
                    A = 1600.0MPa^(-7 // 2) / s,
                    E = 520kJ / mol,
                    V = 11.0e-6m^3 / mol,
                    r = 1.2NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Hirth and Kohlstedt (2003))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Wet Quarzite rheology
        (
            "Wet Quarzite | Kirby (1983)",
            #  Li, Gerya and Burg (2010), table 2
            #  Ranalli (1995), 'Rheology of the Earth' (Book)
            #  referring to Kirby (1983), table 2, first quartzite (wet)
            #  referring to Koch et al. (1981), unpublished manuscript...
            (
                DislocationCreep(;
                    Name = "Wet Quarzite | Kirby (1983)",
                    n = 2.3NoUnits,
                    A = uconvert(Pa^(-23 // 10) / s, 3.2e-4MPa^(-23 // 10) / s),
                    E = 154.0e3J / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Kirby (1983))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Wet Upper Mantle Olivine rheology
        (
            "Wet Upper Mantle Olivine | Afonso and Ranalli (2004)",
            #  Schmalholz, Kaus, Burg (2009), table 1
            #  referring to Afonso and Ranalli (2004), table 1, wet peridotite
            #  referring to Chopra and Paterson papers, but values dont fit the Afonso and Ranalli (2004) ones
            (
                DislocationCreep(;
                    Name = "Wet Upper Mantle Olivine | Afonso and Ranalli (2004)",
                    n = 4.0NoUnits,
                    A = 2.0e3MPa^(-4) / s,
                    E = 471kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Afonso and Ranalli (2004))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Granite rheology
        (
            "Granite | Tirel et al. (2008)",
            #  Tirel et al. (2008), table 1
            #  referring to Kirby and Kronenberg (1987), table 3
            #  different values for n and A in Kirby and Kronenberg (1987) compared with Tirel et al. (2008)
            (
                DislocationCreep(;
                    Name = "Granite | Tirel et al. (2008)",
                    n = 3.2NoUnits,
                    A = 1.25e-9MPa^(-3.2) / s,
                    E = 123kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (Tirel et al. (2008))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Dry olivine rheology
        (
            "Dry Olivine | Gerya (2019)",
            # after Exercise 6.1 of Numerical Geodynamics
            (
                DislocationCreep(;
                    Name = "Dry Olivine | Gerya (2019)",
                    n = 3.5NoUnits,
                    #A = 2.5e-17Pa^(-7//2)/s,
                    A = uconvert(MPa^(-7 // 2) / s, 2.5e-17Pa^(-7 // 2) / s),
                    E = 532.0kJ / mol,
                    V = 0.0m^3 / mol,
                    Apparatus = AxialCompression,   # used in book (according to matlab script)
                    r = 0.0NoUnits,
                ),
                MaterialParamsInfo(;
                    Comment = "This is from Exercise 6.1, indicated to be valid for olivine under upper mantle conditions.",
                    BibTex_Reference = "
                        @book{Gerya_2019, 
                        title={Introduction to Numerical Geodynamic Modelling}, 
                        ISBN={978-1-107-14314-2}, 
                        note={Google-Books-ID: 8XGSDwAAQBAJ}, 
                        publisher={Cambridge University Press}, 
                        author={Gerya, Taras}, 
                        year={2019}, month={May}, 
                        language={en} }
                ",
                ),
            ),
        )

        # Rock salt rheology
        (
            "Rock salt | Li & Urai (2016)",
            #  Li & Urai (2016), table 1
            #  refers to Wawersik & Zeuch (1986), values can not be reproduced!
            (
                DislocationCreep(;
                    Name = "Rock salt | Li & Urai (2016)",
                    n = 5.0NoUnits,
                    A = 7.26e-6MPa^(-5) / s,
                    E = 53.92kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked in (Li & Urai, (2016)) are different from given source (Wawersik & Zeuch, (1986))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Avery Island rock salt rheology
        (
            "Salado rock salt | Li & Urai (2016)",
            #  Li & Urai (2016), table 1
            #  refers to Wawersik & Zeuch (1986), values can not be reproduced!
            (
                DislocationCreep(;
                    Name = "Salado Rock salt | Li & Urai (2016)",
                    n = 5.0NoUnits,
                    A = 7.26e-6MPa^(-5) / s,
                    E = 53.92kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked in (Li & Urai, (2016)) are different from given source (Wawersik & Zeuch, (1986))(NM), plots are not reproduced (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Wet olivine rheology
        (
            "Wet Olivine | Mei & Kohlstedt (2000b)",
            #  Mei & Kohlstedt (2000b), table 1
            (
                DislocationCreep(;
                    Name = "Wet Olivine | Mei & Kohlstedt (2000b)",
                    n = 3.0NoUnits,
                    A = (10^3.2)MPa^(-3) / s,
                    E = 470.0kJ / mol,
                    V = 20.0m^3 / mol,
                    r = 0.98NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "The A value is not exactly the same as in Mei & Kohlstedt (2000b) but is approximated (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Dry olivine rheology
        (
            "Dry Olivine | Karato & Jung (2003)",
            #  Mei & Kohlstedt (2000b), abstract
            (
                DislocationCreep(;
                    Name = "Dry Olivine | Karato & Jung (2003)",
                    n = 3.0NoUnits,
                    A = (10^6.1)MPa^(-3) / s,
                    E = 510.0kJ / mol,
                    V = 14.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        # Wet olivine rheology
        (
            "Wet Olivine | Karato & Jung (2003)",
            #  Karato & Jung (2003), abstract
            (
                DislocationCreep(;
                    Name = "Wet Olivine | Karato & Jung (2003)",
                    n = 3.0NoUnits,
                    A = (10^2.9)MPa^(-3) / s,
                    E = 510.0kJ / mol,
                    V = 24.0m^3 / mol,
                    r = 1.2NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Clinopyroxene rheology
        (
            "Wet Clinopyroxene | Chen et al. (2006)",
            #  Chen et al. (2006), section 4. Discussion, equation (3)
            (
                DislocationCreep(;
                    Name = "Wet Clinopyroxene | Chen et al. (2006)",
                    n = 2.7NoUnits,
                    A = (10^6.7)MPa^(-27 // 10) / s,
                    E = 670.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 3.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Dry Clinopyroxene rheology
        (
            "Dry Clinopyroxene | Bystricky & Mackwell (2001)",
            #  Bystricky & Mackwell (2001), section 4. Discussion, equation (3)
            (
                DislocationCreep(;
                    Name = "Dry Clinopyroxene | Bystricky & Mackwell (2001)",
                    n = 4.7NoUnits,
                    A = (10^9.8)MPa^(-47 // 10) / s,
                    E = 760.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Dry Clinopyroxene rheology
        (
            "Dry Clinopyroxene | Bystricky & Mackwell (2001)",
            #  Bystricky & Mackwell (2001), section 4. Discussion, equation (3)
            (
                DislocationCreep(;
                    Name = "Dry Clinopyroxene | Bystricky & Mackwell (2001)",
                    n = 4.7NoUnits,
                    A = (10^10.8)MPa^(-47 // 10) / s,
                    E = 760.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Dry Diopside rheology
        (
            "Dry Diopside | Dimanov & Dresen (2005)",
            #  Dimanov & Dresen (2005), table 3b
            (
                DislocationCreep(;
                    Name = "Dry Diopside | Dimanov & Dresen (2005)",
                    n = 5.5NoUnits,
                    A = uconvert(MPa^(-55 // 10) / s, 3.01e-28Pa^(-55 // 10) / s),
                    E = 691.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Diopside rheology
        (
            "Wet Diopside | Dimanov & Dresen (2005)",
            #  Dimanov & Dresen (2005), table 3b
            (
                DislocationCreep(;
                    Name = "Wet Diopside | Dimanov & Dresen (2005)",
                    n = 5.5NoUnits,
                    A = uconvert(MPa^(-55 // 10) / s, 5.16e-33Pa^(-55 // 10) / s),
                    E = 534.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Omphacite rheology
        (
            "Wet Omphacite | Zhang et al. (2006)",
            #  Zhang et al. (2006), equation (4)
            (
                DislocationCreep(;
                    Name = "Wet Omphacite | Zhang et al. (2006)",
                    n = 3.5NoUnits,
                    A = (10^-2)MPa^(-7 // 2) / s,
                    E = 310.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Jadeit rheology
        (
            "Wet Jadeit | Orzol et al. (2006)",
            #  Orzol et al. (2006), page 11
            (
                DislocationCreep(;
                    Name = "Wet Jadeit | Orzol et al. (2006)",
                    n = 3.7NoUnits,
                    A = (10^-3.3)MPa^(-37 // 10) / s,
                    E = 326.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Dry Anorthite rheology
        (
            "Dry Anorthite | Rybacki & Dresen (2000)",
            #  Rybacki & Dresen (2000), table 5
            (
                DislocationCreep(;
                    Name = "Dry Anorthite | Rybacki & Dresen (2000)",
                    n = 3.0NoUnits,
                    A = (10^12.7)MPa^(-3) / s,
                    E = 648.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Anorthite rheology
        (
            "Wet Anorthite | Rybacki & Dresen (2000)",
            #  Rybacki & Dresen (2000), table 5
            (
                DislocationCreep(;
                    Name = "Wet Anorthite | Rybacki & Dresen (2000)",
                    n = 3.0NoUnits,
                    A = (10^0.2)MPa^(-3) / s,
                    E = 356.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Quartzite rheology
        (
            "Wet Quartzite | Rutter & Brodie (2004)",
            #  Rutter & Brodie (2004), table 5
            (
                DislocationCreep(;
                    Name = "Wet Quartzite | Rutter & Brodie (2004)",
                    n = 3.0NoUnits,
                    A = (10^-4.9)MPa^(-3) / s,
                    E = 242.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 1.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Quartzite rheology
        (
            "Wet Quartzite | Hirth et al. (2001)",
            #  Hirth et al. (2001), table 5
            (
                DislocationCreep(;
                    Name = "Wet Quartzite | Hirth et al. (2001)",
                    n = 4.0NoUnits,
                    A = (10^-11.2)MPa^(-4) / s,
                    E = 135.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 1.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Dry Quartzite rheology
        (
            "Dry Quartzite | Jaoul et al. (1984)",
            #  Jaoul et al. (1984), table 1, first entry
            (
                DislocationCreep(;
                    Name = "Dry Quartzite | Jaoul et al. (1984)",
                    n = 2.8NoUnits,
                    A = (10^-5.415)MPa^(-14 // 5) / s,
                    E = 184.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Quartzite rheology
        (
            "Wet Quartzite | Jaoul et al. (1984)",
            #  Jaoul et al. (1984), table 1, second entry
            (
                DislocationCreep(;
                    Name = "Wet Quartzite | Jaoul et al. (1984)",
                    n = 2.8NoUnits,
                    A = (10^-5.045)MPa^(-14 // 5) / s,
                    E = 163.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 0.0NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Quartzite rheology
        (
            "Wet Quartzite | Tokle et al. (2019)",
            #  Tokle et al. (2019), table 1, 2nd extrapolated fit
            (
                DislocationCreep(;
                    Name = "Wet Quartzite | Tokle et al. (2019)",
                    n = 3.0NoUnits,
                    A = (10^-11.959)MPa^(-3) / s,
                    E = 115.0kJ / mol,
                    V = 0.0m^3 / mol,
                    r = 1.2NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet Quartzite rheology
        (
            "Wet Quartzite | Lu and Jiang (2019)",
            #  Lu and Jiang (2019), section 3, equation (6)
            (
                DislocationCreep(;
                    Name = "Wet Quartzite | Lu and Jiang (2019)",
                    n = 3.0NoUnits,
                    A = (10^-14.2218)MPa^(-3) / s,
                    E = 132.0kJ / mol,
                    V = 35.3e-6m^3 / mol,
                    r = 2.7NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #"Low pressure" wet quartzite rheology
        (
            "low pressure wet Quartzite | Lusk et al. (2021)",
            #  Lusk et al. (2021), abstract, 1st law
            (
                DislocationCreep(;
                    Name = "low pressure wet Quartzite | Lusk et al. (2021)",
                    n = 3.5NoUnits,
                    A = (10^-9.3)MPa^(-7 // 2) / s,
                    E = 118.0kJ / mol,
                    V = 2.59e-6m^3 / mol,
                    r = 0.49NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #"High pressure" wet quartzite rheology
        (
            "high pressure wet Quartzite | Lusk et al. (2021)",
            #  Lusk et al. (2021), abstract, 2nd law
            (
                DislocationCreep(;
                    Name = "high pressure wet Quartzite | Lusk et al. (2021)",
                    n = 2.1NoUnits,
                    A = (10^-6.36)MPa^(-21 // 10) / s,
                    E = 94.0kJ / mol,
                    V = 1.44e-6m^3 / mol,
                    r = 0.2NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )

        #Wet quartzite rheology
        (
            "Wet Quartzite | Lusk et al. (2021)",
            #  Lusk et al. (2021), table 2, full data set
            (
                DislocationCreep(;
                    Name = "Wet Quartzite | Lusk et al. (2021)",
                    n = 2.0NoUnits,
                    A = (10^-7.9)MPa^(-2) / s,
                    E = 77.0kJ / mol,
                    V = 2.59e-6m^3 / mol,
                    r = 0.49NoUnits,
                    Apparatus = AxialCompression,
                ),
                MaterialParamsInfo(;
                    Comment = "Values checked (NM).",
                    BibTex_Reference = "
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
                ),
            ),
        )
    ],
); # end of setting pre-defined creep laws
