module Diffusion

using GeoParams
export diffusion_database, diffusion_database_info, diffusion_law_list

function diffusion_law_list()
    m = @__MODULE__
    out = string.(names(m; all = true, imported = true))
    wet = Symbol.(filter(x -> startswith(x, "wet"), out))
    dry = Symbol.(filter(x -> startswith(x, "dry"), out))
    all_names = vcat(wet, dry)
    return [getfield(m, fun) for fun in all_names]
end

function TestDiff()
    data = DiffusionCreep(;
        Name = "Test Diff",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = 2.070729911135297e-7MPa^(-1) * m^3 * s^(-1),    # material specific rheological parameter
        E = 375.0kJ / mol,                      # activation energy
        V = 6.0e-6m^3 / mol,                    # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Law to check Thorsten Beckers book rheology and its conversion of the A factor",
        BibTex_Reference = "",
    )
    return data, info
end

function dry_anorthite_Rybacki_2006()
    # Dry Plagioclase rheology
    data = DiffusionCreep(;
        Name = "Dry Anorthite | Rybacki et al. (2006)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^12.1)MPa^(-1) * μm^3 * s^(-1),  # material specific rheological parameter
        E = 460.0kJ / mol,                      # activation energy
        V = 24.0e-6m^3 / mol,                     # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values; not yet plots (BK)",
        BibTex_Reference = "
            @article{Rybacki_Gottschalk_Wirth_Dresen_2006, 
            title={Influence of water fugacity and activation volume on the flow properties of fine-grained anorthite aggregates}, 
            volume={111}, 
            DOI={10.1029/2005JB003663}, 
            number={B3}, 
            journal={Journal of Geophysical Research: Solid Earth}, 
            author={Rybacki, E. and Gottschalk, M. and Wirth, R. and Dresen, G.}, 
            year={2006}, 
            month={Mar}
            }
    ",
    )
    return data, info
end
function wet_olivine_Mei_2000a()
    data = DiffusionCreep(;
        Name = "Wet Olivine | Mei & Kohlstedt (2000a)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 1.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^4.7)MPa^(-1) * μm^3 * s^(-1), # material specific rheological parameter
        E = 295.0kJ / mol,                      # activation energy
        V = 20.0e-6m^3 / mol,                     # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values are not exactly the same as in Mei & Kohlstedt (2000a) but are close to them (NM)",
        BibTex_Reference = "
            @article{mei2000influence,
            title={Influence of water on plastic deformation of olivine aggregates: 1. Diffusion creep regime},
            author={Mei, S and Kohlstedt, DL},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={105},
            number={B9},
            pages={21457--21469},
            year={2000},
            publisher={Wiley Online Library}
            }
    ",
    )

    return data, info
end
function wet_olivine_Hirth_2003()
    data = DiffusionCreep(;
        Name = "Wet Olivine | Hirth & Kohlstedt (2003)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 1.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^7.4)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 375.0kJ / mol,                        # activation energy
        V = 20.0e-6m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "The A value is not exactly the same as in Hirth & Kohlstedt (2003) but is approximated (NM)",
        BibTex_Reference = "
            @article{hirth2004rheology,
            title={Rheology of the Upper Mantle and the Mantle Wedge: A View from the Experimentalists},
            author={Hirth, Greg and Kohlstedt, David},
            journal={Inside the Subduction Factory},
            volume={138},
            pages={83--105},
            year={2004},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end
function dry_olivine_Hirth_2003()
    data = DiffusionCreep(;
        Name = "Dry Olivine | Hirth & Kohlstedt (2003)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^9.2)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 375.0kJ / mol,                        # activation energy
        V = 10.0e-6m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "The A value is not exactly the same as in Hirth & Kohlstedt (2003) but is approximated (NM)",
        BibTex_Reference = "
            @article{hirth2004rheology,
            title={Rheology of the Upper Mantle and the Mantle Wedge: A View from the Experimentalists},
            author={Hirth, Greg and Kohlstedt, David},
            journal={Inside the Subduction Factory},
            volume={138},
            pages={83--105},
            year={2004},
            publisher={Wiley Online Library}
            }
    ",
    )
    return data, info
end
function wet_olivine_Faul_2006()
    data = DiffusionCreep(;
        Name = "Dry Olivine | Faul & Jackson (2006)",
        n = 1.4NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^10.3)MPa^(-7 // 5) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 484.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Paper not found, taking values from Bürgmannn & Dresen (2008), Supplementary Material, Table 1 (NM)",
        BibTex_Reference = "
            @article{burgmann2008rheology,
            title={Rheology of the lower crust and upper mantle: Evidence from rock mechanics, geodesy, and field observations},
            author={B{\"u}rgmann, Roland and Dresen, Georg},
            journal={Annu. Rev. Earth Planet. Sci.},
            volume={36},
            pages={531--567},
            year={2008},
            publisher={Annual Reviews}
            }
    ",
    )
    return data, info
end
function dry_clinopyroxene_HierMajumder_2005()
    data = DiffusionCreep(;
        Name = "Dry Clinopyroxene | Hier-Majumder et al. (2005)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = ((10^25.3) / 64.9e3)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 760.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values checked (NM), dividend at A is the shear modulus of diorite (see paper)",
        BibTex_Reference = "
            @article{hier2005water,
            title={Water weakening of clinopyroxenite in diffusion creep},
            author={Hier-Majumder, Saswata and Mei, Shenghua and Kohlstedt, David L},
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
function wet_clinopyroxene_HierMajumder2005()
    data = DiffusionCreep(;
        Name = "Wet Clinopyroxene | Hier-Majumder et al. (2005)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 1.4NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^7.9 / 64.9e3)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 340.0kJ / mol,                        # activation energy
        V = 14.0e-6m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values checked (NM), dividend at A is the shear modulus of diorite (see paper)",
        BibTex_Reference = "
            @article{hier2005water,
            title={Water weakening of clinopyroxenite in diffusion creep},
            author={Hier-Majumder, Saswata and Mei, Shenghua and Kohlstedt, David L},
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
function dry_clinopyroxene_Bystricky_2001()
    data = DiffusionCreep(;
        Name = "Dry Clinopyroxene | Bystricky & Mackwell (2001)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^15.1)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 560.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values checked (NM)",
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
    )
    return data, info
end

function dry_diopside_Dimanov_2005()
    data = DiffusionCreep(;
        Name = "Dry Diopside | Dimanov & Dresen (2005)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^14)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 528.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values according to Bürgmann & Dresen (2008), Supplementary (NM), cant reproduce Dimanov & Dresen (2005) A value",
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
    )
    return data, info
end

function wet_diopside_Dimanov_2005()
    data = DiffusionCreep(;
        Name = "Wet Diopside | Dimanov & Dresen (2005)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^8.1)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 528.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values according to Bürgmann & Dresen (2008), Supplementary (NM), cant reproduce/find corresponding Dimanov & Dresen (2005) values",
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
    )
    return data, info
end

function wet_clinopyroxene_Chen_2006()
    data = DiffusionCreep(;
        Name = "Wet Clinopyroxene | Chen et al. (2006)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^8.1)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 528.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
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
    )
    return data, info
end

function dry_anorthite_Rybacki_2000()
    data = DiffusionCreep(;
        Name = "Dry Anorthite | Rybacki & Dresen (2000)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^12.1)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 467.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values checked (NM)",
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
    )
    return data, info
end

function wet_anorthite_rybacki_2000()
    data = DiffusionCreep(;
        Name = "Wet Anorthite | Rybacki & Dresen (2000)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^1.7)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 170.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values checked (NM)",
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
    )
    return data, info
end

function wet_anorthite_rybacki_2006()
    data = DiffusionCreep(;
        Name = "Wet Anorthite | Rybacki et al. (2006)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 1.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = (10^-0.7)MPa^(-1) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 159.0kJ / mol,                        # activation energy
        V = 38.0e-6m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values checked (NM)",
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
    )
    return data, info
end

function wet_quartzite_Rutter_2004()
    data = DiffusionCreep(;
        Name = "Wet Quartzite | Rutter & Brodie (2004)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -2.0NoUnits,                        # grain size exponent
        A = (10^-0.4)MPa^(-1) * μm^2 * s^(-1),    # material specific rheological parameter
        E = 220.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Values checked (NM)",
        BibTex_Reference = "
            @article{rutter2004experimental,
            title={Experimental grain size-sensitive flow of hot-pressed Brazilian quartz aggregates},
            author={Rutter, EH and Brodie, KH},
            journal={Journal of Structural Geology},
            volume={26},
            number={11},
            pages={2011--2023},
            year={2004},
            publisher={Elsevier}
            }
    ",
    )
    return data, info
end

@inline diffusion_database(f::F) where {F} = first(f())
@inline diffusion_database_info(f::F) where {F} = last(f())

end
