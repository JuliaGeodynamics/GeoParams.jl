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
NOTE: the printed law is 2.5× too slow vs the paper's own Fig. 2a and data; a data-consistent version is `dry_olivine_Hirth_2003_refit`.
Dislocation creep data for dry olivine after Hirth, G. & Kohlstedt (2003)
"""
function dry_olivine_Hirth_2003()
    # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    # Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters

    data = DislocationCreep(;
        Name = "Dry Olivine | Hirth & Kohlstedt (2003)",
        n = 3.5NoUnits,
        r = 0.0NoUnits,
        A = 1.1e5MPa^(-7 // 2) / s,
        E = 530.0kJ / mol,
        V = 14.0e-6m^3 / mol,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: the printed law is 2.5× too slow vs the paper's own data: -0.40 log10 units vs the Fig. 2a dislocation-creep fit line, -0.30 ± 0.19 vs the Mei & Kohlstedt (2000) dry data (N = 12, diffusion removed) and -0.34/-0.48 vs Karato & Jung (2003), Table 1; a data-consistent version is dry_olivine_Hirth_2003_refit. Hirth & Kohlstedt (2003), Table 1, row 'dry dislocation': A = 1.1e5 MPa^-n/s, n = 3.5, E = 530; V = 14 cm^3/mol from Table 2 (range 13.4-27).",
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
    )
    return data, info
end

"""
NOTE: refit of `dry_olivine_Hirth_2003` (A = 2.8e5 instead of the printed 1.1e5 MPa^-n/s) so that the law reproduces the paper's own Fig. 2a and data.
Dislocation creep data for dry olivine after Hirth, G. & Kohlstedt (2003)
"""
function dry_olivine_Hirth_2003_refit()
    # after Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    # Inside the subduction Factory 83?105. Table 1, "dry dislocation" parameters
    # A refitted to the paper's own data (see Comment); the printed value is in dry_olivine_Hirth_2003

    data = DislocationCreep(;
        Name = "Dry Olivine (refit to data) | Hirth & Kohlstedt (2003)",
        n = 3.5NoUnits,
        r = 0.0NoUnits,
        A = 2.8e5MPa^(-7 // 2) / s,
        E = 530.0kJ / mol,
        V = 14.0e-6m^3 / mol,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: refit of dry_olivine_Hirth_2003 with A = 2.8e5 MPa^-n/s (= 1.1e5 × 10^0.40; not a printed value) so that the law reproduces the paper's own data. The printed A = 1.1e5 is -0.40 log10 units (2.5× too slow) vs the Hirth & Kohlstedt (2003) Fig. 2a dislocation-creep fit line (1250 °C, 300 MPa, σ = 89-250 MPa), -0.30 ± 0.19 vs the Mei & Kohlstedt (2000) dry data it was fitted to (N = 12, σ ≥ 150 MPa, diffusion creep removed) and -0.34/-0.48 vs the dry strengths of Karato & Jung (2003), Table 1. n = 3.5, E = 530 and V = 14 cm^3/mol as in dry_olivine_Hirth_2003, which holds the printed law (Table 1, row 'dry dislocation').",
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
    )
    return data, info
end

"""
Dislocation creep data for wet olivine (1) at constant water content after Hirth, G. & Kohlstedt (2003), Table 1,
row 'wet dislocation (constant C_OH)'. Here `f` is the water content C_OH in H/10^6 Si (dimensionless), not a
fugacity; A is therefore stored in MPa^-n only.
"""
function wet_olivine1_Hirth_2003()
    # After Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists.
    #   Inside the subduction Factory 83?105. Table 1, "wet dislocation (constant C_OH)" parameters
    #  f = C_OH in H/10^6 Si must be passed explicitly (HK03 use C_OH = 1000 in their example, footnote e)

    data = DislocationCreep(;
        Name = "Wet Olivine 1 | Hirth & Kohlstedt (2003)",
        n = 3.5NoUnits,
        A = 90MPa^(-7 // 2) / s,
        E = 480kJ / mol,
        V = 11.0e-6m^3 / mol,
        r = 1.2NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'wet dislocation (constant C_OH)': A = 90, n = 3.5, r = 1.2, E = 480, V = 11 (footnote e example: 2.5e-12 /s). Pass f = C_OH in H/10^6 Si. Reproduces the paper's data: log10 misfit +0.25 ± 0.19 (N = 34, Mei & Kohlstedt 2000 wet data, σ ≥ 120 MPa); -0.01 vs the Fig. 5b line.",
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
    #  Carter & Tsenn (1987), table 4, page 18 in PDF viewer, list n = 2.4, log A = -1.5 (MPa^-n/s), Q = 212 instead

    data = DislocationCreep(;
        Name = "Quartz Diorite | Hansen & Carter (1982)",
        n = 2.25NoUnits,
        A = 3.5e-2MPa^(-9 // 4) / s,
        E = 212kJ / mol,
        V = 0m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hansen (1982, PhD thesis), table 18 (n = 2.25, A = 3.5e-2, Q = 212). Carter & Tsenn (1987), Table 4, list n = 2.4, log A = -1.5 for the same material. Kirby (1983) and Tirel et al. (2008) quote log A = -2.9, Q = 219 (30-55× slower); the experimental ranges of Carter & Tsenn (1987) weakly favour the stored values.",
        BibTex_Reference = "
            @inproceedings{hansen1982creep,
            title={Creep of selected crustal rocks at 1000 MPa},
            author={Hansen, Francis D and Carter, Neville L},
            booktitle={Issues in Rock Mechanics, Proceedings of the 23rd U.S. Symposium on Rock Mechanics},
            address={Berkeley, CA},
            year={1982}
            }
            @phdthesis{hansen1982semibrittle,
            title={Semibrittle creep of selected crustal rocks at 1000 MPa},
            author={Hansen, Francis D},
            school={Texas A\\&M University},
            year={1982}
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
        Name = "Diabase | Caristan (1982)",
        n = 3.05NoUnits,
        A = 6.12e-2MPa^(-61 // 20) / s,
        E = 276kJ / mol,
        V = 0m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Maryland diabase. Caristan (1982), p. 6789 (text beside Fig. 10): ε̇ = A σ^n exp(-Q/RT) with A = 6.12e-2 MPa^-n (s^-1), n = 3.05, Q = 276 kJ/mol, cited from Caristan's thesis; experiments in an argon gas-medium apparatus at 1000 °C, as-received samples, differential stress. Matches Carter & Tsenn (1987), Table 4 (n = 3.05, log A = -1.2 MPa^-n/s, Q = 276 kJ/mol, triaxial compression, 900-1000 °C; listed as dry), and Bremond d'Ars et al. (1999).",
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
    )
    return data, info
end

function serpentinite_Raleigh_1965()
    #  After Raleigh and Paterson (1965), 'Experimental deformation of serpentinite and its tectonic implications'
    #  and, Bremond (1999),
    #  'Hydrothermalism and diapirism in the Archean: gravitational instability constraints'
    #  Bremond (1999), page 5 in text

    data = DislocationCreep(;
        Name = "Tumut Pond Serpentinite | Raleigh and Paterson (1965)",
        n = 2.8NoUnits,
        A = 6.3e-7MPa^(-14 // 5) / s,
        E = 66kJ / mol,
        V = 0m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Not a law from Raleigh & Paterson (1965): the paper reports constant strain-rate strength and ductility tests (up to 700 °C, 1-5 kbar) and dehydration embrittlement, but no power-law fit, activation energy or stress exponent. The values (n = 2.8, A = 6.3e-7 MPa^-n/s, E = 66 kJ/mol) are a later fit taken from Bremond d'Ars et al. (1999), p. 5, which was not available; their origin is unverified. The law agrees with the paper's strengths (at 7e-4 s^-1) only at about 400-475 °C.",
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
    )
    return data, info
end

function strong_diabase_Mackwell_1998()
    #  After Mackwell et al. (1998), 'High-temperature deformation of dry diabase with application to tectonics on Venus'
    #  Mackwell et al. (1998), page 980, flow law for Maryland diabase in the text
    data = DislocationCreep(;
        Name = "Maryland strong diabase | Mackwell et al. (1998)",
        n = 4.7NoUnits,
        A = 8MPa^(-47 // 10) / s,
        E = 485kJ / mol,
        V = 0m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Mackwell et al. (1998), p. 980, flow law for Maryland diabase: ε̇ = (8 ± 4) σ^(4.7 ± 0.6) exp(-(485 ± 30 kJ/mol)/RT), σ (differential) in MPa, A in MPa^-4.7 s^-1; dry, gas-medium creep tests. The weaker Columbia diabase (same page): A = 190 ± 110 MPa^-4.7 s^-1, same n and Q. Reproduces the paper's data: log10 misfit -0.02 ± 0.20 (N = 57, Table 3).",
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
    )
    return data, info
end

function mafic_granulite_Wilks_1990()
    #  After Li, Gerya and Burg (2010), table 2
    #  referring to Ranalli (1995), 'Rheology of the Earth' (Book), page 334, table 10.3
    #  referring to Wilks and Carter (1990), 'Rheology of some continental lower crustal rocks', Fig. 6, Pikwitonei Granulite

    data = DislocationCreep(;
        Name = "Mafic Granulite | Wilks and Carter (1990)",
        n = 4.2NoUnits,
        A = 1.4e4MPa^(-21 // 5) / s,
        E = 445kJ / mol,
        V = 0m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Wilks & Carter (1990), Fig. 6B, Pikwitonei transitional (mafic) granulite: ε̇ = 1.4e4 exp(-445/RT) σ^4.2, A in MPa^-n s^-1 (section 'Steady-state flow'), σ = σ1 - σ3 in MPa, Q in kJ/mol; Griggs-Blacic solid-medium apparatus, 0.8-1.0 GPa. The authors call it preliminary (strain softening above ~800 °C). Same values in Ranalli (1995), Table 10.3, and Afonso & Ranalli (2004), Table 1. Reproduces the paper's data: log10 misfit 0.00 ± 0.24 (N = 10, Table 2, steady-state runs).",
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
    )
    return data, info
end

function wet_quartzite_Ueda_2008()
    #  Ueda et al. (2008), table 1
    data = DislocationCreep(;
        Name = "Wet Quartzite | Ueda et al. (2008)",
        n = 2.3NoUnits,
        A = 1 * exp10(-3.5)MPa^(-23 // 10) / s,
        E = 154kJ / mol,
        V = 0m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Ranalli (1995) wet quartzite (A = 3.2e-4 MPa^-n/s, n = 2.3, E = 154) as listed in Gerya-group tables; Ueda et al. (2008) Table 1 not seen.",
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
    )
    return data, info
end

function granite_Carter_1987()
    #  Huismans et al. (2001), table 2
    #  referring to Carter and Tsenn (1987), 'Flow properties of continental lithosphere', table 4, Westerly Granite (dry)
    #  referring to Hansen and Carter (1983), 'Semibrittle Creep Of Dry And Wet Westerly Granite At 1000 MPa', not accessible
    data = DislocationCreep(;
        Name = "Granite | Carter and Tsenn (1987)",
        n = 3.3NoUnits,
        A = 1.0 * exp10(-5.7)MPa^(-33 // 10) / s,
        E = 186.5kJ / mol,
        V = 0m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Carter & Tsenn (1987), Table 4, 'Westerly granite, dry (Hansen & Carter, 1983)': log A = -5.7 MPa^-n/s, n = 3.3, Q = 186.5.",
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
    )
    return data, info
end

function plagioclase_An75_Ji_1993()
    #  Ranalli (1995), page 334, table 10.3
    #  referring to Ji and Zhao (1993), 'Flow laws of multiphase rocks calculated from experimental data on the constituent phases', table 2 , plagioclase (Ab25An75)
    #  referring to Shelton and Tullis (1981), 'Experimental flow laws for crustal rocks', not accessible
    data = DislocationCreep(;
        Name = "Plagioclase An75 | Ji and Zhao (1993)",
        n = 3.2NoUnits,
        A = 3.27e-4MPa^(-16 // 5) / s,
        E = 238kJ / mol,
        V = 0m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Ji & Zhao (1993), Table 2, row 3 'plagioclase (Ab25An75)', dry: A = 3.27e-4 MPa^-n/s, n = 3.2, Q = 238 kJ/mol, dislocation creep; cited there from Shelton & Tullis (1981) (EOS abstract, not seen). Ranalli (1995), Table 10.3, gives the same law. Rybacki & Dresen (2000, Fig. 6) refit the same Shelton data as wet, as-received material: log A = -3.7, n = 3.2, Q = 239; Ji & Zhao list 3.27e-4 as 'dry'.",
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
    )
    return data, info
end

function dry_anorthite_Rybacki_2006()
    # Rybacki, Gottschalk, Wirth and Dresen (2006), table 5
    data = DislocationCreep(;
        Name = "Dry Anorthite | Rybacki et al. (2006)",
        n = 3.0NoUnits,
        A = exp10(12.7)MPa^(-3) / s,
        E = 641kJ / mol,
        V = 24.0e-6m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Rybacki et al. (2006), Table 5, 'Dislocation creep, dry': log A = 12.7 MPa^-n/s, n = 3, Q = 641, V = 24 cm^3/mol. Reproduces the paper's data: log10 misfit -0.06 ± 0.09 (N = 20, Rybacki & Dresen 2000 data, which Table 5 refits).",
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

function wet_anorthite_Rybacki_2006()
    # Rybacki, Gottschalk, Wirth and Dresen (2006), table 5
    data = DislocationCreep(;
        Name = "Wet Anorthite | Rybacki et al. (2006)",
        n = 3.0NoUnits,
        A = (10^0.2)MPa^(-4) / s,       # MPa^-(n+r), fH2O in MPa
        E = 345kJ / mol,
        V = 38.0e-6m^3 / mol,
        r = 1NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Rybacki et al. (2006), Table 5, 'Dislocation creep, wet': log A = 0.2 MPa^-(n+r)/s, n = 3, r = 1, Q = 345, V = 38 cm^3/mol. Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit -0.04 to -0.08 ± 0.09 (N = 38, Rybacki & Dresen 2000 data, which Table 5 refits; fH2O ≈ 214-251 MPa).",
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

function wet_olivine_Hirth_2003()
    #  Hirth and Kohlstedt (2003), table 1, "wet dislocation" (constant fH2O, not constant C_OH)
    data = DislocationCreep(;
        Name = "Wet Olivine | Hirth and Kohlstedt (2003)",
        n = 3.5NoUnits,
        A = 1600.0MPa^(-47 // 10) / s,  # MPa^-(n+r), fH2O in MPa
        E = 520kJ / mol,
        V = 22.0e-6m^3 / mol,
        r = 1.2NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'wet dislocation': A = 1600 MPa^-(n+r)/s (fH2O in MPa), n = 3.5, r = 1.2, E = 520, V = 22. Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit +0.09 ± 0.19 (N = 34, Mei & Kohlstedt 2000 wet data, σ ≥ 120 MPa); -0.17 to -0.25 vs the Fig. 5a line.",
        BibTex_Reference = "
            @article{hirth2003rheology,
            title={Rheology of the upper mantle and the mantle wedge: A view from the experimentalists},
            author={Hirth, Greg and Kohlstedt, D},
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
        Name = "Wet Quarzite | Kirby (1983)",
        n = 2.4NoUnits,
        A = 1.0e-4MPa^(-12 // 5) / s, # 10^3.2 GPa^-2.4/s in the source = 10^(3.2 - 7.2) MPa^-2.4/s
        E = 160.0e3J / mol,
        V = 0m^3 / mol,
        r = 0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Kirby (1983), Table 2, row 'quartzite (wet), α-quartz, water from talc' (Koch et al. 1983): log10 A = 3.2 GPa^-n/s (column unit GPa^-n s^-1, differential stress), n = 2.4, E = 160 kJ/mol; A converted to 10^(3.2 - 3n) = 1.0e-4 MPa^-n/s. Same law in Kirby & Kronenberg (1987), Table 3. Reproduces the paper's data: log10 misfit +0.08 ± 0.25 (N = 35, Koch et al. 1989 data; indirect).",
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
    )
    return data, info
end

function wet_olivine_Afonso_2004()
    #  Schmalholz, Kaus, Burg (2009), table 1
    #  referring to Afonso and Ranalli (2004), table 1, wet peridotite
    #  referring to Chopra and Paterson papers, but values dont fit the Afonso and Ranalli (2004) ones
    data = DislocationCreep(;
        Name = "Wet Upper Mantle Olivine | Afonso and Ranalli (2004)",
        n = 4.0NoUnits,
        A = 2.0e3MPa^(-4) / s,
        E = 471kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Afonso & Ranalli (2004), Table 1, row 'Peridotite (wet)': A = 2.0e3 MPa^-n s^-1, n = 4.0, E = 471 kJ/mol, used with σ1 - σ3 (their eq. 2). A compilation; it cites Chopra & Paterson (1981, 1984). Chopra & Paterson (1981), Table IV, has no n = 4.0/Q = 471 law, so the values must come from Chopra & Paterson (1984) (not seen). Same values in Gerya (2010), Table 17.1.",
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
    )
    return data, info
end

function granite_Tirel_2008()
    #  Tirel et al. (2008), table 1
    #  referring to Kirby and Kronenberg (1987), table 3
    #  different values for n and A in Kirby and Kronenberg (1987) compared with Tirel et al. (2008)
    data = DislocationCreep(;
        Name = "Granite | Tirel et al. (2008)",
        n = 3.2NoUnits,
        A = 1.25e-9MPa^(-16 // 5) / s,
        E = 123kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Tirel et al. (2008), Table 1, 'granite (anomaly)': A = 1.25e-9 MPa^-n/s, n = 3.2, Ea = 123. These do not match the cited Kirby & Kronenberg (1987) Table 3.",
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
    )
    return data, info
end

function dry_olivine_Gerya_2019()
    # after Exercise 6.1 of Numerical Geodynamics
    data = DislocationCreep(;
        Name = "Dry Olivine | Gerya (2019)",
        n = 3.5NoUnits,
        #A = 2.5e-17Pa^(-7//2)/s,
        A = uconvert(MPa^(-7 // 2) / s, 2.5e-17Pa^(-7 // 2) / s),
        E = 532.0kJ / mol,
        V = 0.0m^3 / mol,
        Apparatus = AxialCompression,   # used in book (according to matlab script)
        r = 0.0NoUnits,
    )
    info = MaterialParamsInfo(;
        Comment = "Gerya, Exercise 6.1 (A = 2.5e-17 Pa^-n/s, n = 3.5, Ea = 532 kJ/mol, axial compression); checked in the 2010 edition. Valid for olivine under upper mantle conditions.",
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
    )
    return data, info
end

function rock_salt_Li_Urai_2016()
    data = DislocationCreep(;
        Name = "Rock salt | Li & Urai (2016)",
        n = 5.0NoUnits,
        A = 7.26e-6MPa^(-5) / s,
        E = 53.92kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: the shear modulus is not given by Wawersik & Zeuch (1986); A0 = 7.26e-6 implies μ = 10.5 GPa, while the lines of W&Z Fig. 2 imply 11-12 GPa (A0 ≈ 3.8e-6 MPa^-5/s). The law fits W&Z's Salado (WIPP7) data (0.00 log10 units, N = 5) but is 10× too fast vs their Asse (BGRC1) data at 60 °C (+1.03, N = 8). Li & Urai (2016), Table 1, 'Salt B: dislocation creep' (A0 = 7.26e-6 MPa^-5/s, Q = 53920 J/mol, n = 5). Origin: Wawersik & Zeuch (1986), Table 1, row 'Asse - BGR' (the published fit of Albrecht & Hunsche 1980): γ̇ = C exp(-Q/RT) (τ/μ)^n with C = 9.27e14 s^-1, n = 5.0, Q = 53.92 kJ/mol, τ = σ1 - σ3, γ = e1 - e3 (eq. 3). μ is not given in the paper; A0 = C/μ^5 implies μ ≈ 10.5 GPa (or 9.7 GPa if γ̇ was converted to axial strain rate, 2/3 γ̇), so A0 could not be reproduced exactly.",
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
    )
    return data, info
end

function wet_olivine_Mei_2000b()
    #  Mei & Kohlstedt (2000b), eq. (2), V = 20 cm^3/mol variant
    data = DislocationCreep(;
        Name = "Wet Olivine | Mei & Kohlstedt (2000b)",
        n = 3.0NoUnits,
        A = (1.5e3)MPa^(-199 // 50) / s,     # MPa^-(n+r), fH2O in MPa
        E = 470.0kJ / mol,
        V = 20.0e-6m^3 / mol,
        r = 0.98NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: eq. (2) fits the total strain rate at σ ≥ 120 MPa, including diffusion creep (+0.20 ± 0.16 log10 units, N = 34, Table 1 wet data); it is 2.4× above the paper's own dislocation-creep component (Fig. 2), and combined with wet_olivine_Mei_2000a in a composite rheology it overpredicts the data about 2.2×. Mei & Kohlstedt (2000b), eq. (2), section 4.1: ε̇ = A σ^3.0 fH2O^r exp(-(470 kJ/mol + PV)/RT), σ (differential), P, fH2O in MPa; the V = 20 cm^3/mol fit, r = 0.98, A = 1.5e3 MPa^-(n+r)/s (printed as 'μm^3 MPa^-4 s^-1'; the μm^3 is a misprint, there is no grain-size term). The other fits are V = 0: r = 0.69, A = 4.6e3; V = 38: r = 1.25, A = 5.0e2. Triaxial compression, gas medium. Pass f = fH2O in Pa.",
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
    )
    return data, info
end

function dry_olivine_Karato_2003()
    #  Karato & Jung (2003), Table 3
    data = DislocationCreep(;
        Name = "Dry Olivine | Karato & Jung (2003)",
        n = 3.0NoUnits,
        A = (10^6.1)MPa^(-3) / s,
        E = 510.0kJ / mol,
        V = 14.0e-6m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Karato & Jung (2003), Table 3 (open system) and abstract: A_d = 10^6.1 s^-1 MPa^-n, n = 3.0, r = 0, E*_d = 510 kJ/mol (taken there from Mei & Kohlstedt 2000), V*_d = 14 cm^3/mol. σ is differential stress (dislocation-density piezometer; Griggs simple-shear runs converted assuming isotropic plasticity). Reproduces the paper's data: log10 misfit -0.07 ± 0.24 (N = 5, Table 1 dry).",
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
    )
    return data, info
end

function wet_olivine_Karato_2003()
    #  Karato & Jung (2003), Table 3
    data = DislocationCreep(;
        Name = "Wet Olivine | Karato & Jung (2003)",
        n = 3.0NoUnits,
        A = (10^2.9)MPa^(-21 // 5) / s,
        E = 470.0kJ / mol,
        V = 24.0e-6m^3 / mol,
        r = 1.2NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: the printed law is 1.6× too fast vs the paper's own data: +0.21 ± 0.07 log10 units vs all 8 wet strengths of Table 1 and +0.19 to +0.20 vs the Fig. 5b fit line, which both correspond to A_w = 10^2.7; a data-consistent version is wet_olivine_Karato_2003_refit. Karato & Jung (2003), Table 3 (open system, constant fH2O) and abstract: A_w = 10^2.9 s^-1 MPa^-(n+r) (fH2O in MPa), n = 3.0, r = 1.2, E*_w = 470 kJ/mol (from Mei & Kohlstedt 2000), V*_w = 24 cm^3/mol. Table 4 gives the closed-system (C_OH) variant: 10^0.56, E = 410, V = 11. Differential stress. Pass f = fH2O in Pa.",
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
    )
    return data, info
end

function wet_olivine_Karato_2003_refit()
    #  Karato & Jung (2003), Table 3
    #  A refitted to the paper's own data (see Comment); the printed value is in wet_olivine_Karato_2003
    data = DislocationCreep(;
        Name = "Wet Olivine (refit to data) | Karato & Jung (2003)",
        n = 3.0NoUnits,
        A = (10^2.7)MPa^(-21 // 5) / s,
        E = 470.0kJ / mol,
        V = 24.0e-6m^3 / mol,
        r = 1.2NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: refit of wet_olivine_Karato_2003 with A = 10^2.7 MPa^-(n+r)/s (not a printed value) so that the law reproduces the paper's own data. The printed A_w = 10^2.9 is +0.21 ± 0.07 log10 units (1.6× too fast) vs the 8 wet strengths of Karato & Jung (2003), Table 1 (all positive), and +0.19 to +0.20 vs their Fig. 5b fit line; 10^2.7 fits both. n = 3, r = 1.2, E = 470 and V = 24 cm^3/mol as in wet_olivine_Karato_2003, which holds the printed law (Table 3). Pass f = fH2O in Pa.",
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
    )
    return data, info
end

function wet_clinopyroxene_Chen_2006()
    #  Chen et al. (2006), eq. (3): σ and fH2O in MPa, so A is in MPa^-(n+r) s^-1
    data = DislocationCreep(;
        Name = "Wet Clinopyroxene | Chen et al. (2006)",
        n = 2.7NoUnits,
        A = (10^6.7)MPa^(-57 // 10) / s,
        E = 670.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 3.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Chen et al. (2006), eq. (3), paragraph [26]: ε̇ = 10^6.7 σ^2.7 fH2O^3.0 exp(-670 kJ/mol / RT), ε̇ in s^-1, σ (differential) and fH2O in MPa, so A = 10^6.7 MPa^-(n+r)/s; V = 0 (paragraph [25]). Triaxial compression, gas medium. Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit +0.11 ± 0.15 (N = 53, Table 2).",
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
    )
    return data, info
end

function dry_clinopyroxene_Bystricky_Mackwell_2001()
    #  Bystricky & Mackwell (2001), eq. (2): natural coarse-grained Sleaford Bay clinopyroxenite.
    #  The hot-pressed fine-grained samples give 10^10.8 (dislocation term of eq. (5), n and Q fixed)
    data = DislocationCreep(;
        Name = "Dry Clinopyroxene | Bystricky & Mackwell (2001)",
        n = 4.7NoUnits,
        A = (10^9.8)MPa^(-47 // 10) / s,
        E = 760.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Bystricky & Mackwell (2001), eq. (2), natural Sleaford Bay clinopyroxenite: ε̇ = 10^9.8 σ^4.7 exp(-760 kJ/mol / RT), ε̇ in s^-1, σ = σ1 - σ3 in MPa. Same values in Bürgmann & Dresen (2008), Suppl. Table 1. Reproduces the paper's data: log10 misfit -0.01 ± 0.09 (N = 28, Table 2, natural samples).",
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
    #  Dimanov & Dresen (2005), table 3b
    data = DislocationCreep(;
        Name = "Dry Diopside | Dimanov & Dresen (2005)",
        n = 5.47NoUnits,
        A = uconvert(MPa^(-547 // 100) / s, 3.01e-28Pa^(-547 // 100) / s),
        E = 691.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Dimanov & Dresen (2005), Table 3b, row DiD: A = 3.01e-28 Pa^-n/s, n = 5.47, Q = 691. Reproduces the paper's data: log10 misfit -0.08 ± 0.07 (N = 14, Table 1, dislocation-dominated).",
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
    #  Dimanov & Dresen (2005), table 3a
    data = DislocationCreep(;
        Name = "Wet Diopside | Dimanov & Dresen (2005)",
        n = 5.52NoUnits,
        A = uconvert(MPa^(-552 // 100) / s, 5.16e-33Pa^(-552 // 100) / s),
        E = 534.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: over-predicts the paper's own high-stress (dislocation-dominated) data by about 0.27 log10 units (N = 30), as do the paper's Fig. 5a curves. Dimanov & Dresen (2005), Table 3a, row DiW: A = 5.16e-33 Pa^-n/s, n = 5.52, Q = 534.",
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

function wet_omphacite_Zhang_2006()
    #  Zhang et al. (2006), equation (4)
    data = DislocationCreep(;
        Name = "Wet Omphacite | Zhang et al. (2006)",
        n = 3.5NoUnits,
        A = (10^-2)MPa^(-7 // 2) / s,
        E = 310.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Zhang, Green & Bozhilov (2006), eq. (4): ε̇ = 10^-2.0 σ^3.5 exp(-310000/RT), polycrystalline omphacite (Di58Jd42), Griggs apparatus at 3 GPa, 1300-1500 K. Units of σ are not printed with eq. (4); MPa follows from Fig. 4 (log stress in MPa), with which the law agrees. Same values in Bürgmann & Dresen (2008), Suppl. Table 1. Reproduces the paper's data: log10 misfit +0.22 ± 0.19 (N = 14, Fig. 4).",
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
    )
    return data, info
end

function wet_jadeit_Orzol_2006()
    #  Orzol et al. (2006), paragraph [36] and abstract: ln A = -3.3 (natural log)
    data = DislocationCreep(;
        Name = "Wet Jadeit | Orzol et al. (2006)",
        n = 3.7NoUnits,
        A = exp(-3.3)MPa^(-37 // 10) / s,
        E = 326.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Orzol et al. (2006), eq. (1) and paragraph [36] (also abstract): ε̇ = A Δσ^n exp(-Q/RT), A in MPa^-n s^-1, ln A = -3.3 ± 2.0, n = 3.7 ± 0.4, Q = 326 ± 27 kJ/mol; synthetic wet jadeitite, Griggs-type apparatus (molten salt), 2.5 GPa. A = exp(-3.3) = 3.7e-2 MPa^-n/s; the former 10^-3.3 (Bürgmann & Dresen 2008, Suppl. Table 1) read ln A as log10 A (factor 74). Reproduces the paper's data: log10 misfit +0.07 ± 0.16 (N = 10, Table 2).",
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
    )
    return data, info
end

function dry_anorthite_Rybacki_2000()
    #  Rybacki & Dresen (2000), table 2
    data = DislocationCreep(;
        Name = "Dry Anorthite | Rybacki & Dresen (2000)",
        n = 3.0NoUnits,
        A = (10^12.7)MPa^(-3) / s,
        E = 648.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Rybacki & Dresen (2000), Table 2, 'Dislocation creep, dry': log A = 12.7 MPa^-n/s, n = 3, Q = 648. Reproduces the paper's data: log10 misfit -0.06 ± 0.09 (N = 20, Table 1).",
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

function wet_anorthite_Rybacki_2000()
    #  Rybacki & Dresen (2000), table 2
    data = DislocationCreep(;
        Name = "Wet Anorthite | Rybacki & Dresen (2000)",
        n = 3.0NoUnits,
        A = (10^2.6)MPa^(-3) / s,
        E = 356.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Rybacki & Dresen (2000), Table 2, 'Dislocation creep, wet': log A = 2.6 MPa^-n/s, n = 3, Q = 356 (water content fixed, no fH2O term). Reproduces the paper's data: log10 misfit -0.02 ± 0.09 (N = 38, Table 1).",
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

function wet_quartzite_Rutter_2004()
    #  Rutter & Brodie (2004a), abstract and section 4.1
    data = DislocationCreep(;
        Name = "Wet Quartzite | Rutter & Brodie (2004)",
        n = 2.97NoUnits,
        A = (10^-4.93)MPa^(-397 // 100) / s,    # MPa^-(n+r), fH2O in MPa
        E = 242.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 1.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Rutter & Brodie (2004a), abstract and section 4.1: ε̇ = 10^-4.93 σ^2.97 fH2O exp(-242 kJ/mol / RT), σ (differential) and fH2O in MPa, so A in MPa^-(n+r)/s with r = 1 (taken from Hirth et al. 2001; fH2O assumed 300 MPa in the fit). Argon gas-medium apparatus, 300 MPa, 1273-1473 K. Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit +0.21 ± 0.31 (N = 18, Table 1); within ±0.08 of the Fig. 5 fit lines.",
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
    )
    return data, info
end

function wet_quartzite_Hirth_2001()
    #  Hirth et al. (2001), table 5
    data = DislocationCreep(;
        Name = "Wet Quartzite | Hirth et al. (2001)",
        n = 4.0NoUnits,
        A = (10^-11.2)MPa^(-5) / s,     # MPa^-(n+r), fH2O in MPa
        E = 135.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 1.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hirth et al. (2001), abstract: log A = -11.2 (fH2O in MPa), n = 4, m = 1, Q = 135. Pass f = fH2O in Pa.",
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
    )
    return data, info
end

function dry_quartzite_Jaoul_1984()
    #  Jaoul et al. (1984), table 1, first entry (W-292)
    data = DislocationCreep(;
        Name = "Dry Quartzite | Jaoul et al. (1984)",
        n = 2.8NoUnits,
        A = (1.37 * 100^-2.8)MPa^(-14 // 5) / s,   # 1.37 s^-1 kbar^-n in the source; 1 kbar = 100 MPa
        E = 184.096kJ / mol,                       # 44.0 kcal/mol
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Jaoul et al. (1984), Table 1, row 1 (W-292, creep in CaCO3 assembly, vacuum heated at 685 °C for 4 days, 700-1050 °C, 1.5 GPa): A = 1.37 s^-1 kbar^-n, n = 2.8 ± 0.2, E = 44.0 ± 1.5 kcal/mol; σ = σ1 - σ3 in kbar (Griggs apparatus). A converted to 1.37 × 100^-2.8 = 3.44e-6 MPa^-n/s (previously 10^-5.415, of unknown origin). Reproduces the paper's data: log10 misfit -0.10 ± 0.09 (N = 15, Fig. 5a).",
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
    )
    return data, info
end

function wet_quartzite_Jaoul_1984()
    #  Jaoul et al. (1984), table 1, second entry (W-291, "as is")
    data = DislocationCreep(;
        Name = "Wet Quartzite | Jaoul et al. (1984)",
        n = 2.4NoUnits,
        A = (0.63 * 100^-2.4)MPa^(-12 // 5) / s,   # 0.63 s^-1 kbar^-n in the source; 1 kbar = 100 MPa
        E = 163.176kJ / mol,                       # 39.0 kcal/mol
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Jaoul et al. (1984), Table 1, row 2 (W-291, creep in CaCO3 assembly, 'as is' = not dried, no water added, 850-1050 °C, 1.5 GPa): A = 0.63 s^-1 kbar^-n, n = 2.4 ± 0.4, E = 39.0 ± 3.0 kcal/mol; σ = σ1 - σ3 in kbar (Griggs apparatus). A converted to 0.63 × 100^-2.4 = 9.98e-6 MPa^-n/s. Previously n = 2.8, A = 10^-5.045 (E matched, n was the dry value). Water-added row W-294: A = 3.32 kbar^-n/s, n = 1.4, E = 35 kcal/mol. The data give n = 2.28 (fit lines 2.25-2.63); the Fig. 5b caption's n = 2.8 is a caption error. Reproduces the paper's data: log10 misfit +0.04 ± 0.04 (N = 14, Fig. 5b).",
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
    )
    return data, info
end

function wet_quartzite_Tokle_2019()
    #  Tokle et al. (2019), table 1, 2nd extrapolated fit
    data = DislocationCreep(;
        Name = "Wet Quartzite | Tokle et al. (2019)",
        n = 3.0NoUnits,
        A = 1.1e-12MPa^(-42 // 10) / s,    # MPa^-(n+r), fH2O in MPa
        E = 115.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 1.2NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Tokle et al. (2019), Table 1, 'Extrapolated fit', row n = 3 (low temperature/high stress, basal <a> fabrics): A = 1.1e-12, n = 3, r = 1.2, Q = 115 kJ/mol. The header gives A in MPa^-n/s, but fH2O is in MPa throughout the paper, so A is in MPa^-(n+r)/s. σ is differential (von Mises equivalent) stress. Other rows: laboratory fits n = 4 (8.0e-12, r = 1, Q = 140) and n = 2.7 (5.4e-12, r = 1.1, Q = 105). Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit -0.11 ± 0.27 (N = 51, supplement, n ≈ 3 regime).",
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
    )
    return data, info
end

function lowstress_wet_quartzite_Tokle_2019()
    #  Tokle et al. (2019), table 1, extrapolated n = 4 (high-T/low-stress) fit. The authors interpret it as
    #  dislocation creep rate-limited by easy prism <a> slip, i.e. dislocation-accommodated GBS (section 4, conclusions)
    data = DislocationCreep(;
        Name = "low stress wet Quartzite | Tokle et al. (2019)",
        n = 4.0NoUnits,
        A = 1.75e-12MPa^(-5) / s,       # 10^-11.757, MPa^-(n+r), fH2O in MPa
        E = 125.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 1.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Tokle et al. (2019), Table 1, 'Extrapolated fit', row n = 4 (high temperature/low stress, prism <a> / Y-max fabrics): A = 1.75e-12, n = 4, r = 1, Q = 125 kJ/mol; A in MPa^-(n+r)/s with fH2O in MPa (header: MPa^-n/s). Listed as a dislocation-creep law; the authors propose that it reflects disGBS rate-limited by easy prism <a> slip. Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit +0.02 ± 0.45 (N = 94, supplement, T > 900 °C).",
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
    )
    return data, info
end

function wet_quartzite_Lu_2019()
    #  Lu and Jiang (2019), section 3, equation (6)
    data = DislocationCreep(;
        Name = "Wet Quartzite | Lu and Jiang (2019)",
        n = 4.0NoUnits,
        A = (6.0e-15)MPa^(-67 // 10) / s,
        E = 132.0kJ / mol,
        V = 35.3e-6m^3 / mol,
        r = 2.7NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: the printed A is 6× too fast vs the paper's own data: +0.79 ± 0.34 log10 units vs the Gleason & Tullis (1995) and Luan & Paterson (1992) data it was fitted to (N = 65) and +0.69 vs the paper's Fig. 4; a data-consistent version is wet_quartzite_Lu_2019_refit. Lu & Jiang (2019), eq. (6): ε̇ = 6.0e-15 fw^2.7 exp(-(132000 + 35.3 P)/RT) σ^4, fw, P and σ (differential) in MPa; A = (6.0 ± 5.0)e-15 MPa^-(n+m)/s, V = 35.3 cm^3/mol. Pass f = fH2O in Pa.",
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
    )
    return data, info
end

function wet_quartzite_Lu_2019_refit()
    #  Lu and Jiang (2019), section 3, equation (6)
    #  A refitted to the paper's own data (see Comment); the printed value is in wet_quartzite_Lu_2019
    data = DislocationCreep(;
        Name = "Wet Quartzite (refit to data) | Lu and Jiang (2019)",
        n = 4.0NoUnits,
        A = (1.0e-15)MPa^(-67 // 10) / s,
        E = 132.0kJ / mol,
        V = 35.3e-6m^3 / mol,
        r = 2.7NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: refit of wet_quartzite_Lu_2019 with A = 1.0e-15 MPa^-(n+m)/s (not a printed value) so that the law reproduces the data the paper fitted. The printed A = 6.0e-15 is +0.79 ± 0.34 log10 units (6× too fast) vs the Gleason & Tullis (1995) and Luan & Paterson (1992) data (N = 65; with Lu's own treatment: G&T stresses × 0.73, P = Pc + σ/3, Pitzer-Sterner fH2O) and +0.69 vs the reduced data of Lu & Jiang's Fig. 4. The best fit is 9.7e-16; 1.0e-15 lies within the paper's stated A = (6.0 ± 5.0)e-15. n = 4, r = 2.7, E = 132 and V = 35.3 cm^3/mol as in wet_quartzite_Lu_2019, which holds the printed law (eq. 6). Pass f = fH2O in Pa.",
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
    )
    return data, info
end

function lowP_wet_quartzite_Lusk_2021()
    #  Lusk et al. (2021), abstract, 1st law
    data = DislocationCreep(;
        Name = "low pressure wet Quartzite | Lusk et al. (2021)",
        n = 3.5NoUnits,
        A = (10^-9.3)MPa^(-399 // 100) / s,
        E = 118.0kJ / mol,
        V = 2.59e-6m^3 / mol,
        r = 0.49NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Lusk et al. (2021), Table 2, row 'Low pressure' (< 560 MPa), and abstract: log A = -9.30 MPa^-(n+r)/s, n = 3.5, r = 0.49, Q = 118, V = 2.59 cm^3/mol (r and V from the 'representative reduced' fit). σ is von Mises equivalent stress = differential stress for axial tests. Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit +0.04 ± 0.27 (N = 49, Fig. 6a).",
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
    )
    return data, info
end

function wet_quartzite_Lusk_2021()
    #  Lusk et al. (2021), table 2, full data set
    data = DislocationCreep(;
        Name = "Wet Quartzite | Lusk et al. (2021)",
        n = 2.1NoUnits,
        A = (10^-6.36)MPa^(-23 // 10) / s,
        E = 94.0kJ / mol,
        V = 1.44e-6m^3 / mol,
        r = 0.2NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: matches Table 2, but is +2.16 ± 0.51 log10 units (about 150×) too fast vs the paper's own natural samples (N = 11, Fig. 1c; lab data +0.04 to +0.24); the authors recommend the lowP/highP pair (lowP_wet_quartzite_Lusk_2021, highP_wet_quartzite_Lusk_2021) instead. Lusk et al. (2021), Table 2, row 'Full': log A = -6.36 (+0.25/-0.69) MPa^-(n+r)/s, n = 2.1, r = 0.20, Q = 94, V = 1.44 cm^3/mol. σ is von Mises equivalent stress = differential stress for axial tests. Pass f = fH2O in Pa.",
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
    )
    return data, info
end


function highP_wet_quartzite_Lusk_2021()
    #  Lusk et al. (2021), abstract, 2nd law
    data = DislocationCreep(;
        Name = "high pressure wet Quartzite | Lusk et al. (2021)",
        n = 2.0NoUnits,
        A = (10^-7.9)MPa^(-249 // 100) / s,
        E = 77.0kJ / mol,
        V = 2.59e-6m^3 / mol,
        r = 0.49NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Lusk et al. (2021), Table 2, row 'High pressure' (700-1600 MPa), and abstract: log A = -7.90 MPa^-(n+r)/s, n = 2.0, r = 0.49, Q = 77, V = 2.59 cm^3/mol (r and V from the 'representative reduced' fit). σ is von Mises equivalent stress = differential stress for axial tests. Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit +0.25 ± 0.46 (N = 158, Fig. 7a).",
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
    )
    return data, info
end

function high_stress_wet_dunite_Chopra_1981()
    # Chopra & Paterson (1981) Table IV, Anita Bay dunite ($$\sigma$$ > ~ 100MPa)
    data = DislocationCreep(;
        Name = "High Stress Wet Dunite | Chopra & Paterson (1981)",
        n = 3.35NoUnits,
        A = (10^3.98)MPa^(-67 // 20) / s,
        E = 444.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Chopra & Paterson (1981), Table IV, row 'Anita Bay dunite, constant ε̇ (σ > ~100 MPa)': log A = 3.98 ± 0.17 s^-1 MPa^-n, n = 3.35, Q = 444 kJ/mol; eq. (1) with differential stress, gas-medium apparatus, 300 MPa confining pressure, 'wet' (water from hydrous minerals). Matches Carter & Tsenn (1987), Table 4, and Kirby & Kronenberg (1987), Table 1. Reproduces the paper's data: log10 misfit -0.01 ± 0.16 (N = 22, Table II).",
        BibTex_Reference = "
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
        Name = "Low Stress Wet Dunite | Chopra & Paterson (1981)",
        n = 2.44NoUnits,
        A = (10^3.88)MPa^(-61 // 25) / s,
        E = 386.0kJ / mol,
        V = 0.0m^3 / mol,
        r = 0.0NoUnits,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Chopra & Paterson (1981), Table IV, row 'Anita Bay dunite, relaxation (σ < ~100 MPa)': log A = 3.88 ± 0.43 s^-1 MPa^-n, n = 2.44, Q = 386 kJ/mol; eq. (1) with differential stress, from stress-relaxation tests, gas-medium apparatus, 300 MPa confining pressure, 'wet'. The authors note Q may be overestimated. Reproduces the paper's Fig. 6 relaxation isotherms within ±0.11 log10 units (N = 4 lines, 1000-1300 °C).",
        BibTex_Reference = "
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
