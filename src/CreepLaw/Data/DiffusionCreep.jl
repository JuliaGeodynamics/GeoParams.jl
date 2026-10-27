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
        Comment = "Rybacki et al. (2006), Table 5, 'Diffusion creep, dry': log A = 12.1 MPa^-n μm^m/s, Q = 460 kJ/mol, V = 24 cm^3/mol, m = 3.",
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
        n = 1.1NoUnits,                         # power-law exponent
        r = 1.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = 4.8e4MPa^(-21 // 10) * μm^3 * s^(-1), # material specific rheological parameter; MPa^-(n+r), fH2O in MPa
        E = 295.0kJ / mol,                      # activation energy
        V = 20.0e-6m^3 / mol,                     # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Mei & Kohlstedt (2000a), eq. (4), section 4.1: ε̇ = A σ^1.1 d^-3 fH2O^r exp(-(295 kJ/mol + PV)/RT), σ (differential), P, fH2O in MPa, d in μm; the V = 20 cm^3/mol fit: r = 1.0, A = 4.8e4 μm^3 MPa^-2.1 s^-1 (V = 0: r = 0.7, A = 1.3e5 μm^3 MPa^-1.8 s^-1). Triaxial compression, gas medium. Pass f = fH2O in Pa.",
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
        A = 2.5e7MPa^(-2) * μm^3 * s^(-1),      # material specific rheological parameter; MPa^-(n+r), fH2O in MPa
        E = 375.0kJ / mol,                        # activation energy
        V = 20.0e-6m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'wet diffusion' (constant fH2O): A = 2.5e7 for r = 1 (footnote d), E = 375, V = 0-20 (upper end used). A includes the fH2O unit (MPa^-(n+r)); pass f = fH2O in Pa.",
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
"""
Diffusion creep data for wet olivine at constant water content after Hirth & Kohlstedt (2003), Table 1,
row 'wet diffusion (constant C_OH)'. Here `f` is the water content C_OH in H/10^6 Si (dimensionless), not a fugacity.
"""
function wet_olivine1_Hirth_2003()
    data = DiffusionCreep(;
        Name = "Wet Olivine 1 | Hirth & Kohlstedt (2003)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 1.0NoUnits,                         # exponent of water content C_OH
        p = -3.0NoUnits,                        # grain size exponent
        A = 1.0e6MPa^(-1) * μm^3 * s^(-1),      # material specific rheological parameter; (H/10^6 Si)^-r not included
        E = 335.0kJ / mol,                        # activation energy
        V = 4.0e-6m^3 / mol,                        # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'wet diffusion (constant C_OH)': A = 1.0e6, r = 1, E = 335, V = 4 (footnote c example: 7.8e-15 /s). Pass f = C_OH in H/10^6 Si.",
        BibTex_Reference = "
            @incollection{HirthKohlstedt_OlivineRheology_2003,
            title = {Rheology of the upper mantle and the mantle wedge: {A} view from the experimentalists},
            booktitle = {Inside the Subduction Factory, Geophysical Monograph Series},
            volume = {138},
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
function dry_olivine_Hirth_2003()
    data = DiffusionCreep(;
        Name = "Dry Olivine | Hirth & Kohlstedt (2003)",
        n = 1.0NoUnits,                         # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = 1.5e9MPa^(-1) * μm^3 * s^(-1),      # material specific rheological parameter
        E = 375.0kJ / mol,                        # activation energy
        V = 10.0e-6m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'dry diffusion': A = 1.5e9 MPa^-1 μm^3/s, E = 375, V = 2-10 (upper end used).",
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
function dry_olivine_Faul_2007()
    data = DiffusionCreep(;
        Name = "Dry Olivine | Faul & Jackson (2007)",
        n = 1.37NoUnits,                        # power-law exponent
        r = 0.0NoUnits,                         # exponent of water-fugacity
        p = -3.0NoUnits,                        # grain size exponent
        A = 2.0e10MPa^(-137 // 100) * μm^3 * s^(-1),    # material specific rheological parameter
        E = 484.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Faul & Jackson (2007), eq. (2) and paragraph [32]: A = 2e10 MPa^-n μm^p/s, E = 484 kJ/mol, n = 1.37, p = 3; dry, melt-free sol-gel olivine, V not determined. Formerly wet_olivine_Faul_2006.",
        BibTex_Reference = "
            @article{faul2007diffusion,
            title={Diffusion creep of dry, melt-free olivine},
            author={Faul, Ulrich H and Jackson, Ian},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={112},
            number={B4},
            pages={B04204},
            year={2007},
            doi={10.1029/2006JB004586}
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
        A = (10^25.3 * 10^3 / 64.9e3)MPa^(-1) * μm^3 * s^(-1),    # ε̇o d_o^3 / G with d_o = 10 μm, G = 64.9e3 MPa
        E = 760.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hier-Majumder et al. (2005), eq. (3) and Table 3, row 'Anhydrous': ε̇ = ε̇o (σ/G)^n (d/do)^p exp(-Q/RT) with ε̇o = 10^25.3 s^-1, n = 1, p = -3, Q = 760 kJ/mol, G = 64.9 GPa, do = 10 μm, σ differential stress. In MPa and μm: A = ε̇o do^3 / G = 10^25.3 × 10^3 / 64.9e3 = 3.07e23 MPa^-1 μm^3/s (equal to Bürgmann & Dresen 2008, 10^25.3/64.9; the former 10^25.3/64.9e3 lacked the do^3 factor).",
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
        A = (10^7.9 * 10^3 / 64.9e3)MPa^(-24 // 10) * μm^3 * s^(-1),    # ε̇o d_o^3 / (G fH2O°^r), d_o = 10 μm, G = 64.9e3 MPa, fH2O° = 1 MPa
        E = 340.0kJ / mol,                        # activation energy
        V = 14.0e-6m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hier-Majumder et al. (2005), eq. (4) and Table 3, row 'Water-saturated': ε̇ = ε̇o (σ/G)^n (d/do)^p (fH2O/fH2O°)^r exp(-(Q + PV)/RT) with ε̇o = 10^7.9 s^-1, n = 1, p = -3, r = 1.4, Q = 340 kJ/mol, V = 14 cm^3/mol, G = 64.9 GPa, do = 10 μm, fH2O° = 1 MPa, σ differential stress. In MPa and μm: A = 10^7.9 × 10^3 / 64.9e3 = 1.22e6 MPa^-(n+r) μm^3/s (equal to Bürgmann & Dresen 2008; the former value lacked the do^3 factor). Pass f = fH2O in Pa.",
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
        Comment = "Bystricky & Mackwell (2001), eq. (5), diffusion term (dry hot-pressed clinopyroxenite, Ni/NiO buffer): ε̇ = 10^15.1 σ d^-3 exp(-560 kJ/mol / RT), ε̇ in s^-1, σ = σ1 - σ3 in MPa; d in μm (unit not printed, implied by the 8 μm normalisation and the strain rates). Same values in Bürgmann & Dresen (2008), Suppl. Table 1.",
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
        Comment = "Dimanov & Dresen (2005), Table 3b, row Di35D (n = 1, Q = 528); A = 10^14 from Bürgmann & Dresen (2008), Suppl. Table 1, since the paper folds d^m into A (m = 3 assumed).",
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
        E = 337.0kJ / mol,                        # activation energy
        V = 0.0m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Bürgmann & Dresen (2008), Suppl. Table 1 (log A = 8.1, Q = 337, m = 3); Q is the mean of Dimanov & Dresen (2005), Table 3a, rows Di25W/Di35W/Di45W (345, 349, 318 kJ/mol). The paper folds d^m into A.",
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
        Comment = "Rybacki & Dresen (2000), Table 2, 'Diffusion creep, dry': log A = 12.1 MPa^-n μm^m/s, Q = 467 kJ/mol, m = 3.",
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
        Comment = "Rybacki & Dresen (2000), Table 2, 'Diffusion creep, wet': log A = 1.7 MPa^-n μm^m/s, Q = 170 kJ/mol, m = 3.",
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
        A = (10^-0.7)MPa^(-2) * μm^3 * s^(-1),    # material specific rheological parameter; MPa^-(n+r), fH2O in MPa
        E = 159.0kJ / mol,                        # activation energy
        V = 38.0e-6m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Rybacki et al. (2006), Table 5, 'Diffusion creep, wet': log A = -0.7 MPa^-(n+r) μm^m/s, Q = 159 kJ/mol, V = 38 cm^3/mol, r = 1, m = 3. Pass f = fH2O in Pa.",
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
        Comment = "Rutter & Brodie (2004b), section 4.1 and abstract: ε̇ = A σ^n exp(-H/RT) / d^m with log A = -0.4 ± 2.1, n = 1 and m = 2 (fixed), H = 220 kJ/mol; σ (differential) in MPa, d in μm, no fH2O term (water fugacity about 20-100 MPa, not in the law). Gas-medium apparatus, 300 MPa confining pressure.",
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
