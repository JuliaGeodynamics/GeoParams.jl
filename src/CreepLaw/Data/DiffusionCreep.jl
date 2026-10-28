module Diffusion

using GeoParams
export diffusion_database, diffusion_database_info

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
        Comment = "Rybacki et al. (2006), Table 5, 'Diffusion creep, dry': log A = 12.1 MPa^-n μm^m/s, Q = 460 kJ/mol, V = 24 cm^3/mol, m = 3. Table 5 refits the Rybacki & Dresen (2000) data (primary data: Rybacki & Dresen 2000, Table 1): log10 misfit -0.04 ± 0.09 (N = 19). Vs the 2006 paper's own dry runs it is about 4× too fast (+0.64 ± 0.26, N = 113); the paper's direct fit to those is log A = 12.4, Q = 482, V = 32.4 cm^3/mol.",
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
        Comment = "Mei & Kohlstedt (2000a), eq. (4), section 4.1: ε̇ = A σ^1.1 d^-3 fH2O^r exp(-(295 kJ/mol + PV)/RT), σ (differential), P, fH2O in MPa, d in μm; the V = 20 cm^3/mol fit: r = 1.0, A = 4.8e4 μm^3 MPa^-2.1 s^-1 (V = 0: r = 0.7, A = 1.3e5 μm^3 MPa^-1.8 s^-1). Triaxial compression, gas medium. Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit +0.11 ± 0.20 (N = 28, Table 1 of Mei & Kohlstedt 2000b, σ ≤ 60 MPa).",
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
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'wet diffusion' (constant fH2O): A = 2.5e7 for r = 1 (footnote d), E = 375, V = 0-20 (upper end used). A includes the fH2O unit (MPa^-(n+r)); pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit -0.06 ± 0.19 (N = 28; primary data: Mei & Kohlstedt 2000b, Table 1, wet runs at σ ≤ 60 MPa).",
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
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'wet diffusion (constant C_OH)': A = 1.0e6, r = 1, E = 335, V = 4 (footnote c example: 7.8e-15 /s). Pass f = C_OH in H/10^6 Si. Reproduces the paper's data: log10 misfit +0.02 ± 0.18 (N = 28; primary data: Mei & Kohlstedt 2000b, Table 1, wet runs at σ ≤ 60 MPa).",
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
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'dry diffusion': A = 1.5e9 MPa^-1 μm^3/s, E = 375, V = 2-10 (upper end used). Reproduces the paper's data: log10 misfit +0.16 ± 0.21 (N = 8; primary data: Mei & Kohlstedt 2000b, Table 1, dry runs at σ ≤ 60 MPa).",
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
        Comment = "Faul & Jackson (2007), eq. (2) and paragraph [32]: A = 2e10 MPa^-n μm^p/s, E = 484 kJ/mol, n = 1.37, p = 3; dry, melt-free sol-gel olivine, V not determined. Formerly wet_olivine_Faul_2006. The 'mm' printed in the PDF for A and the grain size is a misprint for μm. Reproduces the paper's data: log10 misfit +0.02 ± 0.09 (N = 77, Table 1, melt-free runs).",
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
        Comment = "Hier-Majumder et al. (2005), eq. (3) and Table 3, row 'Anhydrous': ε̇ = ε̇o (σ/G)^n (d/do)^p exp(-Q/RT) with ε̇o = 10^25.3 s^-1, n = 1, p = -3, Q = 760 kJ/mol, G = 64.9 GPa, do = 10 μm, σ differential stress. In MPa and μm: A = ε̇o do^3 / G = 10^25.3 × 10^3 / 64.9e3 = 3.07e23 MPa^-1 μm^3/s (equal to Bürgmann & Dresen 2008, 10^25.3/64.9; the former 10^25.3/64.9e3 lacked the do^3 factor). Reproduces the paper's data: log10 misfit +0.03 ± 0.10 (N = 36, Table 2a).",
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
        Comment = "Hier-Majumder et al. (2005), eq. (4) and Table 3, row 'Water-saturated': ε̇ = ε̇o (σ/G)^n (d/do)^p (fH2O/fH2O°)^r exp(-(Q + PV)/RT) with ε̇o = 10^7.9 s^-1, n = 1, p = -3, r = 1.4, Q = 340 kJ/mol, V = 14 cm^3/mol, G = 64.9 GPa, do = 10 μm, fH2O° = 1 MPa, σ differential stress. In MPa and μm: A = 10^7.9 × 10^3 / 64.9e3 = 1.22e6 MPa^-(n+r) μm^3/s (equal to Bürgmann & Dresen 2008; the former value lacked the do^3 factor). Pass f = fH2O in Pa. Reproduces the paper's data: log10 misfit -0.10 ± 0.21 (N = 27, Table 2a).",
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
        Comment = "Bystricky & Mackwell (2001), eq. (5), diffusion term (dry hot-pressed clinopyroxenite, Ni/NiO buffer): ε̇ = 10^15.1 σ d^-3 exp(-560 kJ/mol / RT), ε̇ in s^-1, σ = σ1 - σ3 in MPa; d in μm (unit not printed, implied by the 8 μm normalisation and the strain rates). Same values in Bürgmann & Dresen (2008), Suppl. Table 1. Reproduces the paper's data: log10 misfit 0.00 ± 0.15 (N = 20, Fig. 5 curves, d = 8 μm).",
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
        Comment = "Dimanov & Dresen (2005), Table 3b, row Di35D (n = 1, Q = 528); A = 10^14 from Bürgmann & Dresen (2008), Suppl. Table 1, since the paper folds d^m into A (m = 3 assumed). A implies an effective grain size of about 13-15 μm, which the paper does not state (it gives only particle size classes).",
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
        Comment = "Bürgmann & Dresen (2008), Suppl. Table 1 (log A = 8.1, Q = 337, m = 3); Q is the mean of Dimanov & Dresen (2005), Table 3a, rows Di25W/Di35W/Di45W (345, 349, 318 kJ/mol). The paper folds d^m into A. A implies an effective grain size of about 13-15 μm, which the paper does not state (it gives only particle size classes). Primary data: Dimanov & Dresen (2005), not refitted here.",
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
        Comment = "Rybacki & Dresen (2000), Table 2, 'Diffusion creep, dry': log A = 12.1 MPa^-n μm^m/s, Q = 467 kJ/mol, m = 3. Reproduces the paper's data: log10 misfit -0.03 ± 0.09 (N = 19, Table 1).",
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
        Comment = "Rybacki & Dresen (2000), Table 2, 'Diffusion creep, wet': log A = 1.7 MPa^-n μm^m/s, Q = 170 kJ/mol, m = 3. Reproduces the paper's data: log10 misfit -0.04 ± 0.15 (N = 47, Table 1).",
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
        Comment = "Rybacki et al. (2006), Table 5, 'Diffusion creep, wet': log A = -0.7 MPa^-(n+r) μm^m/s, Q = 159 kJ/mol, V = 38 cm^3/mol, r = 1, m = 3. Pass f = fH2O in Pa. Table 5 refits the Rybacki & Dresen (2000) data (primary data: Rybacki & Dresen 2000, Table 1) (implied fH2O ≈ 250 MPa): log10 misfit -0.06 to -0.11 ± 0.15 (N = 47). Vs the 2006 paper's own wet runs it is about 0.4-0.6 log10 units too slow.",
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
        Comment = "NOTE: the printed law is about 2.5× too slow vs the paper's own data: -0.44 log10 units vs the Fig. 7 grain-size centroids (160 MPa, 1473 K) and median -0.36 vs Table 1 at σ ≤ 200 MPa (N = 10), which imply log A ≈ 0.0; this is within the stated log A = -0.4 ± 2.1. Rutter & Brodie (2004b), section 4.1 and abstract: ε̇ = A σ^n exp(-H/RT) / d^m with log A = -0.4 ± 2.1, n = 1 and m = 2 (fixed), H = 220 kJ/mol; σ (differential) in MPa, d in μm, no fH2O term (water fugacity about 20-100 MPa, not in the law). Gas-medium apparatus, 300 MPa confining pressure.",
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

function olivine_disGBS_Ohuchi_2015()
    # disGBS law with a grain-size and a water-fugacity term; DiffusionCreep is the existing form with both.
    data = DiffusionCreep(;
        Name = "Olivine disGBS | Ohuchi et al. (2015)",
        n = 3.0NoUnits,                                     # power-law exponent
        r = 1.25NoUnits,                                    # exponent of water-fugacity
        p = -1.0NoUnits,                                    # grain size exponent
        A = 10^(-4.89) * MPa^(-17 // 4) * m * s^(-1),       # MPa^-(n+r), grain size in m
        E = 423.0kJ / mol,                                  # activation energy
        V = 17.6e-6m^3 / mol,                               # activation volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Ohuchi et al. (2015), Eq. 2 (p. 5): A = 10^-4.89 MPa^-4.25 m s^-1 (grain size in m), n = 3.0, p = -1 (fixed; best fit 1.1 ± 0.3), r = 1.25 (A includes the fH2O unit; pass f = fH2O in Pa), E = 423 kJ/mol, V = 17.6 cm³/mol; Fo90 olivine aggregates, D-DIA, 1.5-6.4 GPa, 1323-1473 K, axial differential stress. Reproduces the paper's data: log10 misfit -0.02 ± 0.31 (N = 6, digitized from Fig. 3D; tabulated data are in Supplementary Tables S1-S2).",
        BibTex_Reference = "
            @article{ohuchi2015dislocation,
            title={Dislocation-accommodated grain boundary sliding as the major deformation mechanism of olivine in the Earth's upper mantle},
            author={Ohuchi, Tomohiro and Kawazoe, Takaaki and Higo, Yuji and Funakoshi, Ken-ichi and Suzuki, Akio and Kikegawa, Takumi and Irifune, Tetsuo},
            journal={Science Advances},
            volume={1},
            number={9},
            pages={e1500360},
            year={2015},
            doi={10.1126/sciadv.1500360}
            }
    ",
    )
    return data, info
end

function wet_enstatite_Zhang_2017()
    data = DiffusionCreep(;
        Name = "Wet Enstatite | Zhang et al. (2017)",
        n = 1.0NoUnits,                                     # power-law exponent
        r = 0.7NoUnits,                                     # exponent of water-fugacity
        p = -3.0NoUnits,                                    # grain size exponent
        A = 32.0MPa^(-17 // 10) * μm^3 * s^(-1),            # MPa^-(n+r), fH2O in MPa
        E = 200.0kJ / mol,                                  # activation energy
        V = 14.0e-6m^3 / mol,                               # activation volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Zhang et al. (2017), Eq. 2 (p. 7722): A = 32 MPa^-1.7 μm³ s^-1 (A includes the fH2O unit; pass f = fH2O in Pa), n = 1 and m = 3 (fixed), r = 0.7, E = 200 kJ/mol, V = 14 cm³/mol; enstatite aggregates, D-DIA, 3.8-6.3 GPa, 1323-1573 K, differential stress. Reproduces the paper's data: log10 misfit +0.10 ± 0.01 (N = 6, Table 1, fH2O converted from GPa) with d = 6 μm; per-run grain sizes are not given (~6 μm), and d = 6.5 μm removes the offset.",
        BibTex_Reference = "
            @article{zhang2017diffusion,
            title={Diffusion Creep of Enstatite at High Pressures Under Hydrous Conditions},
            author={Zhang, Guinan and Mei, Shenghua and Song, Maoshuang and Kohlstedt, David L.},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={122},
            number={10},
            pages={7718--7728},
            year={2017},
            doi={10.1002/2017JB014400}
            }
    ",
    )
    return data, info
end

function wet_quartz_Fukuda_2018()
    # Mixed diffusion/dislocation law with grain-size and water-fugacity terms.
    data = DiffusionCreep(;
        Name = "Wet fine-grained Quartz | Fukuda et al. (2018)",
        n = 1.7NoUnits,                                     # power-law exponent
        r = 1.0NoUnits,                                     # exponent of water-fugacity
        p = -0.51NoUnits,                                   # grain size exponent
        A = 10^(-2.97) * MPa^(-27 // 10) * μm^(51 // 100) * s^(-1), # MPa^-(n+r), fH2O in MPa
        E = 183.0kJ / mol,                                  # activation energy
        V = 0.0m^3 / mol,                                   # activation volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Fukuda et al. (2018), Eq. 8 (p. 4691): A = 10^-2.97 MPa^-2.7 μm^0.51 s^-1 (A includes the fH2O unit; pass f = fH2O in Pa), n = 1.7, p = -0.51, r = 1.0, E = 183 kJ/mol; wet fine-grained quartz, 800-950 °C, 1.5 GPa, Griggs apparatus, differential stress. Reproduces the paper's data: log10 misfit 0.17 ± 0.29 (N = 45, Table 3, grain sizes from Table 2). The low-T fit (Q = 129 kJ/mol) has no published A.",
        BibTex_Reference = "
            @article{fukuda2018deformation,
            title={Deformation of Fine-Grained Quartz Aggregates by Mixed Diffusion and Dislocation Creep},
            author={Fukuda, Jun-ichi and Holyoke, Caleb W. and Kronenberg, Andreas K.},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={123},
            number={6},
            pages={4676--4696},
            year={2018},
            doi={10.1029/2017JB015133}
            }
    ",
    )
    return data, info
end

function wet_plagioclase_An60_diffusion_Rybacki_2004()
    data = DiffusionCreep(;
        Name = "Wet Plagioclase An60 | Rybacki & Dresen (2004)",
        n = 1.0NoUnits,                                     # power-law exponent
        r = 0.0NoUnits,                                     # exponent of water-fugacity
        p = -3.0NoUnits,                                    # grain size exponent
        A = 10^1.1 * MPa^(-1) * μm^3 * s^(-1),              # material specific rheological parameter
        E = 153.0kJ / mol,                                  # activation energy
        V = 0.0m^3 / mol,                                   # activation volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Rybacki & Dresen (2004), Table 1 (p. 175): A = 10^(1.1 ± 0.6) MPa^-1 μm^3 s^-1, n = 1, p = -3, Q = 153 kJ/mol; synthetic An60 aggregates with 0.3 wt% H2O, data of Dimanov et al. (unpublished); differential stress. No data table; reproduces the diffusion-creep lines of Fig. 2 (p. 177) within about 0.1 dex.",
        BibTex_Reference = "
            @article{rybacki2004deformation,
            title={Deformation mechanism maps for feldspar rocks},
            author={Rybacki, Erik and Dresen, Georg},
            journal={Tectonophysics},
            volume={382},
            number={3-4},
            pages={173--187},
            year={2004},
            doi={10.1016/j.tecto.2004.01.006}
            }
    ",
    )
    return data, info
end

function blueschist_Tokle_2023()
    # Grain size d is the boudin size (microboudinage-limited diffusion creep).
    data = DiffusionCreep(;
        Name = "Blueschist | Tokle et al. (2023)",
        n = 2.0NoUnits,                                     # power-law exponent
        r = 0.0NoUnits,                                     # exponent of water-fugacity
        p = -1.0NoUnits,                                    # boudin size exponent
        # The paper fits σ = 2τ against ε̇ = γ̇/√3; SimpleShear uses ε̇ = γ̇, hence the factor √3.
        A = sqrt(3) * 2.43e11MPa^(-2) * μm * s^(-1),
        E = 384.0kJ / mol,                                  # activation energy
        V = 0.0m^3 / mol,                                   # activation volume
        Apparatus = SimpleShear,
    )
    info = MaterialParamsInfo(;
        Comment = "Tokle et al. (2023), Eq. 5 (p. 10) and p. 12: A = 2.43e11 MPa^-2 μm s^-1 for σ_eqv = 2τ and ε̇_eqv = γ̇/√3, stored as SimpleShear with A = √3 × 2.43e11; n = 2, p = -1 with d = arithmetic-mean microboudin neck length, Q = 384 kJ/mol. Na-amphibole blueschist, ~1 GPa, 650-700 °C, Griggs general shear. Reproduces the paper's data: log10 misfit 0.02 ± 0.07 (N = 20, Table 1), and the Fig. 3b lines within 0.03 dex.",
        BibTex_Reference = "
            @article{tokle2023diffusion,
            title={Diffusion Creep of Sodic Amphibole-Bearing Blueschist Limited by Microboudinage},
            author={Tokle, Leif and Hufford, Lonnie J. and Behr, Whitney M. and Morales, Luiz F. G. and Madonna, Claudio},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={128},
            number={9},
            year={2023},
            doi={10.1029/2023JB026848}
            }
    ",
    )
    return data, info
end

function calcite_diffusion_Gilgannon_2025()
    data = DiffusionCreep(;
        Name = "Calcite (fine-grained, diffusion, X_MgCO3 = 0.001) | Gilgannon & Herwegh (2025)",
        n = 1.0NoUnits,
        r = 0.0NoUnits,
        p = -3.0NoUnits,
        A = 3.40e7 * 0.001^0.8 * MPa^(-1) * μm^3 * s^(-1),   # X_MgCO3^u folded into A
        E = 200.0kJ / mol,
        V = 0.0m^3 / mol,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Gilgannon & Herwegh (2025, preprint), Eq. 5 and Table 5 (pp. 17-18): A = 3.40e7 s^-1 × X_MgCO3^0.8 with σ_ref = 1 MPa, d_ref = 1 μm; X_MgCO3 = 0.001 (the paper's pure-calcite limit, p. 21) is folded into A. This extrapolates below the fine-grained calibration range X_MgCO3 = 0.07-0.17. Q = 200 kJ/mol is assumed, not fitted. Stress is equivalent stress (torsion data converted with σ = √3 τ). Reproduces the Fig. 6a lines at σ = 1 MPa within 0.01 dex; with calcite_fine_dislocation_Gilgannon_2025 at each sample's X_MgCO3, log10 misfit -0.20 ± 0.17 (N = 72, digitized from Fig. 6a).",
        BibTex_Reference = "
            @article{gilgannon2025unifying,
            title={On unifying carbonate rheology},
            author={Gilgannon, James and Herwegh, Marco},
            journal={EGUsphere (preprint)},
            year={2025},
            doi={10.5194/egusphere-2025-1718}
            }
    ",
    )
    return data, info
end

function calcite_fine_dislocation_Gilgannon_2025()
    # Dislocation creep with a grain-size term; DiffusionCreep is the existing form with one.
    data = DiffusionCreep(;
        Name = "Calcite (fine-grained, dislocation, X_MgCO3 = 0.001) | Gilgannon & Herwegh (2025)",
        n = 7.0NoUnits,
        r = 0.0NoUnits,
        p = -3.0NoUnits,
        A = 4.00e-8 * 0.001^(-2.7) * MPa^(-7) * μm^3 * s^(-1), # X_MgCO3^u folded into A
        E = 200.0kJ / mol,
        V = 0.0m^3 / mol,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Gilgannon & Herwegh (2025, preprint), Eq. 6 and Table 5 (pp. 17-18): A = 4.00e-8 s^-1 × X_MgCO3^-2.7 with σ_ref = 1 MPa, d_ref = 1 μm, n = 7, d^-3; X_MgCO3 = 0.001 (the paper's pure-calcite limit, p. 21) is folded into A. This extrapolates below the calibration range X_MgCO3 = 0.07-0.17. Q = 200 kJ/mol is assumed. Stress is equivalent stress. With calcite_diffusion_Gilgannon_2025 at each sample's X_MgCO3, reproduces the paper's data: log10 misfit -0.20 ± 0.17 (N = 72, digitized from Fig. 6a). The n = 7 guide line of Fig. 4a implies A ≈ 9.1e-8 instead (0.36 dex faster); the Table 5 value is kept.",
        BibTex_Reference = "
            @article{gilgannon2025unifying,
            title={On unifying carbonate rheology},
            author={Gilgannon, James and Herwegh, Marco},
            journal={EGUsphere (preprint)},
            year={2025},
            doi={10.5194/egusphere-2025-1718}
            }
    ",
    )
    return data, info
end

@inline diffusion_database(f::F) where {F} = first(f())
@inline diffusion_database_info(f::F) where {F} = last(f())

end
