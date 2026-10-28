module GBS

using GeoParams
export GrainBoundarySliding_database, GrainBoundarySliding_database_info

function cold_dry_olivine_Hirth_2003()
    data = GrainBoundarySliding(;
        Name = "Dry Olivine < 1523K | Hirth and Kohlstedt (2003)",
        n = 3.5NoUnits,                           # power-law exponent
        p = -2.0NoUnits,                          # grain size exponent
        A = 6500.0MPa^(-7 // 2) * μm^(2) * s^(-1),  # material specific rheological parameter
        E = 400.0kJ / mol,                        # activation energy
        V = 18.0e-6m^3 / mol,                     # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'dry GBS, T < 1250 °C' (A in MPa^-3.5 μm^2/s, differential stress). V as for dry dislocation creep (footnote g; 18e-6 is the thermocouple-corrected value of Table 2). HK03 Fig. 7a draws the GBS line with p ≈ 1.4 while Table 1 gives p = 2; the law reproduces the line near d ≈ 15 μm (log10 misfit +0.09, -0.01, -0.34 at d = 10, 15, 60 μm). Primary data: Hirth & Kohlstedt (1995b), Table 1, melt-free samples hot-pressed for 4 h (PI-40, PI-66; 1250 °C, 300 MPa, steps in the dislocation regime): log10 misfit +0.04 ± 0.08 (N = 5) at d = 15 μm (+0.40 at 10 μm, -0.21 at 20 μm; grain size not tabulated, about 12-16 μm from their Fig. 3). All runs are at one temperature, so E is not constrained by them; E = 400 kJ/mol matches the (010)[100] ([110]c, orthopyroxene-buffered) mechanism-1 law of Bai et al. (1991), Table 4, Q = 230 kJ/mol plus 0.36 × the fO2 term along Ni-NiO (≈ 402; HK03 footnote f). The coarser 21 h samples (PI-72, PI-82) are not GBS-dominated; dry_olivine_Hirth_2003 (dislocation) fits them (-0.09 ± 0.08, N = 10).",
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

function hot_dry_olivine_Hirth_2003()
    data = GrainBoundarySliding(;
        Name = "Dry Olivine >= 1523K | Hirth and Kohlstedt (2003)",
        n = 3.5NoUnits,                         # power-law exponent
        p = -2.0NoUnits,                        # grain size exponent
        A = 4.7e10MPa^(-7 // 2) * μm^(2) * s^(-1), # material specific rheological parameter
        E = 600.0kJ / mol,                      # activation energy
        V = 18.0e-6m^3 / mol,                     # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'dry GBS, T > 1250 °C' (A in MPa^-3.5 μm^2/s, differential stress). V as for dry dislocation creep (footnote g; 18e-6 is the thermocouple-corrected value of Table 2). No GBS data exist above 1250 °C, so the law cannot be checked against data. HK03, Table 1, footnote f: E is assumed equal to that of slip on (010)[100] of Bai et al. (1991) single crystals (dislocation creep, 0.1 MPa), including the effect of temperature on fO2; E = 600 kJ/mol is consistent with the apparent activation energy of the Bai et al. (1991), Table 4, [110]c orthopyroxene-buffered law along Ni-NiO between 1250 and 1500 °C (599 kJ/mol on average; locally 390-770). A follows from continuity with cold_dry_olivine_Hirth_2003 at 1250 °C (ratio 1.00), where both reproduce Hirth & Kohlstedt (1995b), Table 1 (+0.05 ± 0.08, N = 5, d = 15 μm).",
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

function test_GBS()
    data = GrainBoundarySliding(;
        Name = "Test GBS",
        n = 3.5NoUnits,                         # power-law exponent
        p = -2.0NoUnits,                        # grain size exponent
        A = 1.506190693026593e2MPa^(-7 // 2) * m^(2) * s^(-1),    # material specific rheological parameter
        E = 600.0kJ / mol,                        # activation energy
        V = 18.0e-6m^3 / mol,                       # activation Volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values; not yet plots (NM)", BibTex_Reference = ""
    )
    return data, info
end

function dry_olivine_disGBS_Hansen_2011()
    data = GrainBoundarySliding(;
        Name = "Dry Olivine disGBS | Hansen et al. (2011)",
        n = 2.9NoUnits,                                     # power-law exponent
        p = -0.7NoUnits,                                    # grain size exponent
        A = 10^4.8 * MPa^(-29 // 10) * μm^(7 // 10) * s^(-1), # material specific rheological parameter
        E = 445.0kJ / mol,                                  # activation energy
        V = 0.0m^3 / mol,                                   # activation volume (not given)
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Hansen et al. (2011), Eq. 6 (p. 9) and Table A2: A = 10^4.8 MPa^-2.9 μm^0.7 s^-1, n = 2.9, p = -0.7 (d = linear intercept × 1.5), E = 445 kJ/mol; dry San Carlos olivine aggregates, gas-medium apparatus, 300 MPa, 1373-1523 K, differential stress. Reproduces the paper's data: log10 misfit 0.01 ± 0.19 (N = 71, Table 1), and the Fig. 5 lines within 0.1 dex.",
        BibTex_Reference = "
            @article{hansen2011grain,
            title={Grain boundary sliding in San Carlos olivine: Flow law parameters and crystallographic-preferred orientation},
            author={Hansen, L. N. and Zimmerman, M. E. and Kohlstedt, D. L.},
            journal={Journal of Geophysical Research},
            volume={116},
            number={B8},
            year={2011},
            doi={10.1029/2011JB008220}
            }
    ",
    )
    return data, info
end

function ice_GBS_warm_Goldsby_2001()
    data = GrainBoundarySliding(;
        Name = "Ice Ih GBS, T > 255 K | Goldsby & Kohlstedt (2001)",
        n = 1.8NoUnits,
        p = -1.4NoUnits,
        A = 3.0e26 * MPa^(-9 // 5) * m^(7 // 5) * s^(-1),
        E = 192.0kJ / mol,
        V = -13.0e-6m^3 / mol,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "NOTE: the printed law is about 2× too slow vs the paper's own data: -0.28 ± 0.12 log10 units (N = 3, Fig. 4, d = 73 μm), about 0.25 below the drawn 190 kJ/mol line; A ≈ 5.6e26 would fit, but three points do not justify a refit. At 255 K it lies 0.41 log10 units off ice_GBS_cold_Goldsby_2001. Goldsby & Kohlstedt (2001), Table 5 (p. 11,024): A = 3.0e26 MPa^-1.8 m^1.4 s^-1 (d in m), n = 1.8, p = 1.4, Q = 192 kJ/mol, V = -13 cm³/mol; ice Ih, T > 255 K, differential stress. In the paper, GBS acts in series with basal slip, which GeoParams cannot express.",
        BibTex_Reference = "
            @article{goldsby2001superplastic,
            title={Superplastic deformation of ice: Experimental observations},
            author={Goldsby, D. L. and Kohlstedt, D. L.},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={106},
            number={B6},
            pages={11017--11030},
            year={2001},
            doi={10.1029/2000JB900336}
            }
    ",
    )
    return data, info
end

function ice_GBS_cold_Goldsby_2001()
    data = GrainBoundarySliding(;
        Name = "Ice Ih GBS, T < 255 K | Goldsby & Kohlstedt (2001)",
        n = 1.8NoUnits,
        p = -1.4NoUnits,
        A = 3.9e-3 * MPa^(-9 // 5) * m^(7 // 5) * s^(-1),
        E = 49.0kJ / mol,
        V = -13.0e-6m^3 / mol,
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Goldsby & Kohlstedt (2001), Table 5 (p. 11,024): A = 3.9e-3 MPa^-1.8 m^1.4 s^-1 (d in m), n = 1.8, p = 1.4, Q = 49 kJ/mol, V = -13 cm³/mol; ice Ih, T < 255 K, differential stress. Reproduces the paper's data: log10 misfit 0.07 ± 0.09 (N = 10, Fig. 4) and -0.07 ± 0.10 (N = 19, n = 1.8 regime of Fig. 2a), and the drawn lines of Figs. 2a and 4.",
        BibTex_Reference = "
            @article{goldsby2001superplastic,
            title={Superplastic deformation of ice: Experimental observations},
            author={Goldsby, D. L. and Kohlstedt, D. L.},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={106},
            number={B6},
            pages={11017--11030},
            year={2001},
            doi={10.1029/2000JB900336}
            }
    ",
    )
    return data, info
end

@inbounds GrainBoundarySliding_database(f::F) where {F} = first(f())
@inbounds GrainBoundarySliding_database_info(f::F) where {F} = last(f())

end
