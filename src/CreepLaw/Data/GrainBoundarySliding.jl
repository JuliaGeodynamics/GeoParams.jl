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

@inbounds GrainBoundarySliding_database(f::F) where {F} = first(f())
@inbounds GrainBoundarySliding_database_info(f::F) where {F} = last(f())

end
