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
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'dry GBS, T < 1250 °C' (A in MPa^-3.5 μm^2/s, differential stress). V as for dry dislocation creep (footnote g; 18e-6 is the thermocouple-corrected value of Table 2). HK03 Fig. 7a draws the GBS line with p ≈ 1.4 while Table 1 gives p = 2; the law reproduces the line near d ≈ 15 μm (log10 misfit +0.09, -0.01, -0.34 at d = 10, 15, 60 μm).",
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
        Comment = "Hirth & Kohlstedt (2003), Table 1, row 'dry GBS, T > 1250 °C' (A in MPa^-3.5 μm^2/s, differential stress). V as for dry dislocation creep (footnote g; 18e-6 is the thermocouple-corrected value of Table 2).",
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
