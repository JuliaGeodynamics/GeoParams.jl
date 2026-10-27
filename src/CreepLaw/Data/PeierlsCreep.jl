module Peierls

using GeoParams
export peierls_database, peierls_database_info


function dry_olivine_Goetze_1979()
    data = PeierlsCreep(;
        Name = "Dry Olivine | Goetze and Evans (1979)",
        n = 1.0NoUnits,                           # power-law exponent (not used)
        o = 1.0NoUnits,                           # inner stress exponent p in (1 - (σ/σP)^p)^q
        q = 2.0NoUnits,                           # outer exponent q
        TauP = 8.5e9Pa,                           # Peierls stress
        A = (5.7e11)s^(-1),                       # material specific rheological parameter
        E = 536.0kJ / mol,                        # activation energy
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Goetze & Evans (1979), eq. (1a), p. 469: ε̇ = 5.7e11 exp(-Q/RT (1 - σ/σP)^2), Q = 128 kcal/mol, σP = 8.5 GPa, σ = σ1 - σ3 (p = 1, q = 2). Reproduces the Goetze curve of Mei et al. (2010), Fig. 6, within 10 MPa (N = 8, 650-1350 K) and the paper's Fig. 4 Peierls curve up to about 1200 °C. Primary data: Evans & Goetze (1979), microindentation hardness of olivine single crystals (not checked here).",
        BibTex_Reference = "
            @article{goetze1979stress,
            title={Stress and temperature in the bending lithosphere as constrained by experimental rock mechanics},
            author={Goetze, Christopher and Evans, Brian},
            journal={Geophysical Journal International},
            volume={59},
            number={3},
            pages={463--478},
            year={1979},
            publisher={Blackwell Publishing Ltd Oxford, UK}
            }
    ",
    )
    return data, info
end

function dry_olivine_Demouchy_2013()
    data = PeierlsCreep(;
        Name = "Dry Olivine | Demouchy (2013)",
        n = 1.0NoUnits,                           # power-law exponent (not used)
        o = 0.5NoUnits,                           # inner stress exponent p in (1 - (σ/σP)^p)^q
        q = 2.0NoUnits,                           # outer exponent q
        TauP = 15.0e9Pa,                          # Peierls stress
        A = (1.0e6)s^(-1),                        # material specific rheological parameter
        E = 450.0kJ / mol,                        # activation energy
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Demouchy et al. (2013), eq. (4) and Table 2 'Model 3': A = 1e6 1/s, Q = 450 kJ/mol, σP = 15 GPa, p = 1/2, q = 2; σ is the differential stress. Reproduces the paper's Fig. 10b 'model 3' curve within 10 MPa (N = 13, 820-1540 K).",
        BibTex_Reference = "
            @article{demouchy2013low,
            title={Low strength of Earth’s uppermost mantle inferred from tri-axial deformation experiments on dry olivine crystals},
            author={Demouchy, Sylvie and Tommasi, Andr{\'e}a and Ballaran, Tiziana Boffa and Cordier, Patrick},
            journal={Physics of the Earth and Planetary Interiors},
            volume={220},
            pages={37--49},
            year={2013},
            publisher={Elsevier}
            }
    ",
    )
    return data, info
end

function dry_olivine_Idrissi_2016()
    data = PeierlsCreep(;
        Name = "Dry Olivine | Idrissi et al. (2016)",
        n = 1.0NoUnits,                           # power-law exponent (not used)
        o = 0.5NoUnits,                           # inner stress exponent p in (1 - (σ/σP)^p)^q
        q = 2.0NoUnits,                           # outer exponent q
        TauP = 3.8e9Pa,                           # Peierls stress
        A = (1.0e6)s^(-1),                        # material specific rheological parameter
        E = 566.0kJ / mol,                        # activation energy
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Idrissi et al. (2016), eq. (1) and Table 1: A = 1e6 1/s, Q = 566 kJ/mol, σ0 = 3.8 GPa, p = 0.5, q = 2 (uniaxial stress). Renamed from dry_olivine_Idrissei_2016. Reproduces the paper's Fig. 1 law at 1e-5 s^-1 within 10 MPa (N = 10, 208-1103 K).",
        BibTex_Reference = "
            @article{idrissi2016low,
            title={Low-temperature plasticity of olivine revisited with in situ TEM nanomechanical testing},
            author={Idrissi, Hosni and Bollinger, Caroline and Boioli, Francesca and Schryvers, Dominique and Cordier, Patrick},
            journal={Science Advances},
            volume={2},
            number={3},
            pages={e1501671},
            year={2016},
            doi={10.1126/sciadv.1501671}
            }
    ",
    )
    return data, info
end

@inbounds peierls_database(f::F) where {F} = first(f())
@inbounds peierls_database_info(f::F) where {F} = last(f())

end
