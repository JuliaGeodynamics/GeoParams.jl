"""
    Rt_O_Arita1979_para_c

Diffusion data of O in rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 1150-1450K.
From Arita et al. (1979) (https://doi.org/10.1111/j.1151-2916.1979.tb19101.x).
"""
function Rt_O_Arita1979_para_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Arita et al. (1979)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        D0 = 3.4e-7m^2 / s,
        Ea = 251kJ / mol,
        Ea_1σ = 75.3kJ / mol,
        Charge = -2,
        T_range_min = 1150K,
        T_range_max = 1450K,
        P0 = 6000Pa
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Arita1979,
            author = {M. Arita and M. Hosoya and M. Kobayashi and M. Someno},
            doi = {10.1111/j.1151-2916.1979.tb19101.x},
            issn = {15512916},
            issue = {9-10},
            journal = {Journal of the American Ceramic Society},
            pages = {443-446},
            title = {Depth Profile Measurement by Secondary Ion Mass Spectrometry for Determining the Tracer Diffusivity of Oxygen in Rutile},
            volume = {62},
            year = {1979},
        }

          ",
    )

    return data, info
end

"""
    Rt_O_Arita1979_para_c_Cr

Diffusion data of O in Cr-doped rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 1150-1450K.
From Arita et al. (1979) (https://doi.org/10.1111/j.1151-2916.1979.tb19101.x).
"""
function Rt_O_Arita1979_para_c_Cr()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Arita et al. (1979)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Doping = "Cr",
        D0 = 2.0e-8m^2 / s,
        Ea = 204kJ / mol,
        Ea_1σ = 61.2kJ / mol,
        Charge = -2,
        T_range_min = 1150K,
        T_range_max = 1450K,
        P0 = 6000Pa
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Arita1979,
            author = {M. Arita and M. Hosoya and M. Kobayashi and M. Someno},
            doi = {10.1111/j.1151-2916.1979.tb19101.x},
            issn = {15512916},
            issue = {9-10},
            journal = {Journal of the American Ceramic Society},
            pages = {443-446},
            title = {Depth Profile Measurement by Secondary Ion Mass Spectrometry for Determining the Tracer Diffusivity of Oxygen in Rutile},
            volume = {62},
            year = {1979},
        }

          ",
    )

    return data, info
end


"""
    Rt_O_Dennis1993_perp_c

Diffusion data of O in rutile. With hydrothermal conditions and perpendicular to c-axis.
Calibrated between 873-1373K. From Dennis et al. (1993) (https://doi.org/10.1007/BF00414275).
"""
function Rt_O_Dennis1993_perp_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (perpendicular to c-axis and hydrothermal conditions) | Dennis et al. (1993)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "hydrothermal",
        D0 = 2.41e-12m^2 / s,
        log_D0_1σ = 25.4NoUnits,
        Ea = 172.5kJ / mol,
        Ea_1σ = 23.6kJ / mol,
        Charge = -2,
        T_range_min = 873K,
        T_range_max = 1373K,
        P0 = 100MPa
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Dennis1993,
            author = {P. F. Dennis and R. Freer},
            doi = {10.1007/BF00414275},
            issn = {00222461},
            issue = {17},
            journal = {Journal of Materials Science},
            pages = {4804-4810},
            title = {Oxygen self-diffusion in rutile under hydrothermal conditions},
            volume = {28},
            year = {1993},
        }
          ",
    )

    return data, info
end


"""
    Rt_O_Dennis1993_perp_c_nat

Diffusion data of O in natural rutile. With hydrothermal conditions and perpendicular to c-axis.
Calibrated between 873-1373K. From Dennis et al. (1993) (https://doi.org/10.1007/BF00414275).
"""
function Rt_O_Dennis1993_perp_c_nat()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (perpendicular to c-axis and hydrothermal conditions) | Dennis et al. (1993)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "hydrothermal",
        Doping = "natural",
        D0 = 1.14e-11m^2 / s,
        log_D0_1σ = 23.1NoUnits,
        Ea = 168.8kJ / mol,
        Ea_1σ = 32.9kJ / mol,
        Charge = -2,
        T_range_min = 873K,
        T_range_max = 1373K,
        P0 = 100MPa
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Dennis1993,
            author = {P. F. Dennis and R. Freer},
            doi = {10.1007/BF00414275},
            issn = {00222461},
            issue = {17},
            journal = {Journal of Materials Science},
            pages = {4804-4810},
            title = {Oxygen self-diffusion in rutile under hydrothermal conditions},
            volume = {28},
            year = {1993},
        }
          ",
    )

    return data, info
end


"""
    Rt_O_Derry1981_para_c

Diffusion data of O in natural rutile. With parallel to c-axis.
Calibrated between 1173-1673K. From Derry et al. (1981) (https://doi.org/10.1016/0022-3697(81)90011-1).
"""
function Rt_O_Derry1981_para_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis) | Derry et al. (1981)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        D0 = 2.4e-6m^2 / s,
        Ea = 282.6kJ / mol,
        Ea_1σ = 5kJ / mol,
        Charge = -2,
        T_range_min = 1173K,
        T_range_max = 1673K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Derry1981,
            author = {D. J. Derry and D. G. Lees and J. M. Calvert},
            doi = {https://doi.org/10.1016/0022-3697(81)90011-1},
            journal = {J. Phys. Chem. Solids},
            pages = {57-64},
            title = {A study of oxygen self-diffusion in the c-direction of rutile using a nuclear technique},
            volume = {42},
            year = {1981},
            }
          ",
    )

    return data, info
end

"""
    Rt_O_Haul1965_unor

Diffusion data of O in rutile. Unoriented, anhydrous.
Calibrated between 710-1300C. From Haul et al. (1965) (https://doi.org/10.1016/0022-3697(65)90066-1).
"""
function Rt_O_Haul1965_unor()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (Unoriented) | Haul et al. (1965)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Unoriented",
        Crystallography = "Tetragonal",
        Fluid = "anhydrous",
        D0 = 2.0e-3cm^2 / s,
        Ea = (60.0e3)u"J / mol",
        Ea_1σ = (1.5e3)u"J / mol",
        Charge = -2,
        T_range_min = 710C,
        T_range_max = 1300C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{haul1965sauerstoff,
            title={Sauerstoff-selbstdiffusion in rutilkristallen},
            author={Haul, R and D{\"u}mbgen, G},
            journal={Journal of Physics and Chemistry of Solids},
            volume={26},
            number={1},
            pages={1--10},
            year={1965},
            publisher={Elsevier}
            }
          ",
    )

    return data, info
end


"""
    Rt_O_Lundy1973_para_c

Diffusion data of O in natural rutile. With parallel to c-axis.
Calibrated between 1200-1500K. From Lundy et al. (1973) (https://doi.org/10.1051/jphyscol:1973953).
"""
function Rt_O_Lundy1973_para_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis) | Lundy et al. (1973)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        D0 = 0.046cm^2 / s,
        Ea = 59.9u"cal/mol",
        Charge = -2,
        T_range_min = 1200K,
        T_range_max = 1500K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{lundy1973processus,
            title={Processus de transport dans les oxydes cation self diffusion in rutile},
            author={Lundy, TS and Coghlan, WA},
            journal={Le Journal de Physique Colloques},
            volume={34},
            number={C9},
            pages={C9--299},
            year={1973}
            }
          ",
    )

    return data, info
end


"""
    Rt_O_Lundy1973_perp_c

Diffusion data of O in natural rutile. With perpendicular to c-axis.
Calibrated between 1200-1500K. From Lundy et al. (1973) (https://doi.org/10.1051/jphyscol:1973953).
"""
function Rt_O_Lundy1973_perp_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (perpendicular to c-axis) | Lundy et al. (1973)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        D0 = 0.0024cm^2 / s,
        Ea = 48.5u"cal/mol",
        Charge = -2,
        T_range_min = 1200K,
        T_range_max = 1500K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{lundy1973processus,
            title={Processus de transport dans les oxydes cation self diffusion in rutile},
            author={Lundy, TS and Coghlan, WA},
            journal={Le Journal de Physique Colloques},
            volume={34},
            number={C9},
            pages={C9--299},
            year={1973}
            }
          ",
    )

    return data, info
end


"""
    Rt_O_Moore1998_para_c_fast

Diffusion data of O in natural rutile. With hydrous and anhydrous conditions, parallel to c-axis, at 0.1 to 1000 MPa.
Calibrated between 1200-1500K. From Moore et al. (1998) (https://doi.org/10.2138/am-1998-7-803).
"""
function Rt_O_Moore1998_para_c_fast()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis) | Moore et al. (1998)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Buffer = "1atm to NNO",
        D0 = 4.7e-7m^2 / s,
        Ea = 258.0e3J / mol,
        Ea_1σ = 22.0e3J / mol,
        Charge = -2,
        T_range_min = 1200K,
        T_range_max = 1500K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{Moore1998,
            author = {D. K. Moore and Daniele J. Cherniak and E. B. Watson},
            doi = {10.2138/am-1998-7-803},
            issn = {0003004X},
            journal = {American Mineralogist},
            pages = {700-711},
            title = {Oxygen diffusion in rutile from 750 to 1000 °C and 0.1 to 1000 MPa},
            volume = {83},
            year = {1998},
            }
          ",
    )

    return data, info
end


"""
    Rt_O_Moore1998_para_c_slow

Diffusion data of O in natural rutile. With hydrous and anhydrous conditions, parallel to c-axis, at 0.1 to 1000 MPa.
Calibrated between 1200-1500K. From Moore et al. (1998) (https://doi.org/10.2138/am-1998-7-803).
"""
function Rt_O_Moore1998_para_c_slow()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis) | Moore et al. (1998)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Buffer = "1atm to NNO",
        D0 = 5.9e-5m^2 / s,
        Ea = 330.0e3J / mol,
        Ea_1σ = 15.0e3J / mol,
        Charge = -2,
        T_range_min = 1200K,
        T_range_max = 1500K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{Moore1998,
            author = {D. K. Moore and Daniele J. Cherniak and E. B. Watson},
            doi = {10.2138/am-1998-7-803},
            issn = {0003004X},
            journal = {American Mineralogist},
            pages = {700-711},
            title = {Oxygen diffusion in rutile from 750 to 1000 °C and 0.1 to 1000 MPa},
            volume = {83},
            year = {1998},
            }
          ",
    )

    return data, info
end

"""
    Rt_O_Venkatu1970

Diffusion data of O in natural rutile.
Calibrated between 900-1300C. From Venkatu et al. (1970).
"""
function Rt_O_Venkatu1970()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile | Venkatu et al. (1970)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Crystallography = "Tetragonal",
        D0 = 6.4e-2cm^2 / s,
        Ea = 61400u"cal/mol",
        Ea_1σ = 22.0e3J / mol,
        Charge = -2,
        T_range_min = 900C,
        T_range_max = 1300C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{Venkatu1970,
            author = {D.A. Venkatu and L.E. Poteat},
            journal = {Material Science and Engineering},
            pages = {258-262},
            title = {Diffusion of Titanium in Single Crystal Rutile},
            volume = {5},
            year = {1970},
            }
          ",
    )

    return data, info
end
