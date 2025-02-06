"""
    Rt_3H_Caskey1974_para_c

Diffusion data of Tritium in rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 155-300°C.
From Caskey et al. (1974) (https://doi.org/10.1016/0025-5416(74)90003-2).
"""
function Rt_3H_Caskey1974_para_c()
    data = DiffusionData(
        Name = "Tritium diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Caskey et al. (1974)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "3H",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Fluid = "Anhydrous",
        D0 = 7.5e-6cm^2 / s,
        Ea = 9040u"cal / mol",
        Charge = 1,  # charge of the cation
        T_range_min = 155C,
        T_range_max = 300C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{caskey1974diffusion,
            title={Diffusion of tritium in rutile (TiO2)},
            author={Caskey Jr, GR},
            journal={Materials Science and Engineering},
            volume={14},
            number={2},
            pages={109--114},
            year={1974},
            publisher={Elsevier}
            }
          ",
    )

    return data, info
end


"""
    Rt_3H_Caskey1974_perp_c

Diffusion data of Tritium in rutile. With anhydrous conditions and perpendicular to c-axis. Calibrated between 155-300°C.
From Caskey et al. (1974) (https://doi.org/10.1016/0025-5416(74)90003-2).
"""
function Rt_3H_Caskey1974_perp_c()
    data = DiffusionData(
        Name = "Tritium diffusion in Rutile (perpendicular to c-axis and anhydrous conditions) | Caskey et al. (1974)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "3H",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "Anhydrous",
        D0 = 2.7e-6cm^2 / s,
        Ea = 13100u"cal / mol",
        Charge = 1,  # charge of the cation
        T_range_min = 155C,
        T_range_max = 300C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{caskey1974diffusion,
            title={Diffusion of tritium in rutile (TiO2)},
            author={Caskey Jr, GR},
            journal={Materials Science and Engineering},
            volume={14},
            number={2},
            pages={109--114},
            year={1974},
            publisher={Elsevier}
            }
          ",
    )

    return data, info
end


"""
    Rt_3H_Cathcart1979_perp_c

Diffusion data of Tritium in rutile. With anhydrous conditions and perpendicular to c-axis. Calibrated between 250-900°C.
From Cathcart et al. (1979) (https://doi.org/110.1063/1.326490).
"""
function Rt_3H_Cathcart1979_perp_c()
    data = DiffusionData(
        Name = "Tritium diffusion in Rutile (perpendicular to c-axis and anhydrous conditions) | Cathcart et al. (1979)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "3H",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "Anhydrous",
        D0 = 1.77e-2cm^2 / s,
        log_D0_1σ = 4.3NoUnits,
        Ea = 107kJ / mol,
        Ea_1σ = 4kJ / mol,
        Charge = 1,  # charge of the cation
        T_range_min = 250C,
        T_range_max = 900C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{cathcart1979tritium,
            title={Tritium diffusion in rutile (TiO2)},
            author={Cathcart, JV and Perkins, RA and Bates, JB and Manley, LC},
            journal={Journal of Applied Physics},
            volume={50},
            number={6},
            pages={4110--4119},
            year={1979},
            publisher={American Institute of Physics}
            }
          ",
    )

    return data, info
end

"""
    Rt_3H_Cathcart1979_para_c

Diffusion data of Tritium in rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 250-900°C.
From Cathcart et al. (1979) (https://doi.org/110.1063/1.326490).
"""
function Rt_3H_Cathcart1979_para_c()
    data = DiffusionData(
        Name = "Tritium diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Cathcart et al. (1979)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "3H",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Fluid = "Anhydrous",
        D0 = 0.0085cm^2 / s,
        log_D0_1σ = 4.1NoUnits,
        Ea = 72.2kJ / mol,
        Ea_1σ = 6.2kJ / mol,
        Charge = 1,  # charge of the cation
        T_range_min = 250C,
        T_range_max = 900C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{cathcart1979tritium,
            title={Tritium diffusion in rutile (TiO2)},
            author={Cathcart, JV and Perkins, RA and Bates, JB and Manley, LC},
            journal={Journal of Applied Physics},
            volume={50},
            number={6},
            pages={4110--4119},
            year={1979},
            publisher={American Institute of Physics}
            }
          ",
    )

    return data, info
end


"""
    Rt_H_Johnson1975_para_c

Diffusion data of H in rutile. With anhydrous conditions and parallel to c-axis. 
Calibrated between 125-750°C.
From Johnson et al. (1975) (https://doi.org/10.1063/1.322206).
"""
function Rt_H_Johnson1975_para_c()
    data = DiffusionData(
        Name = "H diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Johnson et al. (1975)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "H",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Fluid = "Anhydrous",
        D0 = 1.8e-3cm^2 / s,
        log_D0_1σ = log(0.8e-3)NoUnits,
        Ea = (0.59u"eV" * Unitful.Na),
        Ea_1σ = (0.02u"eV" * Unitful.Na),
        Charge = 1,  # charge of the cation
        T_range_min = 125C,
        T_range_max = 750C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Johnson1975,
            author = {O. W. Johnson and S. H. Paek and J. W. Deford},
            doi = {10.1063/1.322206},
            issn = {00218979},
            issue = {3},
            journal = {Journal of Applied Physics},
            pages = {1026-1033},
            title = {Diffusion of H and D in TiO2: Suppression of internal fields by isotope exchange},
            volume = {46},
            year = {1975},
            }
          ",
    )

    return data, info
end


"""
    Rt_H_Johnson1975_perp_c

Diffusion data of H in rutile. With anhydrous conditions and perpendicular to c-axis. 
Calibrated between 125-750°C. From Johnson et al. (1975) (https://doi.org/10.1063/1.322206).
"""
function Rt_H_Johnson1975_perp_c()
    data = DiffusionData(
        Name = "H diffusion in Rutile (perpendicular to c-axis and anhydrous conditions) | Johnson et al. (1975)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "H",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "Anhydrous",
        D0 = 3.8e-1cm^2 / s,
        log_D0_1σ = log(2.0e-1)NoUnits,
        Ea = (1.28u"eV" * Unitful.Na),
        Ea_1σ = (0.05u"eV" * Unitful.Na),
        Charge = 1,  # charge of the cation
        T_range_min = 125C,
        T_range_max = 750C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Johnson1975,
            author = {O. W. Johnson and S. H. Paek and J. W. Deford},
            doi = {10.1063/1.322206},
            issn = {00218979},
            issue = {3},
            journal = {Journal of Applied Physics},
            pages = {1026-1033},
            title = {Diffusion of H and D in TiO2: Suppression of internal fields by isotope exchange},
            volume = {46},
            year = {1975},
            }
          ",
    )

    return data, info
end
