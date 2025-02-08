"""
    Rt_Li_Johnson1964_perp_c

Diffusion data of Li in rutile. With anhydrous conditions and perpendicular to c-axis.
Calibrated between 80-360°C.
From Johnson et al. (1964) (https://doi.org/10.1103/PhysRev.136.A284).
"""
function Rt_Li_Johnson1964_perp_c()
    data = DiffusionData(
        Name = "Li diffusion in Rutile (perpendicular to c-axis and anhydrous conditions) | Johnson et al. (1964)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "Li",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "Anhydrous",
        D0 = 0.295cm^2 / s,
        log_D0_1σ = log(0.028)NoUnits,
        Ea = (0.33u"eV" * Unitful.Na),
        Ea_1σ = (0.003u"eV" * Unitful.Na),
        Charge = 1,  # charge of the cation
        T_range_min = 80C,
        T_range_max = 360C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Johnson1964,
            author = {O. W. Johnson},
            doi = {https://doi.org/10.1103/PhysRev.136.A284},
            issue = {1A},
            journal = {Physical Review},
            pages = {284-290},
            title = {One-Dimensional Diffusion of Li in Rutile},
            volume = {136},
            year = {1964},
            }
          ",
    )

    return data, info
end
