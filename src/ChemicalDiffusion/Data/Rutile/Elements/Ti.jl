"""
    Rt_Ti_Hoshino1985_para_c

Diffusion data of Ti self-diffusion in rutile. Parallel c-axis, anhydrous.
Calibrated between 1000-1500C.From Hoshino et al. (1985) (https://doi.org/10.1016/0022-3697(85)90079-4).
"""
function Rt_Ti_Hoshino1985_para_c()
    data = DiffusionData(
        Name = "Ti self-diffusion in Rutile (parallel to c-axis) | Hoshino et al. (1985)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "Ti",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Fluid = "anhydrous",
        D0 = 0.0006500000000000001m^2 / s,
        log_D0_1σ = log(1.33e-4)NoUnits,
        Ea = 276604.24J / mol,
        Ea_1σ = 2343.04J / mol,
        Charge = -2,
        T_range_min = 1000C,
        T_range_max = 1500C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{hoshino1985diffusion,
            title={Diffusion and point defects in TiO2- x},
            author={Hoshino, K and Peterson, NL and Wiley, CL},
            journal={Journal of Physics and Chemistry of Solids},
            volume={46},
            number={12},
            pages={1397--1411},
            year={1985},
            publisher={Elsevier}
            }
          ",
    )

    return data, info
end


"""
    Rt_Ti_Hoshino1985_perp_c

Diffusion data of Ti self-diffusion in rutile. Perpendicular c-axis, anhydrous.
Calibrated between 1000-1500C.From Hoshino et al. (1985) (https://doi.org/10.1016/0022-3697(85)90079-4).
"""
function Rt_Ti_Hoshino1985_perp_c()
    data = DiffusionData(
        Name = "Ti self-diffusion in Rutile (Perpendicular to c-axis) | Hoshino et al. (1985)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "Ti",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "anhydrous",
        D0 = 0.000455m^2 / s,
        log_D0_1σ = log(1.78e-4)NoUnits,
        Ea = 268110.72J / mol,
        Ea_1σ = 4142.16J / mol,
        Charge = -2,
        T_range_min = 1000C,
        T_range_max = 1500C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{hoshino1985diffusion,
            title={Diffusion and point defects in TiO2- x},
            author={Hoshino, K and Peterson, NL and Wiley, CL},
            journal={Journal of Physics and Chemistry of Solids},
            volume={46},
            number={12},
            pages={1397--1411},
            year={1985},
            publisher={Elsevier}
            }
          ",
    )

    return data, info
end
