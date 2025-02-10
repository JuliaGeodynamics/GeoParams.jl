"""
    Rt_Mn_Sasaki1985_para_c

Diffusion data of Mn in rutile. With anhydrous conditions, in air and parallel to c-axis.
Calibrated between 600-1000°C. From Sasaki et al. (1985) (https://doi.org/10.1016/0022-3697(85)90129-5).
"""
function Rt_Mn_Sasaki1985_para_c()
    data = DiffusionData(
        Name = "Mn diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Sasaki et al. (1985)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Mn",  # element or species being diffused
        Orientation = "Ξc",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "air",
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 5.82e-7m^2 / s,  # pre-exponential factor
        Ea = 137kJ / mol,  # activation energy
        Charge = 3,  # charge of the cation
        T_range_min = 600C,  # temperature min of the experiment
        T_range_max = 1000C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Values re-fitted by ML (10.02.25)",
        BibTex_Reference = "
            @article{sasaki1985tracer,
            title={Tracer impurity diffusion in single-crystal rutile (TiO2- x)},
            author={Sasaki, Jun and Peterson, NL and Hoshino, K},
            journal={Journal of Physics and Chemistry of Solids},
            volume={46},
            number={11},
            pages={1267--1283},
            year={1985},
            publisher={Elsevier}
            }
          ",
    )

    return data, info
end
