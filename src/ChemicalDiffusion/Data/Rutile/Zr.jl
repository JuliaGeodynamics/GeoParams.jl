# ------------------------------------- Zr -------------------------------------

"""
    Rt_Zr_Cherniak2007_para_c

Diffusion data of Zr in rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 750-1100°C.
From Cherniak et al. (2007) (https://doi.org/10.1016/j.epsl.2007.06.027).
"""
function Rt_Zr_Cherniak2007_para_c()
    data = DiffusionData(
        Name = "Zr diffusion in Rutile (Ξc and anhydrous conditions) | Cherniak et al. (2007)",
        Mineral = "Rutile",
        Formula = "TiO2",
        Species = "Hf",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Buffer = "QFM, NNO",
        Fluid = "Anhydrous",
        D0 = 9.8e-15m^2 / s,
        Ea = 170kJ / mol,
        Ea_1σ = 30kJ / mol,
        Charge = 4,  # charge of the cation
        T_range_min = 750C,
        T_range_max = 1100C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (20.01.25)",
        BibTex_Reference = "
            @article{cherniak2007zr,
            title={Zr and Hf diffusion in rutile},
            author={Cherniak, DJ and Manchester, J and Watson, EB},
            journal={Earth and Planetary Science Letters},
            volume={261},
            number={1-2},
            pages={267--279},
            year={2007},
            publisher={Elsevier}
            }
            ",
    )

    return data, info
end
