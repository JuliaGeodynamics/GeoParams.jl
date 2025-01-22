# ------------------------------------- Hf -------------------------------------

"""
    Rt_Hf_Cherniak2007_perp_c

Diffusion data of Hf in rutile. With anhydrous conditions and perpendicular to c-axis. Calibrated between 750-1050°C. From Cherniak et al. (2007) (https://doi.org/10.1016/j.epsl.2007.06.027).
"""
function Rt_Hf_Cherniak2007_perp_c()
    data = DiffusionData(
        Name = "Hf diffusion in Rutile (perpendicular to c-axis and anhydrous conditions) | Cherniak et al. (2007)",
        Mineral = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Hf",  # element or species being diffused
        D0 = 2.5e-12m^2 / s,  # pre-exponential factor
        Ea = 227kJ / mol,  # activation energy
        Ea_1σ = 62kJ / mol,  # uncertainty at 1σ of the activation energy
        T_range_min = 750C,  # temperature min of the experiment
        T_range_max = 1050C,  # temperature max of the experiment
        Orientation = "⊥c",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "QFM, NNO",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Anhydrous"  # Fluid condition (e.g., anhydrous) during the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (15.01.25)",
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


"""
    Rt_Hf_Cherniak2007_para_c

Diffusion data of Hf in rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 800-1000°C.
From Cherniak et al. (2007) (https://doi.org/10.1016/j.epsl.2007.06.027).
"""
function Rt_Hf_Cherniak2007_para_c()
    data = DiffusionData(
        Name = "Hf diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Cherniak et al. (2007)",
        Mineral = "Rutile",
        Formula = "TiO2",
        Species = "Hf",
        D0 = 9.1e-15m^2 / s,
        Ea = 169kJ / mol,
        Ea_1σ = 36kJ / mol,
        T_range_min = 800C,
        T_range_max = 1000C,
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Buffer = "QFM, NNO",
        Fluid = "Anhydrous"
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (15.01.25)",
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
