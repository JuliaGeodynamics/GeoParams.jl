# ------------------------------------- Zr -------------------------------------

"""
    Rt_Zr_Cherniak2007_para_c

Diffusion data of Zr in rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 750-1100°C.
From Cherniak et al. (2007) (https://doi.org/10.1016/j.epsl.2007.06.027).
"""
function Rt_Zr_Cherniak2007_para_c()
    data = DiffusionData(
        Name = "Zr diffusion in Rutile (Ξc and anhydrous conditions) | Cherniak et al. (2007)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "Zr",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Buffer = "QFM, NNO, Air",
        Fluid = "Anhydrous",
        D0 = 9.8e-15m^2 / s,
        log_D0_1σ = 1.25NoUnits,
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

"""
    Rt_Zr_Sasaki1985_para_c

Diffusion data of Zr in rutile. With anhydrous conditions, in air and parallel to c-axis.
Calibrated between 1100-1500°C. From Sasaki et al. (1985) (https://doi.org/10.1016/0022-3697(85)90129-5).
"""
function Rt_Zr_Sasaki1985_para_c()
    data = DiffusionData(
        Name = "Zr diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Sasaki et al. (1985)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Zr",  # element or species being diffused
        Orientation = "Ξc",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "air",
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 2.29e-7m^2 / s,  # pre-exponential factor
        Ea = 291kJ / mol,  # activation energy
        Charge = 4,  # charge of the cation
        T_range_min = 1100C,  # temperature min of the experiment
        T_range_max = 1500C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Values re-fitted by ML (10.02.25)",
        BibTex_Reference = "
            @article{Sasaki1985,
            author = {Jun Sasaki and N. L. Peterson and K. Hoshino},
            doi = {10.1016/0022-3697(85)90129-5},
            issn = {00223697},
            issue = {11},
            journal = {Journal of Physics and Chemistry of Solids},
            keywords = {diffusion mechanisms,interstitialcy diffusion mechanism,point defects,rutile (TiO2),tracer impurity diffusion},
            pages = {1267-1283},
            title = {Tracer impurity diffusion in single-crystal rutile (TiO2-x)},
            volume = {46},
            year = {1985},
            }
          ",
    )

    return data, info
end


"""
    Rt_Zr_Sasaki1985_perp_c

Diffusion data of Zr in rutile. With anhydrous conditions, in air and perpendicular to c-axis.
Calibrated between 1100-1500°C. From Sasaki et al. (1985) (https://doi.org/10.1016/0022-3697(85)90129-5).
"""
function Rt_Zr_Sasaki1985_perp_c()
    data = DiffusionData(
        Name = "Zr diffusion in Rutile (perpendicular to c-axis and anhydrous conditions) | Sasaki et al. (1985)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Zr",  # element or species being diffused
        Orientation = "⊥c",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "air",
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 8.38e-7m^2 / s,  # pre-exponential factor
        Ea = 288kJ / mol,  # activation energy
        Charge = 4,  # charge of the cation
        T_range_min = 1100C,  # temperature min of the experiment
        T_range_max = 1500C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Values re-fitted by ML (10.02.25)",
        BibTex_Reference = "
            @article{Sasaki1985,
            author = {Jun Sasaki and N. L. Peterson and K. Hoshino},
            doi = {10.1016/0022-3697(85)90129-5},
            issn = {00223697},
            issue = {11},
            journal = {Journal of Physics and Chemistry of Solids},
            keywords = {diffusion mechanisms,interstitialcy diffusion mechanism,point defects,rutile (TiO2),tracer impurity diffusion},
            pages = {1267-1283},
            title = {Tracer impurity diffusion in single-crystal rutile (TiO2-x)},
            volume = {46},
            year = {1985},
            }
          ",
    )

    return data, info
end
