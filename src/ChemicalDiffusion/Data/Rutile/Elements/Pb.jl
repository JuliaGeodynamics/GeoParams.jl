"""
    Rt_Pb_Cherniak2000_unor

Diffusion data of Pb in rutile. With anhydrous conditions and unoriented. 
Calibrated between 700-1100°C. From Cherniak et al. (2000) 
(https://doi.org/10.1007/PL00007671).
"""
function Rt_Pb_Cherniak2000_unor()
    data = DiffusionData(
        Name = "Pb diffusion in Rutile (unoriented and anhydrous conditions) | Cherniak et al. (2000)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Pb",  # element or species being diffused
        Orientation = "unoriented",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "QFM, NNO, self-buffered",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 3.9e-10m^2 / s,  # pre-exponential factor
        log_D0_1σ = 1.04NoUnits,
        Ea = 250kJ / mol,  # activation energy
        Ea_1σ = 12kJ / mol,  # uncertainty at 1σ of the activation energy
        Charge = 2,  # charge of the cation
        T_range_min = 700C,  # temperature min of the experiment
        T_range_max = 1100C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (15.01.25)",
        BibTex_Reference = "
          @article{Cherniak2000,
            author = {D. J. Cherniak},
            doi = {10.1007/PL00007671},
            issn = {00107999},
            issue = {2},
            journal = {Contributions to Mineralogy and Petrology},
            pages = {198-207},
            title = {Pb diffusion in rutile},
            volume = {139},
            year = {2000},
            }
          ",
    )

    return data, info
end
