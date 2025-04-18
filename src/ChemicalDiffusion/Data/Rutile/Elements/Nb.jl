"""
    Rt_Nb_Marschall2013_para_c

Diffusion data of Nb in rutile. With anhydrous conditions and parallel to c-axis.
Calibrated between 850-1250°C. From Marschall et al. (2013) (https://doi.org/10.1016/j.epsl.2013.05.055).
"""
function Rt_Nb_Marschall2013_para_c()
    data = DiffusionData(
        Name = "Nb diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Marschall et al. (2013)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Nb",  # element or species being diffused
        Orientation = "Ξc",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "FMQ",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 5.3e-3m^2 / s,  # pre-exponential factor
        log_D0_1σ = 0.4NoUnits,
        Ea = 377.5kJ / mol,  # activation energy
        Ea_1σ = 9.8kJ / mol,  # uncertainty at 1σ of the activation energy
        Charge = 5,  # charge of the cation
        T_range_min = 850C,  # temperature min of the experiment
        T_range_max = 1250C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{marschall2013diffusion,
            title={Diffusion-induced fractionation of niobium and tantalum during continental crust formation},
            author={Marschall, Horst R and Dohmen, Ralf and Ludwig, Thomas},
            journal={Earth and Planetary Science Letters},
            volume={375},
            pages={361--371},
            year={2013},
            publisher={Elsevier}
            }
          ",
    )

    return data, info
end
