"""
    Rt_Ta_Marschall2013_para_c

Diffusion data of Ta in rutile. With anhydrous conditions and parallel to c-axis. 
Calibrated between 850-1250°C. From Cherniak et al. (2019) 
(https://doi.org/10.1016/j.epsl.2013.05.055).
"""
function Rt_Ta_Marschall2013_para_c()
    data = DiffusionData(
        Name = "Ta diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Cherniak et al. (2019)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Ta",  # element or species being diffused
        Orientation = "Ξc",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "FMQ",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 6.3e-3m^2 / s,  # pre-exponential factor
        log_D0_1σ = 1.42NoUnits,
        Ea = 392kJ / mol,  # activation energy
        Ea_1σ = 36kJ / mol,  # uncertainty at 1σ of the activation energy
        Charge = 5,  # charge of the cation
        T_range_min = 850C,  # temperature min of the experiment
        T_range_max = 1250C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (15.01.25)",
        BibTex_Reference = "
            @article{Marschall2013,
            author = {Horst R. Marschall and Ralf Dohmen and Thomas Ludwig},
            doi = {10.1016/j.epsl.2013.05.055},
            issn = {0012821X},
            journal = {Earth and Planetary Science Letters},
            keywords = {Diffusion,Nb-Ta,Partial melting,Rutile,TiO2},
            month = {8},
            pages = {361-371},
            publisher = {Elsevier},
            title = {Diffusion-induced fractionation of niobium and tantalum during continental crust formation},
            volume = {375},
            url = {http://dx.doi.org/10.1016/j.epsl.2013.05.055 https://linkinghub.elsevier.com/retrieve/pii/S0012821X13003166},
            year = {2013},
            }
          ",
    )

    return data, info
end
