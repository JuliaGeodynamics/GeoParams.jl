"""
    Rt_Nb_Marschall2013_para_c

Diffusion data of Nb in rutile. With anhydrous conditions and parallel to c-axis. 
Calibrated between 850-1250°C. From Cherniak et al. (2019) 
(https://doi.org/10.1016/j.epsl.2013.05.055).
"""
function Rt_Nb_Marschall2013_para_c()
    data = DiffusionData(
        Name = "Nb diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Cherniak et al. (2019)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Nb",  # element or species being diffused
        Orientation = "Ξc",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "FMQ",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 5.3e-3m^2 / s,  # pre-exponential factor
        log_D0_1σ = 0.40NoUnits,
        Ea = 377.5kJ / mol,  # activation energy
        Ea_1σ = 9.8kJ / mol,  # uncertainty at 1σ of the activation energy
        Charge = 5,  # charge of the cation
        T_range_min = 850C,  # temperature min of the experiment
        T_range_max = 1250C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (15.01.25)",
        BibTex_Reference = "
            @article{Marschall2013,
            abstract = {Differentiation of the Earth into its major spheres - crust, mantle and core - has proceeded dominantly through magmatic processes involving melting and melt separation. Models that describe these differentiation processes are guided by elemental abundances in the different reservoirs. Elements are fractionated between coexisting phases during partial melting, and geochemical models are generally based on the fundamental assumption that trace-element equilibrium is established between the partial melts and the restitic minerals. The element pair niobium and tantalum is key to the distinction of different melting regimes involved in crustal differentiation, but equilibrium partition models have largely failed to reproduce the Nb/Ta patterns observed in nature, posing a long-standing geochemical conundrum. Here we demonstrate that kinetic fractionation of Nb and Ta by diffusion may have produced the low Nb/Ta observed in the continental crust. On the basis of the diffusivities of Nb and Ta in rutile (TiO2) determined experimentally in this study, we conclude that equilibrium cannot be expected for the natural range of grain sizes, temperatures and time scales involved in partial melting of crustal rocks. Instead, the observed fractionation of the geochemical twins, Nb and Ta, in the silicate Earth most likely proceeds by partial - as opposed to complete - equilibration of rutile and melt. Hence, the assumption of bulk equilibrium during partial melting for the processes of crustal differentiation may not be justified, as is demonstrated here for Nb/Ta. The concept presented here is based on kinetic fractionation melting and explains the observed low Nb/Ta ratio of the continental crust. © 2013 Elsevier B.V.},
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