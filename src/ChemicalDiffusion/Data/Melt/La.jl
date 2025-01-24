
"""
    Melt_La_Holycross2018_rhyolitic_highH2O()

Diffusion data of La in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-1250°C at 1 GPa with Ni and Ag capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function Melt_La_Holycross2018_rhyolitic_highH2O()
    data = DiffusionData(
        Name = "La diffusion in rhyolitic melt (6.2 wt%) | Holycross and Watson (2018)",
        Phase = "Melt",  # name of the mineral
        Formula = "",  # chemical formula of the mineral
        Species = "REE",  # element or species being diffused
        Orientation = "Amorphous",  # Crystal orientation from the diffusion experiment
        Crystallography = "Amorphous",  # Crystallographic system of the mineral
        Buffer = "NNO, non-buffered",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Hydrated",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = (3.571284964163521e-5)u"m^2 / s",  # pre-exponential factor (exp(-9.28) in the original paper)
        Ea = (221057/2.303)u"J/mol",  # activation energy
        Ea_1σ = (4284/2.303)u"J/mol",  # uncertainty at 1σ of the activation energy
        Charge = 3,  # charge of the cation
        T_range_min = 950C,  # temperature min of the experiment
        T_range_max = 1050C,  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (23.01.25). Values originally reported in log10 but converted here to ln.",
        BibTex_Reference = "
            @article{bloch2020multispecies,
            title={Multispecies diffusion of yttrium, rare earth elements and hafnium in garnet},
            author={Bloch, EM and Jollands, MC and Devoir, A and Bouvier, A-S and Ibañez-Mejia, M and Baumgartner, LP},
            journal={Journal of Petrology},
            volume={61},
            number={7},
            pages={egaa055},
            year={2020},
            publisher={Oxford University Press}
            }
            ",
    )

    return data, info
end