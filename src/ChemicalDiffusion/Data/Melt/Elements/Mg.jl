"""
    Melt_Mg_Sheng1992_basaltic()

Diffusion data of Mg self-diffusion in basaltic melt. Calibrated with experiments with Spinel and anhydrous synthetic glass between 1250 to 1550°C at ambient pressure. From Sheng et al. (1992) (https://doi.org/10.1016/0016-7037(92)90207-Y).
"""
function Melt_Mg_Sheng1992_basaltic()
    data = DiffusionData(
        Name = "Mg self-diffusion in anhydrous basaltic melt | Sheng and Watson (1992)",
        Phase = "Melt",
        Species = "Mg",
        Buffer = "Air",
        D0 = 7791.9cm^2 / s,
        Ea = 343u"kJ / mol",
        Ea_1σ = (25 / 2)u"kJ / mol",
        Charge = 2,
        T_range_min = 1250C,
        T_range_max = 1500C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (19.03.25)",
        BibTex_Reference = "
            @article{sheng1992self,
            title={Self-diffusion of magnesium in spinel and in equilibrium melts: Constraints on flash heating of silicates},
            author={Sheng, YJ and Wasserburg, GJ and Hutcheon, ID},
            journal={Geochimica et Cosmochimica Acta},
            volume={56},
            number={6},
            pages={2535--2546},
            year={1992},
            publisher={Elsevier}
            }
          ",
    )

    return data, info
end
