"""
    Hf_Melt_Holycross2018_rhyolitic_mediumH2O()

Diffusion data of Hf in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function Hf_Melt_Holycross2018_rhyolitic_mediumH2O()
    return create_Melt_Holycross2018_data(
        Name = "Hf diffusion in rhyolitic melt (4.1 wt% H2O) | Holycross and Watson (2018)",
        Species = "Hf",
        Buffer = "NNO",
        D0 = (10^(-3.52))u"m^2 / s",
        log_D0_1σ = (0.2 * 2.303)NoUnits,
        Ea = (231.09)u"kJ/mol",
        Ea_1σ = (5.31)u"kJ/mol",
        T_range_min = 960C,
        T_range_max = 1250C
    )
end
