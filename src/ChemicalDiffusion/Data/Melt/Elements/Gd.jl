"""
    Gd_Melt_Holycross2018_rhyolitic_highH2O()

Diffusion data of Gd in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function Gd_Melt_Holycross2018_rhyolitic_highH2O()
    return create_Melt_Holycross2018_data(
        Name = "Gd diffusion in rhyolitic melt (6.2 wt% H2O) | Holycross and Watson (2018)",
        Species = "Gd",
        Buffer = "non-buffered",
        D0 = (10^(-4.06))u"m^2 / s",
        log_D0_1σ = (0.97 * 2.303)NoUnits,
        Ea = (187.25)u"kJ/mol",
        Ea_1σ = (21.66)u"kJ/mol",
        T_range_min = 850C,
        T_range_max = 935C
    )
end


"""
    Gd_Melt_Holycross2018_rhyolitic_mediumH2O()

Diffusion data of Gd in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function Gd_Melt_Holycross2018_rhyolitic_mediumH2O()
    return create_Melt_Holycross2018_data(
        Name = "Gd diffusion in rhyolitic melt (4.1 wt% H2O) | Holycross and Watson (2018)",
        Species = "Gd",
        Buffer = "NNO",
        D0 = (10^(-4.57))u"m^2 / s",
        log_D0_1σ = (0.41 * 2.303)NoUnits,
        Ea = (188.24)u"kJ/mol",
        Ea_1σ = (10.74)u"kJ/mol",
        T_range_min = 960C,
        T_range_max = 1250C
    )
end
