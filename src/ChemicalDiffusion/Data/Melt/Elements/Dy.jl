"""
    Melt_Dy_Holycross2018_rhyolitic_highH2O()

Diffusion data of Dy in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function Melt_Dy_Holycross2018_rhyolitic_highH2O()
    return create_Melt_Holycross2018_data(
        Name = "Dy diffusion in rhyolitic melt (6.2 wt% H2O) | Holycross and Watson (2018)",
        Species = "Dy",
        Buffer = "non-buffered",
        D0 = (10^(-4.6))u"m^2 / s",
        log_D0_1σ = (1.04 * 2.303)NoUnits,
        Ea = (176.47)u"kJ/mol",
        Ea_1σ = (23.27)u"kJ/mol",
        T_range_min = 850C,
        T_range_max = 935C
    )
end


"""
    Melt_Dy_Holycross2018_rhyolitic_mediumH2O()

Diffusion data of Dy in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function Melt_Dy_Holycross2018_rhyolitic_mediumH2O()
    return create_Melt_Holycross2018_data(
        Name = "Dy diffusion in rhyolitic melt (4.1 wt% H2O) | Holycross and Watson (2018)",
        Species = "Dy",
        Buffer = "NNO",
        D0 = (10^(-4.9))u"m^2 / s",
        log_D0_1σ = (0.32 * 2.303)NoUnits,
        Ea = (181.21)u"kJ/mol",
        Ea_1σ = (8.23)u"kJ/mol",
        T_range_min = 960C,
        T_range_max = 1250C
    )
end
