"""
    Zr_Melt_Holycross2018_rhyolitic_highH2O()

Diffusion data of Zr in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function Zr_Melt_Holycross2018_rhyolitic_highH2O()
    return create_Melt_Holycross2018_data(;
        Name = "Zr diffusion in rhyolitic melt (6.2 wt% H2O) | Holycross and Watson (2018)",
        Species = "Zr",
        Buffer = "non-buffered",
        D0 = (10^(-6.36))u"m^2 / s",
        log_D0_1σ = (1.78)NoUnits,
        Ea = (155.16)u"kJ/mol",
        Ea_1σ = (39.92)u"kJ/mol",
        T_range_min = 850C,
        T_range_max = 935C
    )
end

"""
    Zr_Melt_Holycross2018_rhyolitic_mediumH2O()

Diffusion data of Zr in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function Zr_Melt_Holycross2018_rhyolitic_mediumH2O()
    return create_Melt_Holycross2018_data(;
        Name = "Zr diffusion in rhyolitic melt (4.1 wt% H2O) | Holycross and Watson (2018)",
        Species = "Zr",
        Buffer = "NNO",
        D0 = (10^(-5.45))u"m^2 / s",
        log_D0_1σ = (0.26 * 2.303)NoUnits,
        Ea = (182.38)u"kJ/mol",
        Ea_1σ = (6.55)u"kJ/mol",
        T_range_min = 960C,
        T_range_max = 1250C
    )
end
