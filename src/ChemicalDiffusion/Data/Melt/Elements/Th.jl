"""
    Melt_Th_Holycross2018_rhyolitic_highH2O()

Diffusion data of Th in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function Melt_Th_Holycross2018_rhyolitic_highH2O()
    return create_Melt_Holycross2018_data(;
        Name = "Th diffusion in rhyolitic melt (6.2 wt% H2O) | Holycross and Watson (2018)",
        Species = "Th",
        Buffer = "non-buffered",
        D0 = (10^(-5.1))u"m^2 / s",
        log_D0_1σ = (1.11 * 2.303)NoUnits,
        Ea = (176.71)u"kJ/mol",
        Ea_1σ = (24.86)u"kJ/mol",
        T_range_min = 850C,
        T_range_max = 935C
    )
end
