"""
    U_Melt_Holycross2018_rhyolitic_highH2O()

Diffusion data of U in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function U_Melt_Holycross2018_rhyolitic_highH2O()
    return create_Melt_Holycross2018_data(;
        Name = "U diffusion in rhyolitic melt (6.2 wt% H2O) | Holycross and Watson (2018)",
        Species = "U",
        Buffer = "non-buffered",
        D0 = (10^(-1.06))u"m^2 / s",
        log_D0_1σ = (1.45 * 2.303)NoUnits,
        Ea = (267.76)u"kJ/mol",
        Ea_1σ = (32.4)u"kJ/mol",
        T_range_min = 850C,
        T_range_max = 935C
    )
end

"""
    U_Melt_Holycross2018_rhyolitic_mediumH2O()

Diffusion data of U in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass.
From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).
"""
function U_Melt_Holycross2018_rhyolitic_mediumH2O()
    return create_Melt_Holycross2018_data(;
        Name = "U diffusion in rhyolitic melt (4.1 wt% H2O) | Holycross and Watson (2018)",
        Species = "U",
        Buffer = "NNO",
        D0 = (10^(-3.42))u"m^2 / s",
        log_D0_1σ = (0.23 * 2.303)NoUnits,
        Ea = (228.82)u"kJ/mol",
        Ea_1σ = (6.04)u"kJ/mol",
        T_range_min = 960C,
        T_range_max = 1250C
    )
end
