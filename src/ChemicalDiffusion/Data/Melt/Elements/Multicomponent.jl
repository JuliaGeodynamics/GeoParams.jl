"""
    Melt_multicomponent_major_Guo2020_SiO2_basaltic()

Multicomponent diffusion data of major elements in basaltic melt. Calibrated using 18 diffusion couple experiments for SiO2–TiO2–Al2O3–FeO–MgO–CaO–Na2O–K2O basaltic melt compositions and data from Guo and Zhang (2018). SiO2 is here the dependent variable. This contains the invariant eigenvectors and the diagonal matrices containing the pre-exponential factors and activation energies of the eigenvalues to compute the diffusion matrix at a given temperature. It is supposed here that the eigenvectors are not dependent on temperature.
From Guo and Zhang (2020) (https://10.1016/j.chemgeo.2020.119700).
"""
function Melt_multicomponent_major_Guo2020_SiO2_basaltic()
    data = MeltMulticompDiffusionData(
        Name = "Multicomponent diffusion of major elements in basaltic melt, with SiO2 as dependent variable | Guo and Zhang (2020)",
        Phase = "Melt",
        Species = "TiO2–Al2O3–FeO–MgO–CaO–Na2O–K2O",
        Dependent_Species = "SiO2",
        n = 8NoUnits,
        λD0 = SMatrix{7, 7}(exp(13.752), 0, 0, 0, 0, 0, 0, #=
                         =# 0, exp(14.737), 0, 0, 0, 0, 0, #=
                         =# 0, 0, exp(14.897), 0, 0, 0, 0, #=
                         =# 0, 0, 0, exp(12.375), 0, 0, 0, #=
                         =# 0, 0, 0, 0, exp(15.063), 0, 0, #=
                         =# 0, 0, 0, 0, 0, exp(15.083), 0, #=
                         =# 0, 0, 0, 0, 0, 0, exp(12.18))u"µm^2/s",
        λEa = SMatrix{7, 7}(19_636 * ustrip(Unitful.R), 0, 0, 0, 0, 0, 0, #=
                         =# 0, 20_912 * ustrip(Unitful.R), 0, 0, 0, 0, 0, #=
                         =# 0, 0, 19_987 * ustrip(Unitful.R), 0, 0, 0, 0, #=
                         =# 0, 0, 0, 13_880 * ustrip(Unitful.R), 0, 0, 0, #=
                         =# 0, 0, 0, 0, 18_569 * ustrip(Unitful.R), 0, 0, #=
                         =# 0, 0, 0, 0, 0, 18_279 * ustrip(Unitful.R), 0, #=
                         =# 0, 0, 0, 0, 0, 0, 10_808 * ustrip(Unitful.R))u"J/mol",
        w = SMatrix{7, 7}(-0.76, -0.20, -0.18, -0.02, -0.02, -0.02, -0.02, #=
                       =# -0.18,  0.97, -0.47, -0.15, -0.01, -0.07, -0.10, #=
                       =# -0.51,  0.00,  0.66,  0.86,  0.06, -0.41, -0.36, #=
                       =# -0.17, -0.03,  0.41, -0.14, -0.71, -0.32, -0.15, #=
                       =# -0.22,  0.12,  0.33, -0.33,  0.70,  0.79, -0.08, #=
                       =#  0.17, -0.04, -0.18, -0.12, -0.04, -0.19,  0.91, #=
                       =#  0.13, -0.02, -0.09,  0.32, -0.10,  0.25,  0.06)NoUnits,
        T_range_min = 1260C,
        T_range_max = 1500C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (22.03.25)",
        BibTex_Reference = "
            @article{guo2020multicomponent,
            title={Multicomponent diffusion in a basaltic melt: temperature dependence},
            author={Guo, Chenghuan and Zhang, Youxue},
            journal={Chemical Geology},
            volume={549},
            pages={119700},
            year={2020},
            publisher={Elsevier}
            }
          ",
    )

    return data, info
end
