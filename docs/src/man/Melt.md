## Melt

To obtain the list of implemented diffusion parameters for melt, use:

```@docs
Melt.chemical_diffusion_list
```

For major element diffusion, the multicomponent diffusion matrix for basaltic melt can be calculated using the following function:

```@docs
Melt.Melt_multicomponent_major_Guo2020_SiO2_basaltic
```

For this particular case, to initiate the parameters, use `SetMulticompChemicalDiffusion`:

```julia
    melt_multicomponent = Melt.Melt_multicomponent_major_Guo2020_SiO2_basaltic
    melt_multicomponent = SetMulticompChemicalDiffusion(melt_multicomponent)
```

The function `compute_D` can be used as usual to calculate the diffusion coefficients but will in this case return a [static](https://juliaarrays.github.io/StaticArrays.jl/dev/) matrix of diffusion coefficients for the major elements in the melt:

```julia
    D = compute_D(melt_multicomponent, T=1000C)
7×7 SMatrix{7, 7, Quantity{Float64, 𝐋²·⁰ 𝐓⁻¹·⁰, Unitful.FreeUnits{(m²·⁰, s⁻¹·⁰), 𝐋²·⁰ 𝐓⁻¹·⁰, nothing}}, 49} with indices SOneTo(7)×SOneTo(7):
  1.03519e-13 m²·⁰ s⁻¹·⁰   8.93006e-13 m²·⁰ s⁻¹·⁰  -1.95924e-12 m²·⁰ s⁻¹·⁰  …  -5.84922e-12 m²·⁰ s⁻¹·⁰   8.53056e-12 m²·⁰ s⁻¹·⁰   2.43836e-12 m²·⁰ s⁻¹·⁰
  4.41504e-14 m²·⁰ s⁻¹·⁰   5.07365e-15 m²·⁰ s⁻¹·⁰   1.91304e-13 m²·⁰ s⁻¹·⁰      1.16129e-12 m²·⁰ s⁻¹·⁰  -1.31307e-12 m²·⁰ s⁻¹·⁰  -4.16234e-13 m²·⁰ s⁻¹·⁰
 -2.91173e-13 m²·⁰ s⁻¹·⁰  -1.74973e-13 m²·⁰ s⁻¹·⁰   2.62782e-12 m²·⁰ s⁻¹·⁰      2.77994e-12 m²·⁰ s⁻¹·⁰  -5.25369e-12 m²·⁰ s⁻¹·⁰  -1.49885e-12 m²·⁰ s⁻¹·⁰
 -4.8815e-13 m²·⁰ s⁻¹·⁰    2.53181e-12 m²·⁰ s⁻¹·⁰  -3.84208e-12 m²·⁰ s⁻¹·⁰     -1.54711e-11 m²·⁰ s⁻¹·⁰   2.16066e-11 m²·⁰ s⁻¹·⁰   5.19873e-12 m²·⁰ s⁻¹·⁰
  6.67519e-13 m²·⁰ s⁻¹·⁰  -1.47328e-12 m²·⁰ s⁻¹·⁰  -5.84858e-13 m²·⁰ s⁻¹·⁰      8.55322e-12 m²·⁰ s⁻¹·⁰  -7.44916e-12 m²·⁰ s⁻¹·⁰  -2.12405e-12 m²·⁰ s⁻¹·⁰
 -4.89738e-13 m²·⁰ s⁻¹·⁰   2.07724e-12 m²·⁰ s⁻¹·⁰  -2.8383e-12 m²·⁰ s⁻¹·⁰   …  -1.10541e-11 m²·⁰ s⁻¹·⁰   1.78501e-11 m²·⁰ s⁻¹·⁰   4.10209e-12 m²·⁰ s⁻¹·⁰
 -2.59574e-13 m²·⁰ s⁻¹·⁰   7.36022e-13 m²·⁰ s⁻¹·⁰  -5.65377e-13 m²·⁰ s⁻¹·⁰     -2.34514e-12 m²·⁰ s⁻¹·⁰   3.67054e-12 m²·⁰ s⁻¹·⁰   3.13608e-12 m²·⁰ s⁻¹·⁰
```

Additionally, the function `compute_λ` can be used directly to compute the diagonal matrix of the eigenvalues. This can be useful in the case when the diffusion matrix wants to be diagonalized:

```@docs
GeoParams.MaterialParameters.ChemicalDiffusion.compute_λ
GeoParams.MaterialParameters.ChemicalDiffusion.compute_λ!
```

For trace-, self- or interdiffusion parameters for melt, the following functions are implemented:

```@docs
Melt.Melt_Mg_Sheng1992_basaltic
Melt.Melt_Sc_Holycross2018_rhyolitic_highH2O
Melt.Melt_Sc_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_V_Holycross2018_rhyolitic_highH2O
Melt.Melt_V_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Y_Holycross2018_rhyolitic_highH2O
Melt.Melt_Y_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Zr_Holycross2018_rhyolitic_highH2O
Melt.Melt_Zr_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Hf_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Nb_Holycross2018_rhyolitic_highH2O
Melt.Melt_Nb_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_La_Holycross2018_rhyolitic_highH2O
Melt.Melt_La_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Ce_Holycross2018_rhyolitic_highH2O
Melt.Melt_Ce_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Pr_Holycross2018_rhyolitic_highH2O
Melt.Melt_Pr_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Nd_Holycross2018_rhyolitic_highH2O
Melt.Melt_Nd_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Sm_Holycross2018_rhyolitic_highH2O
Melt.Melt_Sm_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Eu_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Gd_Holycross2018_rhyolitic_highH2O
Melt.Melt_Gd_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Tb_Holycross2018_rhyolitic_highH2O
Melt.Melt_Tb_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Dy_Holycross2018_rhyolitic_highH2O
Melt.Melt_Dy_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Ho_Holycross2018_rhyolitic_highH2O
Melt.Melt_Ho_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Er_Holycross2018_rhyolitic_highH2O
Melt.Melt_Er_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Yb_Holycross2018_rhyolitic_highH2O
Melt.Melt_Yb_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Lu_Holycross2018_rhyolitic_highH2O
Melt.Melt_Lu_Holycross2018_rhyolitic_mediumH2O
Melt.Melt_Th_Holycross2018_rhyolitic_highH2O
Melt.Melt_U_Holycross2018_rhyolitic_highH2O
Melt.Melt_U_Holycross2018_rhyolitic_mediumH2O
```
