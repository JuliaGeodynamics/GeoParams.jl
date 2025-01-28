# Chemical Diffusion

Some routines and experimental data for chemical diffusion in minerals and phases are implemented in this module. Contributions are welcome to extend the database.

## Methods
The following diffusion parameters are implemented:

- Rutile:
```@docs
Rutile.Rt_Zr_Cherniak2007_para_c
Rutile.Rt_Hf_Cherniak2007_para_c
Rutile.Rt_Hf_Cherniak2007_perp_c
```

- Garnet:
```@docs
Garnet.Grt_Mg_Chakraborty1992
Garnet.Grt_Mn_Chakraborty1992
Garnet.Grt_Fe_Chakraborty1992
Garnet.Grt_REE_Bloch2020_slow
Garnet.Grt_REE_Bloch2020_fast
```

- Melt:
```@docs
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
Melt.Melt_Sm_Holycross2018_rhyolitic_mediumH2O
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

## Computational routines
To compute the diffusion coefficient from the parameters, use this:
```@docs
GeoParams.MaterialParameters.ChemicalDiffusion.compute_D
GeoParams.MaterialParameters.ChemicalDiffusion.compute_D!
```

## List of diffusion parameters

To know the available diffusion parameters for each phase, you can use `X.chemical_diffusion_list()`, where `X` is the phase of interest.
For instance, for garnet, one can use `Garnet.chemical_diffusion_list()`. In addition, this function also takes an argument to search for a specific term, i.e. an element ("La") or an author.

```@docs
Rutile.chemical_diffusion_list
Garnet.chemical_diffusion_list
Melt.chemical_diffusion_list
```
