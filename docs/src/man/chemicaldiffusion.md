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
Melt.Sc_Melt_Holycross2018_rhyolitic_highH2O
Melt.Sc_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.V_Melt_Holycross2018_rhyolitic_highH2O
Melt.V_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Y_Melt_Holycross2018_rhyolitic_highH2O
Melt.Y_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Zr_Melt_Holycross2018_rhyolitic_highH2O
Melt.Zr_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Hf_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Nb_Melt_Holycross2018_rhyolitic_highH2O
Melt.Nb_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.La_Melt_Holycross2018_rhyolitic_highH2O
Melt.La_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Ce_Melt_Holycross2018_rhyolitic_highH2O
Melt.Ce_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Pr_Melt_Holycross2018_rhyolitic_highH2O
Melt.Pr_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Nd_Melt_Holycross2018_rhyolitic_highH2O
Melt.Nd_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Sm_Melt_Holycross2018_rhyolitic_highH2O
Melt.Sm_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Eu_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Gd_Melt_Holycross2018_rhyolitic_highH2O
Melt.Gd_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Tb_Melt_Holycross2018_rhyolitic_highH2O
Melt.Tb_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Dy_Melt_Holycross2018_rhyolitic_highH2O
Melt.Dy_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Ho_Melt_Holycross2018_rhyolitic_highH2O
Melt.Ho_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Er_Melt_Holycross2018_rhyolitic_highH2O
Melt.Er_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Yb_Melt_Holycross2018_rhyolitic_highH2O
Melt.Yb_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Lu_Melt_Holycross2018_rhyolitic_highH2O
Melt.Lu_Melt_Holycross2018_rhyolitic_mediumH2O
Melt.Th_Melt_Holycross2018_rhyolitic_highH2O
Melt.U_Melt_Holycross2018_rhyolitic_highH2O
Melt.U_Melt_Holycross2018_rhyolitic_mediumH2O
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
