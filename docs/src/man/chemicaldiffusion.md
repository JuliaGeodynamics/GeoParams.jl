# Chemical Diffusion

Some routines and experimental data for chemical diffusion in minerals and phases are implemented in this module. Contributions are welcome to extend the database.

## Methods
The following diffusion parameters are implemented:

- Rutile:
```@docs
GeoParams.Rutile.Rt_Hf_Cherniak2007_para_c
GeoParams.Rutile.Rt_Hf_Cherniak2007_perp_c
```

- Garnet:
```@docs
GeoParams.Garnet.Grt_Fe_Chakraborty1992
GeoParams.Garnet.Grt_Mg_Chakraborty1992
GeoParams.Garnet.Grt_Mn_Chakraborty1992
```

## Computational routines
To compute the diffusion coefficient from the parameters, use this:
```@docs
GeoParams.MaterialParameters.ChemicalDiffusion.compute_D
GeoParams.MaterialParameters.ChemicalDiffusion.compute_D!
```
