# Chemical Diffusion

Some routines and experimental data for chemical diffusion in minerals and phases are implemented in this module. Contributions are welcome to extend the database.

## Methods
The following diffusion parameters are implemented:

- Rutile:
```@docs
Rutile.Rt_Hf_Cherniak2007_para_c
Rutile.Rt_Hf_Cherniak2007_perp_c
Rutile.Rt_Zr_Cherniak2007_para_c
```

- Garnet:
```@docs
Garnet.Grt_Fe_Chakraborty1992
Garnet.Grt_Mg_Chakraborty1992
Garnet.Grt_Mn_Chakraborty1992
```

## Computational routines
To compute the diffusion coefficient from the parameters, use this:
```@docs
GeoParams.MaterialParameters.ChemicalDiffusion.compute_D
GeoParams.MaterialParameters.ChemicalDiffusion.compute_D!
```

## List of diffusion parameters

To know the available diffusion parameters for each phase, you can use `X.chemical_diffusion_list()`, where `X` is the phase of interest.
For instance, for garnet, one can use `Garnet.chemical_diffusion_list()`.

```@docs
Rutile.chemical_diffusion_list
Garnet.chemical_diffusion_list
```
