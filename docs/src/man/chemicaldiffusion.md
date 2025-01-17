# Chemical Diffusion

Some routines and experimental data for chemical diffusion in minerals and phases are implemented in this module. Contributions are welcome to extend the database.

## Methods
The following diffusion parameters are implemented:

- Hf in Rutile:
```@docs
GeoParams.Rutile.Rt_Hf_Cherniak2007_Ξc
GeoParams.Rutile.Rt_Hf_Cherniak2007_⊥c
```

## Computational routines
To compute the diffusion coefficient from the parameters, use this:
```@docs
GeoParams.MaterialParameters.ChemicalDiffusion.compute_D
GeoParams.MaterialParameters.ChemicalDiffusion.compute_D!
```
