# MaterialParameters

Material properties for a given phase can be set with `SetMaterialParams`, whereas all properties are 
store in the `MaterialParams` structure.

```@docs
SetMaterialParams
MaterialParams
```

# CreepLaws 



### Implemented creep laws

The following viscous creep laws are implemented:
```@docs
GeoParams.MaterialParameters.CreepLaw.LinearViscous
GeoParams.MaterialParameters.CreepLaw.PowerlawViscous
GeoParams.MaterialParameters.CreepLaw.DislocationCreep
```

### Computational routines for creep laws
Once a creep rheology is defined, we can use the following routines to perform computations within the solvers
```@docs
CreepLawVariables
GeoParams.MaterialParameters.CreepLaw.ComputeCreepLaw_EpsII
GeoParams.MaterialParameters.CreepLaw.ComputeCreepLaw_TauII
```

# Density 
The density equation of state can be specified in a number of ways
```@docs
GeoParams.MaterialParameters.Density.ConstantDensity
GeoParams.MaterialParameters.Density.PT_Density
```

To evaluate density within a user routine, use this:
```@docs
GeoParams.MaterialParameters.Density.ComputeDensity
```
Note that density values are usually not used in itself in the governing PDE's, but usually in combination with other parameters, such as $\rho g$ or $\rho c_p$. the non-dimensional value of $\rho$ may thus have very large or small values, but multiplied with the other values one often obtains numbers that are closer to one.

# Gravitational acceleration 
Gravitational acceleration is defined as 
```@docs
GeoParams.MaterialParameters.GravitationalAcceleration.ConstantGravity
```
To compute, use this:
```@docs
GeoParams.MaterialParameters.GravitationalAcceleration.ComputeGravity
```

# Heat capacity
Heat capacity is defined as 
```@docs
GeoParams.MaterialParameters.HeatCapacity.ConstantHeatCapacity
GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Whittacker
```
To compute, use this:
```@docs
GeoParams.MaterialParameters.HeatCapacity.ComputeHeatCapacity
```


# Conductivity
Thermal conductivity is defined as 
```@docs
GeoParams.MaterialParameters.Conductivity.ConstantConductivity
GeoParams.MaterialParameters.Conductivity.T_Conductivity_Whittacker
GeoParams.MaterialParameters.Conductivity.TP_Conductivity
GeoParams.MaterialParameters.Conductivity.Set_TP_Conductivity
```
To compute, use this:
```@docs
GeoParams.MaterialParameters.Conductivity.ComputeConductivity
```


# Latent heat
Latent heat (of crystallisation) is defined as 
```@docs
GeoParams.MaterialParameters.LatentHeat.ConstantLatentHeat
```
To compute, use this:
```@docs
GeoParams.MaterialParameters.LatentHeat.ComputeLatentHeat
```
# Radioactive heat
Radioactive heat sources are defined as 
```@docs
GeoParams.MaterialParameters.RadioactiveHeat.ConstantRadioactiveHeat
```
To compute, use this:
```@docs
GeoParams.MaterialParameters.RadioactiveHeat.ComputeRadioactiveHeat
```

# Shear heating 
Heat caused by non-recoverable deformation can be specified in 
```@docs
GeoParams.MaterialParameters.Shearheating.ConstantShearheating
```
To compute, use this:
```@docs
GeoParams.MaterialParameters.Shearheating.ComputeShearheating
GeoParams.MaterialParameters.Shearheating.ComputeShearheating!
```