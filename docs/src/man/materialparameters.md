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




