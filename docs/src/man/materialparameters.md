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
```

### Computational routines for creep laws
Once a creep rheology is defined, we can use the following routines to perform computations within the solvers
```@docs
GeoParams.MaterialParameters.CreepLaw.CreepLaw_EpsII
GeoParams.MaterialParameters.CreepLaw.CreepLaw_TauII
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

