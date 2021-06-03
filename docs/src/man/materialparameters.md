# MaterialParameters

The material parameters structure holds all data for a given phase.

```@docs
MaterialParams
```

# CreepLaws 



### Implemented creep laws

The following viscous creep laws are implemented:
```@docs
GeoParams.MaterialParameters.LinearViscous
```

Once a creep rheology is defined, we can use the following routines to perform computations within the solvers
```@docs
GeoParams.MaterialParameters.CreepLaw_EpsII
GeoParams.MaterialParameters.CreepLaw_TauII
```

### Computing 

# Density 
The density equation of state can be specified in a number of ways
