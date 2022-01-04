# CreepLaws 

# Implemented creep laws

The following viscous creep laws are implemented:
```@docs
GeoParams.MaterialParameters.CreepLaw.LinearViscous
GeoParams.MaterialParameters.CreepLaw.PowerlawViscous
GeoParams.MaterialParameters.CreepLaw.DislocationCreep
```

# Computational routines for creep laws
Once a creep rheology is defined, we can use the following routines to perform computations within the solvers
```@docs
CreepLawVariables
GeoParams.MaterialParameters.CreepLaw.ComputeCreepLaw_EpsII
GeoParams.MaterialParameters.CreepLaw.ComputeCreepLaw_TauII
```
