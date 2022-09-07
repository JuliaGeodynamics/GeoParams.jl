# CreepLaws 

# Implemented creep laws

The following viscous creep laws are implemented:
```@docs
GeoParams.LinearViscous
GeoParams.PowerlawViscous
GeoParams.DislocationCreep
GeoParams.SetDislocationCreep
GeoParams.DiffusionCreep
GeoParams.SetDiffusionCreep
```

# Computational routines for creep laws
Once a creep rheology is defined, we can use the following routines to perform computations within the solvers
```@docs
CreepLawVariables
GeoParams.compute_εII
GeoParams.compute_εII!
GeoParams.compute_τII
GeoParams.compute_τII!
GeoParams.CorrectionFactor
GeoParams.RemoveTensorCorrection
```
