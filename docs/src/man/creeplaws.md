# CreepLaws

The following viscous creep laws are implemented:
```@docs
GeoParams.LinearViscous
GeoParams.PowerlawViscous
GeoParams.DislocationCreep
GeoParams.MaterialParameters.ConstitutiveRelationships.ViscosityPartialMelt_Costa_etal_2009
GeoParams.MaterialParameters.ConstitutiveRelationships.LinearMeltViscosity
GeoParams.MaterialParameters.ConstitutiveRelationships.GiordanoMeltViscosity
GeoParams.MaterialParameters.ConstitutiveRelationships.DiffusionCreep
GeoParams.MaterialParameters.ConstitutiveRelationships.SetDiffusionCreep
```

## Computational routines for creep laws

Once a creep rheology is defined, we can use the following routines to perform computations within the solvers
```@docs
GeoParams.compute_εII
GeoParams.compute_εII!
GeoParams.compute_τII
GeoParams.compute_τII!
```
