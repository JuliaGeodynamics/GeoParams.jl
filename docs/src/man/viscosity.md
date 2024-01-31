# Effective viscosity computations

# Methods
A number of viscosity parameterisations are implemented, which can be set with:

```@docs
LinearViscous
LinearMeltViscosity, 
DiffusionCreep, 
DislocationCreep, 
ConstantElasticity, 
DruckerPrager, 
ArrheniusType, 
GrainBoundarySliding, 
PeierlsCreep, 
NonLinearPeierlsCreep
```

# Computational routines
To compute the viscosity given stress of strainrate use:
```@docs
GeoParams.compute_viscosity_εII
GeoParams.compute_viscosity_τII
```