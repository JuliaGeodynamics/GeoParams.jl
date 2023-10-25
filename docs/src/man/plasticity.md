# Plasticity 

Plasticity is a non-linear rheology that is activated once stresses exceed a certain yield criteria.
# Implemented laws
The following plastic law are implemented:
```@docs
GeoParams.MaterialParameters.Plasticity.DruckerPrager
GeoParams.MaterialParameters.Plasticity.DruckerPrager_regularised
```

# Computational routines 
Usually, plasticity should be defined as part of a `CompositeRheology` structure and calculations can be done as with all other rheology computations by using `compute_τII`.
