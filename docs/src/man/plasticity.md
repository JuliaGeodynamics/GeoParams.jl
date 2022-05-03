# Plasticity 

Plasticity is a non;linear rheology bthat is actrivated once stresses exceed a certain yield criteria.
# Implemented  laws
The following plastic law are implemented:
```@docs
GeoParams.MaterialParameters.Plasticity.DruckerPrager
```

# Computational routines 
Once a plastic rheology is defined, we can use the following routines to perform computations within the solvers
```@docs
GeoParams.MaterialParameters.Plasticity.compute_yieldfunction
GeoParams.MaterialParameters.Plasticity.compute_yieldfunction!
```
