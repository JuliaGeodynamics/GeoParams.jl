# Elasticity 

Elasticity is, in geodynamics, often used in combination with viscous and plastic rheology.

# Implemented laws
We provide the following elastic constitutive relationships:
```@docs
GeoParams.ConstantElasticity
GeoParams.SetConstantElasticity
```

# Computational routines 
We can compute the elastic strainrate with:
```@docs
GeoParams.compute_εII
GeoParams.compute_εII!
```
