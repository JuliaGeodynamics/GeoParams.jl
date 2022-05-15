# Elasticity 

Elasticity is, in geodynamics, often used in combination with viscous and plastic rheology.

# Implemented laws
We provide the following elastic constitutive relationships:
```@docs
GeoParams.MaterialParameters.Elasticity.ConstantElasticity
```

# Computational routines 
We can compute the elastic strainrate with:
```@docs
GeoParams.MaterialParameters.Elasticity.compute_elastic_shear_strainrate
GeoParams.MaterialParameters.Elasticity.compute_elastic_shear_strainrate!
```
