# Density 
# Methods
The density equation of state can be specified in a number of ways
```@docs
GeoParams.MaterialParameters.Density.ConstantDensity
GeoParams.MaterialParameters.Density.PT_Density
```
# Computational routines
To evaluate density within a user routine, use this:
```@docs
GeoParams.MaterialParameters.Density.ComputeDensity
GeoParams.MaterialParameters.Density.ComputeDensity!
```
Note that density values are usually not used in itself in the governing PDE's, but usually in combination with other parameters, such as $\rho g$ or $\rho c_p$. the non-dimensional value of $\rho$ may thus have very large or small values, but multiplied with the other values one often obtains numbers that are closer to one.
