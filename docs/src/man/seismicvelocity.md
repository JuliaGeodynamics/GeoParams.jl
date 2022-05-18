# Seismic velocity 

# Methods
Seismic velocity can specified in a number of ways
```@docs
GeoParams.MaterialParameters.SeismicVelocity.ConstantSeismicVelocity
```
In addition, you can use phase diagram lookup tables to compute seismic velocities as a function of pressure and temperature.

# Computational routines
To evaluate seismic velocity within a user routine, use this:
```@docs
GeoParams.compute_swave_velocity!
GeoParams.compute_pwave_velocity!
GeoParams.compute_pwave_velocity
GeoParams.compute_swave_velocity
```

# Seismic velocities correction for melt fraction

# Methods
The routine uses the reduction formulation of Clark et al., (2017) and is based on the equilibrium geometry model for the solid skeleton of Takei et al., 1998.

# Computational routines
To compute melt-content based correction for seismic waves velocities, you can use:
```@docs
GeoParams.melt_correction
```