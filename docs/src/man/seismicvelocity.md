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
