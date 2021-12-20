# Seismic velocity 

# Methods
Seismic velocity can specified in a number of ways
```@docs
GeoParams.ConstantSeismicVelocity
```
In addition, you can use phase diagram lookup tables to compute seismic velocities as a function of pressure and temperature.

# Computational routines
To evaluate seismic velocity within a user routine, use this:
```@docs
GeoParams.ComputeSwaveVelocity!
GeoParams.ComputePwaveVelocity!
GeoParams.ComputePwaveVelocity
GeoParams.ComputeSwaveVelocity
```
