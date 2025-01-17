# Seismic velocity 

## Methods
Seismic velocity can specified in a number of ways
```@docs
GeoParams.MaterialParameters.SeismicVelocity.ConstantSeismicVelocity
```
In addition, you can use phase diagram lookup tables to compute seismic velocities as a function of pressure and temperature.

# Seismic velocity correction for partial melt

## Methods
The routine uses the reduction formulation of Clark et al., (2017) and is based on the equilibrium geometry model for the solid skeleton of Takei et al., 1998.

# Seismic S-wave velocity correction for (shallow depth) porosity

## Methods
This routine is based on the equilibrium geometry model for the solid skeleton of Takei et al. (1998) and the porosity-depth empirical relationship of Chen et al. (2020)

# Seismic velocity correction for anelasticity

## Methods
The routine uses the reduction formulation of karato (1993), using the quality factor formulation from Behn et al. (2009)

## Computational routines
To compute a correction of S-wave velocity for anelasticity, use this:
```@docs
GeoParams.anelastic_correction
```
