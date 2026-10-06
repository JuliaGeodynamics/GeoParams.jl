# Seismic velocity 

## Methods
Seismic velocity can specified in a number of ways
```@docs
GeoParams.MaterialParameters.SeismicVelocity.ConstantSeismicVelocity
```
In addition, you can use phase diagram lookup tables to compute seismic velocities as a function of pressure and temperature.

# Seismic velocity correction for partial melt

## Methods
Both routines use the velocity reduction of Clark & Lesher (2017). `melt_correction_Takei` models the melt as oblate-spheroidal inclusions of a given aspect ratio (self-consistent model of Dean, 1983, and Phani, 1996); `melt_correction` uses the equilibrium geometry model for the solid skeleton of Takei (1998), parameterized by contiguity.
```@docs
GeoParams.melt_correction_Takei
GeoParams.melt_correction
```

# Seismic S-wave velocity correction for (shallow depth) porosity

## Methods
The porosity follows the empirical porosity-depth relationship of Chen et al. (2020); the pores are treated as fluid-filled inclusions, as in `melt_correction_Takei`.
```@docs
GeoParams.porosity_correction
```

# Correcting the velocities of a phase diagram
```@docs
GeoParams.correct_wavevelocities_phasediagrams
```

# Seismic velocity correction for anelasticity

## Methods
The routine uses the reduction formulation of karato (1993), using the quality factor formulation from Behn et al. (2009)

## Computational routines
To compute a correction of S-wave velocity for anelasticity, use this:
```@docs
GeoParams.anelastic_correction
```
