---
---

# Seismic velocity {#Seismic-velocity}

## Methods {#Methods}

Seismic velocity can specified in a number of ways
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.SeismicVelocity.ConstantSeismicVelocity' href='#GeoParams.MaterialParameters.SeismicVelocity.ConstantSeismicVelocity'><span class="jlbinding">GeoParams.MaterialParameters.SeismicVelocity.ConstantSeismicVelocity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstantSeismicVelocity(Vp=8.1 km/s, Vs=4.5km/s)
```


Set a constant seismic P and S-wave velocity:

$$    V_p = cst$$

$$    V_s = cst$$

where $V_p, V_s$ are the P-wave and S-wave velocities [$km/s$].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/SeismicVelocity/SeismicVelocity.jl#L35-L46" target="_blank" rel="noreferrer">source</a></Badge>

</details>


In addition, you can use phase diagram lookup tables to compute seismic velocities as a function of pressure and temperature.

# Seismic velocity correction for partial melt {#Seismic-velocity-correction-for-partial-melt}

## Methods {#Methods-2}

The routine uses the reduction formulation of Clark et al., (2017) and is based on the equilibrium geometry model for the solid skeleton of Takei et al., 1998.

# Seismic S-wave velocity correction for (shallow depth) porosity {#Seismic-S-wave-velocity-correction-for-shallow-depth-porosity}

## Methods {#Methods-3}

This routine is based on the equilibrium geometry model for the solid skeleton of Takei et al. (1998) and the porosity-depth empirical relationship of Chen et al. (2020)

# Seismic velocity correction for anelasticity {#Seismic-velocity-correction-for-anelasticity}

## Methods {#Methods-4}

The routine uses the reduction formulation of karato (1993), using the quality factor formulation from Behn et al. (2009)

## Computational routines {#Computational-routines}

To compute a correction of S-wave velocity for anelasticity, use this:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.SeismicVelocity.anelastic_correction' href='#GeoParams.MaterialParameters.SeismicVelocity.anelastic_correction'><span class="jlbinding">GeoParams.MaterialParameters.SeismicVelocity.anelastic_correction</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
    Vs_anel = anelastic_correction(water::Integer, Vs0, P, T)
```


This routine computes a correction of S-wave velocity for anelasticity

**Input:**
- `water`: water flag, 0 = dry; 1 = dampened; 2 = water saturated
  
- `Vs0`  : S-wave velocitiy of the solid phase (with or without melt correction)
  
- `P`    : pressure given in Pa
  
- `T`    : temperature given in °K
  

**Output:**
- `Vs_anel` : corrected S-wave velocity for anelasticity
  

The routine uses the reduction formulation of Karato (1993), using the quality factor formulation from Behn et al. (2009)

**References:**
- Karato, S. I. (1993). Importance of anelasticity in the interpretation of seismic tomography. Geophysical research letters, 20(15), 1623-1626.
  
- Behn, M. D., Hirth, G., & Elsenbeck II, J. R. (2009). Implications of grain size evolution on the seismic structure of the oceanic upper mantle. Earth and Planetary Science Letters, 282(1-4), 178-189.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/SeismicVelocity/SeismicVelocity.jl#L343-L370" target="_blank" rel="noreferrer">source</a></Badge>

</details>

