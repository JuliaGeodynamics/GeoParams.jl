# Latent heat

# Methods
Latent heat (of crystallisation) can be implemented as a source term (usually numerically not very stable):
```@docs
ConstantLatentHeat
```
Alternatively, you can implement it by modifying the heat capacity, which is often numerically better.
```@docs
Latent_HeatCapacity
```

# Computational routines
To compute with the source term, use this:
```@docs
GeoParams.MaterialParameters.LatentHeat.compute_latent_heat
```

If you modify the heat capacity, you simply use this in your thermal computations:
To compute with the source term, use this:
```@docs
compute_heatcapacity
```
