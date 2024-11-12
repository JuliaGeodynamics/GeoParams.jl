# Latent heat

## Methods
Latent heat (of crystallisation) can be implemented as a source term (usually numerically not very stable):
```@docs
ConstantLatentHeat
```
Alternatively, you can implement it by modifying the heat capacity, which is often numerically better.
```@docs
Latent_HeatCapacity
```
