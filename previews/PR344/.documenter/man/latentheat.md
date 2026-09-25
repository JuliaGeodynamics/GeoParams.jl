---
---

# Latent heat {#Latent-heat}

## Methods {#Methods}

Latent heat (of crystallisation) can be implemented as a source term (usually numerically not very stable):
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.LatentHeat.ConstantLatentHeat' href='#GeoParams.MaterialParameters.LatentHeat.ConstantLatentHeat'><span class="jlbinding">GeoParams.MaterialParameters.LatentHeat.ConstantLatentHeat</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstantLatentHeat(Q_L=400kJ/kg)
```


Set a constant latent heat:

$$    Q_L  = cst$$

where $Q_L$ is the latent heat [$J/kg$].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Energy/LatentHeat.jl#L27-L35" target="_blank" rel="noreferrer">source</a></Badge>

</details>


Alternatively, you can implement it by modifying the heat capacity, which is often numerically better.
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.HeatCapacity.Latent_HeatCapacity' href='#GeoParams.MaterialParameters.HeatCapacity.Latent_HeatCapacity'><span class="jlbinding">GeoParams.MaterialParameters.HeatCapacity.Latent_HeatCapacity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Latent_HeatCapacity(Cp=ConstantHeatCapacity(), Q_L=400kJ/kg)
```


This takes the effects of latent heat into account by modifying the heat capacity in the temperature equation:

$$\rho C_p^{\textrm{new}} \frac{\partial T}{\partial t}  = \frac{\partial }{\partial x_i} \left( k \frac{\partial T}{\partial x_i} \right)  + H_s$$

with

$$C_p^{\textrm{new}}  = C_p + \frac{\partial \phi}{\partial T} Q_L$$

where $Q_L$ is the latent heat [$J/kg$], and $\frac{\partial \phi}{\partial T}$ is the derivative of the melt fraction with respect to temperature


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Energy/HeatCapacity.jl#L115-L130" target="_blank" rel="noreferrer">source</a></Badge>

</details>

