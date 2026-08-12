---
---

# Radioactive heat {#Radioactive-heat}

## Methods {#Methods}

Radioactive heat sources are defined as 
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.RadioactiveHeat.ConstantRadioactiveHeat' href='#GeoParams.MaterialParameters.RadioactiveHeat.ConstantRadioactiveHeat'><span class="jlbinding">GeoParams.MaterialParameters.RadioactiveHeat.ConstantRadioactiveHeat</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstantRadioactiveHeat(H_r=1e-6Watt/m^3)
```


Set a constant radioactive heat:

$$    H_r  = cst$$

where $H_r$ is the radioactive heat source [$Watt/m^3$].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Energy/RadioactiveHeat.jl#L23-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.RadioactiveHeat.ExpDepthDependentRadioactiveHeat' href='#GeoParams.MaterialParameters.RadioactiveHeat.ExpDepthDependentRadioactiveHeat'><span class="jlbinding">GeoParams.MaterialParameters.RadioactiveHeat.ExpDepthDependentRadioactiveHeat</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ExpDepthDependent(H_0=1e-6Watt/m^3, h_r=10e3m, z_0=0m)
```


Sets an exponential depth-dependent radioactive

$$    H_r  = H_0 \exp \left( -\frac{z - z_0}{h_r} \right)$$

where $H_0$ is the radioactive heat source [$Watt/m^3$] at $z=z_0$ which decays with depth over a characteristic distance $h_r$.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Energy/RadioactiveHeat.jl#L63-L71" target="_blank" rel="noreferrer">source</a></Badge>

</details>

