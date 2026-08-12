---
---

# Gravitational acceleration {#Gravitational-acceleration}

## Methods {#Methods}

Gravitational acceleration is defined as 
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.GravitationalAcceleration.ConstantGravity' href='#GeoParams.MaterialParameters.GravitationalAcceleration.ConstantGravity'><span class="jlbinding">GeoParams.MaterialParameters.GravitationalAcceleration.ConstantGravity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GravityConstant(g=9.81m/s^2)
```


Set a constant value for the gravitational acceleration:

$$    g  = 9.81 m s^{-2}$$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/GravitationalAcceleration/GravitationalAcceleration.jl#L19-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.GravitationalAcceleration.DippingGravity' href='#GeoParams.MaterialParameters.GravitationalAcceleration.DippingGravity'><span class="jlbinding">GeoParams.MaterialParameters.GravitationalAcceleration.DippingGravity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DippingGravity(g=9.81m/s^2)
```


Set a constant value for the gravitational acceleration with dip and strike angles:

$$    g  = R_z  R_y  9.81 m s^{-2} $$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/GravitationalAcceleration/GravitationalAcceleration.jl#L50-L57" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Computational routines {#Computational-routines}

To compute, use this:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.GravitationalAcceleration.compute_gravity' href='#GeoParams.MaterialParameters.GravitationalAcceleration.compute_gravity'><span class="jlbinding">GeoParams.MaterialParameters.GravitationalAcceleration.compute_gravity</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



compute_gravity(s:&lt;AbstractGravity)

Returns the gravitational acceleration 


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/GravitationalAcceleration/GravitationalAcceleration.jl#L120-L125" target="_blank" rel="noreferrer">source</a></Badge>

</details>

