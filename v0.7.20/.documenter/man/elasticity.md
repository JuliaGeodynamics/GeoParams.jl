---
---

# Elasticity {#Elasticity}

Elasticity is, in geodynamics, often used in combination with viscous and plastic rheology. We provide the following elastic constitutive relationships:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.ConstantElasticity' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.ConstantElasticity'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.ConstantElasticity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstantElasticity(G, ν, K, E)
```


Structure that holds parameters for constant, isotropic, linear elasticity.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/Elasticity/Elasticity.jl#L34-L38" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.SetConstantElasticity' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.SetConstantElasticity'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.SetConstantElasticity</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
SetConstantElasticity(; G=nothing, ν=nothing, E=nothing, Kb=nothing)
```


This allows setting elastic parameters by specifying 2 out of the 4 elastic parameters `G` (Elastic shear modulus), `ν` (Poisson's ratio), `E` (Young's modulus), or `Kb` (bulk modulus).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/Elasticity/Elasticity.jl#L50-L54" target="_blank" rel="noreferrer">source</a></Badge>

</details>

