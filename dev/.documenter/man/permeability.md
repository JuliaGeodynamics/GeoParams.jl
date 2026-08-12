---
---

# Permeability Parameterizations {#Permeability-Parameterizations}

## Methods {#Methods}

A number of permeability parameterisations are implemented, which can be set with:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Permeability.ConstantPermeability' href='#GeoParams.MaterialParameters.Permeability.ConstantPermeability'><span class="jlbinding">GeoParams.MaterialParameters.Permeability.ConstantPermeability</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstantPermeability(k = 1e-12m^2)
```


Defines a constant permeability value for a given material.

**Arguments**
- `k`: The permeability value in square meters (m^2). Default is `1e-12 m^2`.
  

**Example**

```julia
rheology = SetMaterialParams(;
                      Phase=1,
                      CreepLaws=(PowerlawViscous(), LinearViscous(; η=1e21Pa * s)),
                      Gravity=ConstantGravity(; g=9.81.0m / s^2),
                      Density= MeltDependent_Density(),
                      Permeability = ConstantPermeability(; k=1e-12m^2),
                      )
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Permeability/Permeability.jl#L37-L55" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Permeability.HazenPermeability' href='#GeoParams.MaterialParameters.Permeability.HazenPermeability'><span class="jlbinding">GeoParams.MaterialParameters.Permeability.HazenPermeability</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
HazenPermeability(C = 1.0 * NoUnits, D10 = 1e-4m)
```


Defines the Hazen permeability equation for a given material.

$$    k = C \cdot D_{10}^2$$

**Arguments**
- `C`: The Hazen constant. Default is `1.0`.
  
- `D10`: The effective grain size. Default is `1e-4 m`.
  

**Example**

```julia
rheology = SetMaterialParams(;
                      Phase=1,
                      CreepLaws=(PowerlawViscous(), LinearViscous(; η=1e21Pa * s)),
                      Gravity=ConstantGravity(; g=9.81.0m / s^2),
                      Density= MeltDependent_Density(),
                      Permeability = HazenPermeability(; C=1.0, D10=1e-4m),
                      )
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Permeability/Permeability.jl#L86-L108" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Permeability.PowerLawPermeability' href='#GeoParams.MaterialParameters.Permeability.PowerLawPermeability'><span class="jlbinding">GeoParams.MaterialParameters.Permeability.PowerLawPermeability</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PowerLawPermeability(c = 1.0, k0 = 1e-12m^2, ϕ = 1e-2, n = 3)
```


Defines the power-law permeability equation for a given material.

$$    c * k_0 * \phi^{n}$$

**Arguments**
- `c`: The power-law constant. Default is `1.0`.
  
- `k0`: The reference permeability. Default is `1e-12 m^2`.
  
- `ϕ`: The reference porosity. Default is `1e-2`.
  
- `n`: The exponent. Default is `3`.
  

**Example**

```julia
rheology = SetMaterialParams(;
                      Phase=1,
                      CreepLaws=(PowerlawViscous(), LinearViscous(; η=1e21Pa * s)),
                      Gravity=ConstantGravity(; g=9.81.0m / s^2),
                      Density= MeltDependent_Density(),
                      Permeability = PowerLawPermeability(; c=1.0, k0=1e-12m^2, ϕ=1e-2, n=3),
                      )
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Permeability/Permeability.jl#L139-L163" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Permeability.CarmanKozenyPermeability' href='#GeoParams.MaterialParameters.Permeability.CarmanKozenyPermeability'><span class="jlbinding">GeoParams.MaterialParameters.Permeability.CarmanKozenyPermeability</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
CarmanKozenyPermeability(c = 1.0m^2, ϕ0 = 0.01, n = 3)
```


Defines the Carman-Kozeny permeability equation for a given material.

$$    k = c \left(\frac{\phi}{\phi_0}\right)^n$$

**Arguments**
- `c`: The Carman-Kozeny constant. Default is `1.0 m^2`.
  
- `ϕ0`: The reference porosity. Default is `0.01`.
  
- `n`: The exponent. Default is `3`.
  

**Example**

```julia
rheology = SetMaterialParams(;
                      Phase=1,
                      CreepLaws=(PowerlawViscous(), LinearViscous(; η=1e21Pa * s)),
                      Gravity=ConstantGravity(; g=9.81.0m / s^2),
                      Density= MeltDependent_Density(),
                      Permeability = CarmanKozenyPermeability(; c=1.0m^2, ϕ0=0.01, n=3),
                      )
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Permeability/Permeability.jl#L195-L218" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Computational routines {#Computational-routines}

To compute the permeability, use:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Permeability.compute_permeability!' href='#GeoParams.MaterialParameters.Permeability.compute_permeability!'><span class="jlbinding">GeoParams.MaterialParameters.Permeability.compute_permeability!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_permeability!(k::AbstractArray{_T, N}, MatParam::NTuple{K,AbstractMaterialParamsStruct}, PhaseRatios::AbstractArray{_T, M}, P=nothing, T=nothing)
```


In-place computation of permeability `k` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays. This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Permeability/Permeability.jl#L260-L265" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Permeability.compute_permeability' href='#GeoParams.MaterialParameters.Permeability.compute_permeability'><span class="jlbinding">GeoParams.MaterialParameters.Permeability.compute_permeability</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_permeability(k::AbstractPermeability, args)
```


Computation of permeability `k` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Permeability/Permeability.jl#L250-L254" target="_blank" rel="noreferrer">source</a></Badge>

</details>

