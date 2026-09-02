---
---

# Heat capacity {#Heat-capacity}

## Methods {#Methods}

Heat capacity is defined as 
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.HeatCapacity.ConstantHeatCapacity' href='#GeoParams.MaterialParameters.HeatCapacity.ConstantHeatCapacity'><span class="jlbinding">GeoParams.MaterialParameters.HeatCapacity.ConstantHeatCapacity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstantHeatCapacity(Cp=1050J/kg/K)
```


Set a constant heat capacity:

$$    Cp  = cst$$

where $Cp$ is the thermal heat capacity [$J/kg/K$].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Energy/HeatCapacity.jl#L29-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Whittington' href='#GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Whittington'><span class="jlbinding">GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Whittington</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
T_HeatCapacity_Whittington()
```


Sets a temperature-dependent heat capacity following the parameterization of Whittington et al. (2009), Nature:

$$    Cp = (a + b T - c/T^2)/m$$

where $Cp$ is the heat capacity [$J/kg/K$], and $a,b,c$ are parameters that dependent on the temperature $T$ [$K$]:
- a = 199.50 J/mol/K    if T&lt;= 846 K
  
- a = 229.32 J/mol/K    if T&gt; 846 K
  
- b = 0.0857J/mol/K^2   if T&lt;= 846 K
  
- b = 0.0323J/mol/K^2   if T&gt; 846 K
  
- c = 5e6J/mol*K        if T&lt;= 846 K
  
- c = 47.9e-6J/mol*K    if T&gt; 846 K
  
- molmass =   0.22178kg/mol
  

Note that this is slightly different than the equation in the manuscript, as Cp is in J/kg/K (rather than $J/mol/K$ as in eq.3/4 of the paper)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Energy/HeatCapacity.jl#L60-L78" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Computational routines {#Computational-routines}

To compute, use this:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.HeatCapacity.compute_heatcapacity' href='#GeoParams.MaterialParameters.HeatCapacity.compute_heatcapacity'><span class="jlbinding">GeoParams.MaterialParameters.HeatCapacity.compute_heatcapacity</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_heatcapacity(a::Vector_HeatCapacity; index::Int64, kwargs...)
```


Pointwise calculation of heat capacity from a vector where `index` is the index of the point


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Energy/HeatCapacity.jl#L181-L185" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_heatcapacity(P,T, s::AbstractPhaseDiagramsStruct)
```


Interpolates heat capacity as a function of `T,P` from a lookup table


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Energy/HeatCapacity.jl#L226-L229" target="_blank" rel="noreferrer">source</a></Badge>



```julia
Cp = compute_heatcapacity(s:<AbstractHeatCapacity, P, T)
```


Returns the heat capacity `Cp` at any temperature `T` and pressure `P` using any of the heat capacity laws implemented.

Currently available:
- ConstantHeatCapacity
  
- T_HeatCapacity_Whittington
  
- Latent_HeatCapacity
  

**Example**

Using dimensional units

```julia
julia> T = 250.0:100:1250
julia> Cp2 = T_HeatCapacity_Whittington()
julia> Cp = similar(T)
julia> args = (; T=T)
julia> Cp =compute_heatcapacity!(Cp, Cp2, args)
11-element Vector nitful.Quantity{Float64, 𝐋² 𝚯⁻¹ 𝐓⁻², Unitful.FreeUnits{(kg⁻¹, J, K⁻¹), 𝐋² 𝚯⁻¹ 𝐓⁻², nothing}}}:
635.4269997294616 J kg⁻¹ K⁻¹
850.7470171764261 J kg⁻¹ K⁻¹
962.0959598489883 J kg⁻¹ K⁻¹
1037.542043377064 J kg⁻¹ K⁻¹
1097.351792196648 J kg⁻¹ K⁻¹
1149.274556367170 J kg⁻¹ K⁻¹
1157.791505094840 J kg⁻¹ K⁻¹
1172.355487419726 J kg⁻¹ K⁻¹
1186.919469744596 J kg⁻¹ K⁻¹
1201.483452069455 J kg⁻¹ K⁻¹
1216.0474343943067 J kg⁻¹ K⁻¹
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Energy/HeatCapacity.jl#L236-L269" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.HeatCapacity.compute_heatcapacity!' href='#GeoParams.MaterialParameters.HeatCapacity.compute_heatcapacity!'><span class="jlbinding">GeoParams.MaterialParameters.HeatCapacity.compute_heatcapacity!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_heatcapacity!(Cp::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat})
```


In-place computation of heat capacity `Cp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays. This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Energy/HeatCapacity.jl#L291-L296" target="_blank" rel="noreferrer">source</a></Badge>

</details>

