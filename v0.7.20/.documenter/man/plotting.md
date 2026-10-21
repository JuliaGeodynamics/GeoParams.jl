---
---

# Plotting {#Plotting}

We provide a number of plotting routines. Note that these plotting routines become available only when a backend of [Makie.jl](https://docs.makie.org/stable/) is loaded.

::: warning Missing docstring.

Missing docstring for `GeoParams.PlotStressStrainrate_CreepLaw`. Check Documenter's build log for details.

:::
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PlotHeatCapacity' href='#GeoParams.PlotHeatCapacity'><span class="jlbinding">GeoParams.PlotHeatCapacity</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
PlotHeatCapacity(Cp::AbstractHeatCapacity; kwargs...)
```


Plots heat capacity as a function of temperature for the parameterization `Cp`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/GeoParams.jl#L529-L533" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PlotConductivity' href='#GeoParams.PlotConductivity'><span class="jlbinding">GeoParams.PlotConductivity</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
PlotConductivity(k::AbstractConductivity; kwargs...)
```


Plots thermal conductivity as a function of temperature for the parameterization `k`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/GeoParams.jl#L536-L540" target="_blank" rel="noreferrer">source</a></Badge>

</details>


::: warning Missing docstring.

Missing docstring for `GeoParams.PlotMeltFraction`. Check Documenter's build log for details.

:::
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PlotPhaseDiagram' href='#GeoParams.PlotPhaseDiagram'><span class="jlbinding">GeoParams.PlotPhaseDiagram</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
PlotPhaseDiagram(p::AbstractPhaseDiagramsStruct, fieldname::Symbol; kwargs...)
```


Plots the field `fieldname` of a phase diagram as a function of temperature (x-axis) and pressure (y-axis).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/GeoParams.jl#L550-L554" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PlotDeformationMap' href='#GeoParams.PlotDeformationMap'><span class="jlbinding">GeoParams.PlotDeformationMap</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
PlotDeformationMap(v; kwargs...)
```


Plots a deformation-mechanism map (deformation regime as a function of temperature and stress or strain rate) for the rheology `v`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/GeoParams.jl#L571-L575" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PlotDiffusionCoefArrhenius' href='#GeoParams.PlotDiffusionCoefArrhenius'><span class="jlbinding">GeoParams.PlotDiffusionCoefArrhenius</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
PlotDiffusionCoefArrhenius(x; kwargs...)
```


Arrhenius plot of the diffusion coefficient (`log(D)` versus `10⁴/T`) for one or more `ChemicalDiffusionData` structures `x`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/GeoParams.jl#L595-L599" target="_blank" rel="noreferrer">source</a></Badge>

</details>

