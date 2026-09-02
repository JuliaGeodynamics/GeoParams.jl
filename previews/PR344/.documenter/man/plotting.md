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
 fig,ax,T,Cp_vec = PlotHeatCapacity(Cp::AbstractHeatCapacity; T=nothing, plt=nothing, lbl=nothing)
```


Creates a plot of temperature `T` vs. heat capacity, as specified in Cp (which can be temperature-dependent).

**Optional parameters**
- T: temperature range
  
- plt: a previously generated plotting object
  
- lbl: label of the curve
  

**Example**

```
julia> Cp = T_HeatCapacity_Whittacker()
julia> fig,ax,T,Cp_vec = PlotHeatCapacity(Cp)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Plotting/Plotting.jl#L536-L551" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PlotConductivity' href='#GeoParams.PlotConductivity'><span class="jlbinding">GeoParams.PlotConductivity</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fig, ax, T, Cond = PlotConductivity(k::AbstractConductivity; T=nothing, P=nothing, fig=nothing, ax=nothing, lbl=nothing, filename=nothing)
```


Creates a plot of temperature `T` vs. thermal conductivity, as specified in `k` (which can be temperature-dependent). Note: if you want to create plots you need to install and load a `Makie.jl` backend (e.g. `GLMakie.jl` or, for headless use, `CairoMakie.jl`).

**Optional parameters**
- `T`: temperature range
  
- `P`: pressure (scalar or array matching `T`)
  
- `fig`/`ax`: a previously generated figure/axis to plot into
  
- `lbl`: label of the curve
  
- `filename`: if provided, the figure is saved to this file instead of being displayed
  

**Example**

```
julia> using CairoMakie, GeoParams
julia> k = T_Conductivity_Whittington()
julia> fig, ax, T, Cond = PlotConductivity(k)
julia> save("Tdependent_conductivity.png", fig)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Plotting/Plotting.jl#L641-L661" target="_blank" rel="noreferrer">source</a></Badge>

</details>


::: warning Missing docstring.

Missing docstring for `GeoParams.PlotMeltFraction`. Check Documenter's build log for details.

:::
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PlotPhaseDiagram' href='#GeoParams.PlotPhaseDiagram'><span class="jlbinding">GeoParams.PlotPhaseDiagram</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fig, ax, Z, Tvec_K, Pvec_Pa = PlotPhaseDiagram(p::AbstractPhaseDiagramsStruct, fieldname::Symbol; Tvec=nothing, Pvec=nothing, fig=nothing, filename=nothing)
```


Plots a phase diagram field `fieldname` as a function of `T` (x-axis, in K) and `P` (y-axis, in Pa). We either use the default ranges of the diagram, or you can specify the temperature and pressure ranges (while specifying units). The return arguments are the figure/axis, the gridded data `Z` and the temperature/pressure vectors used. Note: if you want to create plots you need to install and load a `Makie.jl` backend (e.g. `GLMakie.jl` or, for headless use, `CairoMakie.jl`).

**Example**

```julia
julia> using CairoMakie, GeoParams
julia> PD_data = PerpleX_LaMEM_Diagram("./test/test_data/Peridotite.in")
julia> PlotPhaseDiagram(PD_data, :meltFrac, Tvec=(100:1:1400).*C, Pvec=(.1:.1:30).*kbar)
```


You can also use the default pressure/temperature ranges in the diagrams:

```julia
julia> PlotPhaseDiagram(PD_data, :Rho)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Plotting/Plotting.jl#L770-L789" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PlotDeformationMap' href='#GeoParams.PlotDeformationMap'><span class="jlbinding">GeoParams.PlotDeformationMap</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fig = PlotDeformationMap(v;    args=(P=0.0, d=1e-3, f=1.0),
                            σ = (1e-2, 1e8),                # in MPa
                            T = (10, 1000),                 # in C
                            ε = (1e-22, 1e-8),              # in 1/s
                            n = 400,                        # number of points
                            viscosity_cutoff = (0, 25),     # viscosity cutoff in log10()
                            rotate_axes = false,            # flip x & y axes
                            strainrate = true,              # strainrate (otherwise stress)
                            log_stress = true,              # log10-scale stress: axis on strainrate map, τII color on stress map; false = linear MPa
                            viscosity = false,              # plot viscosity instead of strainrate/stress
                            boundaries = true,              # plot deformation boundaries
                            levels = 20,                    # number of contour levels
                            colormap=:viridis,              # colormap
                            filename=nothing,               # if you want to save this to file
                            fontsize=40,                    # fontsize of labels
                            res=(1200, 900))                # resolution in pixels
```


Creates a deformation mechanism map (T/εII vs. stress/viscosity or T/τII vs. strainrate/viscosity) for given (composite) rheology `v`

**Example**

```julia
julia> import GeoParams.Diffusion, GeoParams.Dislocation
julia> v1 = SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006);
julia> v2 = SetDislocationCreep(Dislocation.dry_anorthite_Rybacki_2006);
julia> v = CompositeRheology(v1,v2)
julia> PlotDeformationMap(v, levels=100, colormap=:roma)
```


Next, let's plot viscosity and flip x & y axis:

```julia
julia> PlotDeformationMap(v, viscosity=true, rotate_axes=true)
```


Instead of plotting stress vs. T and computing strainrate, we can also provide strainrate/T and compute stress:

```julia
julia> PlotDeformationMap(v, strainrate=false)
```


Or plot viscosity but only add contours in a certain range:

```julia
julia> PlotDeformationMap(v,  strainrate=false, viscosity=true, levels=Vector(18:.25:24))
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Plotting/Plotting.jl#L1026-L1067" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.PlotDiffusionCoefArrhenius' href='#GeoParams.PlotDiffusionCoefArrhenius'><span class="jlbinding">GeoParams.PlotDiffusionCoefArrhenius</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fig, ax = PlotDiffusionCoefPlotDiffusionCoefArrhenius(x::Union{Tuple{Vararg{AbstractChemicalDiffusion}}, NTuple{N, AbstractChemicalDiffusion} where N, AbstractChemicalDiffusion};
                            P=1u"GPa", fO2=1.0e-25NoUnits, log_type=:log10, linestyle=:solid, linewidth=1, color=nothing, label=nothing,
                            title="", fig=nothing, filename=nothing, res=(1200, 1200), legend=true, legendsize=15, position=:rt,
                            labelsize=35, xlims=(nothing, nothing), ylims=(nothing, nothing))
```


Creates a plot of log(D) versus 10^4/T for one or a tuple of `ChemicalDiffusionData` structures.

**Optional parameters**
- `P`: Pressure (default: 1 GPa)
  
- `fO2`: Oxygen fugacity (default: 1.0e-25 NoUnits (Graphite buffer))
  
- `log_type`: Logarithm type for `D` (default: :log10, options: :log10, :ln)
  
- `linestyle`: Line style for the plot (default: :solid)
  
- `linewidth`: Line width for the plot (default: 1)
  
- `color`: Line color for the plot (default: nothing)
  
- `label`: Label for the plot (default: nothing)
  
- `title`: Title for the plot (default: "")
  
- `fig`: Existing figure to plot on (default: nothing)
  
- `filename`: Filename to save the plot (default: nothing)
  
- `res`: Resolution of the plot (default: (1200, 1200))
  
- `legend`: Whether to display the legend (default: true)
  
- `legendsize`: Size of the legend text (default: 15)
  
- `position`: Position of the legend (default: :rt)
  
- `labelsize`: Size of the axis labels (default: 35)
  
- `xlims`: Limits for the x-axis (default: (nothing, nothing))
  
- `ylims`: Limits for the y-axis (default: (nothing, nothing))
  
- `ticklabelsize`: Size of the tick labels (default: 35)
  

**Example**

```julia
using GeoParams
using CairoMakie

# obtain diffusion data
Fe_Grt = Garnet.Grt_Fe_Chakraborty1992
Fe_Grt = SetChemicalDiffusion(Fe_Grt)
Mg_Grt = Garnet.Grt_Mg_Chakraborty1992
Mg_Grt = SetChemicalDiffusion(Mg_Grt)
Mn_Grt = Garnet.Grt_Mn_Chakraborty1992
Mn_Grt = SetChemicalDiffusion(Mn_Grt)

fig, ax = PlotDiffusionCoefArrhenius((Fe_Grt, Mg_Grt, Mn_Grt), P= 1u"GPa", linewidth=3,
                                     label= ("Fe from Chakraborty and Ganguly (1992)",
                                             "Mg from Chakraborty and Ganguly (1992)",
                                             "Mn from Chakraborty and Ganguly (1992)"))
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Plotting/Plotting.jl#L1378-L1425" target="_blank" rel="noreferrer">source</a></Badge>

</details>

