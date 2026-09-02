---
---

# Chemical Diffusion {#Chemical-Diffusion}

Some routines and experimental data to compute chemical diffusion coefficients in minerals and phases are implemented in [GeoParams.jl](https://juliageodynamics.github.io/GeoParams.jl/stable/).

## Diffusion parameters from the literature {#Diffusion-parameters-from-the-literature}

Currently, four phases are implemented in independent modules: `Rutile`, `Olivine`, `Garnet`, and `Melt`. Each phase has a list of diffusion parameters from the literature for different chemical elements. Other phases will be implemented in the future and contributions are welcome to extend the database.

To initiate the diffusion parameters of an element of a phase, call the function `SetChemicalDiffusion`. For instance, to obtain the diffusion parameters of Hf in rutile, use:

```julia
using GeoParams

Hf_Rt_para = Rutile.Rt_Hf_Cherniak2007_para_c
Hf_Rt_para = SetChemicalDiffusion(Hf_Rt_para)
```


`Hf_Rt_para` is in this case a structure of type `DiffusionData` containing the diffusion parameters for Hf in rutile, from Cherniak et al. (2007).

To know the available diffusion parameters for each phase, one may use `X.chemical_diffusion_list()`, where `X` is the phase of interest. This will return a vector containing the functions available.

For instance, for garnet, using:

```julia
Garnet.chemical_diffusion_list()
```


will return:

```julia
13-element Vector{Function}:
Grt_Ca_Carlson2006 (generic function with 1 method)
 Grt_Ca_Chu2015 (generic function with 1 method)
 Grt_Fe_Carlson2006 (generic function with 1 method)
 Grt_Fe_Chakraborty1992 (generic function with 1 method)
 Grt_Fe_Chu2015 (generic function with 1 method)
 ⋮
 Grt_Mn_Carlson2006 (generic function with 1 method)
 Grt_Mn_Chakraborty1992 (generic function with 1 method)
 Grt_Mn_Chu2015 (generic function with 1 method)
 Grt_REE_Bloch2020_fast (generic function with 1 method)
 Grt_REE_Bloch2020_slow (generic function with 1 method)
```


## Diffusion parameters {#Diffusion-parameters}

A diffusion coefficient $D$ can be computed assuming the following general equation:

$$    D = D_0 \left(\frac{{fO_2}}{d_{f{O_2}}}^{n_{f{O_2}}}\right) \exp(a_X * (X + b_X)) \exp \left(-\frac{Ea + (P - P0) \Delta V} {RT}\right),$$

where $D_0$ is the pre-exponential factor, $f{O_2}$ is the oxygen fugacity, $d_{fO_2}$ and $n_{fO_2}$ are the parameters on the $f{O_2}$ dependency, $X$ is a composition or activity, $a_X$ and $b_X$ are parameters, $Ea$ is the activation energy, $P$ is the pressure, $P0$ is the reference pressure, $\Delta V$ is the activation volume, $R$ is the gas constant, and $T$ is the temperature.

To compute the diffusion coefficient in-place, use the following functions:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.compute_D' href='#GeoParams.MaterialParameters.ChemicalDiffusion.compute_D'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.compute_D</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_D(data::DiffusionData; T=1K, P=1GPa, fO2 = 1e-25NoUnits, X = 0 NoUnits, kwargs...)
```


Computes the diffusion coefficient `D` [m^2/s] from the diffusion data `data` at temperature `T` [K], pressure `P` [Pa], oxygen fugacity `fO2` [NoUnits] and composition dependency `X` [NoUnits] from a structure of type `DiffusionData`. If `T` and `P` are provided without unit, the function assumes the units are in Kelvin and Pascal, respectively, and outputs the diffusion coefficient without unit based on the value in m^2/s.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/ChemicalDiffusion/ChemicalDiffusion.jl#L321-L326" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_D(data::MeltMulticompDiffusionData; T=1K, kwargs...)
```


Computes the diffusion coefficient matrix [m^2/s] from the diffusion data `data` at temperature `T` [K] from a structure of type `MeltMulticompDiffusionData`. If `T` is provided without unit, the function assumes the unit is in Kelvin and outputs the diffusion coefficient without unit based on the value in m^2/s. The output is a static matrix of size `n-1` x `n-1` where `n` is the number of components.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/ChemicalDiffusion/ChemicalDiffusion.jl#L346-L352" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.compute_D!' href='#GeoParams.MaterialParameters.ChemicalDiffusion.compute_D!'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.compute_D!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_D!(D, data::AbstractChemicalDiffusion; T=ones(size(D))K, P=ones(size(D))Pa, fO2 = ones(size(D)), X = zeros(size(D)), kwargs...)
```


In-place version of `compute_D(data::AbstractChemicalDiffusion; T=1K, P=1GPa, fO2=0NoUnits, kwargs...)`. `D` should be an array of the same size as T, P and fO2.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/ChemicalDiffusion/ChemicalDiffusion.jl#L372-L376" target="_blank" rel="noreferrer">source</a></Badge>

</details>


For instance, taking the previous example with Hf in Rutile:

```julia
D = compute_D(Hf_Rt_para, T=1000C)
```


## Plotting routines {#Plotting-routines}

It is often convenient to plot the diffusion coefficients vs 1000/T to compare different diffusion parameters from the literature and to check the implementation of the diffusion parameters.

A function (`GeoParams.PlotDiffusionCoefArrhenius`) for that purpose is implemented using [Makie.jl](https://docs.makie.org/stable/).

For instance, plotting the diffusion coefficients of Fe, Mg and Mn from Chakraborty and Ganguly (1992) can be done with:

```julia
using GeoParams
using CairoMakie  # need to import explicitly a Makie backend for the extension to be loaded

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


which produces the following plot:


![](assets/img/ChemicalDiffusion_garnetdiffusion.png)


More details on the function `GeoParams.PlotDiffusionCoefArrhenius` and its options can be found in the section [Plotting](/man/plotting#Plotting).
