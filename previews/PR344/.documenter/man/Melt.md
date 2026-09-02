---
---

## Melt {#Melt}

To obtain the list of implemented diffusion parameters for melt, use:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.chemical_diffusion_list' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.chemical_diffusion_list'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.chemical_diffusion_list</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
chemical_diffusion_list(search::String="")
```


List all available chemical diffusion data for melt. Includes an argument to search for a specific term, i.e. an element ("La") or an author.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Melt.jl#L33-L38" target="_blank" rel="noreferrer">source</a></Badge>

</details>


For major element diffusion, the multicomponent diffusion matrix for basaltic melt can be calculated using the following function:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_multicomponent_major_Guo2020_SiO2_basaltic' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_multicomponent_major_Guo2020_SiO2_basaltic'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_multicomponent_major_Guo2020_SiO2_basaltic</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_multicomponent_major_Guo2020_SiO2_basaltic()
```


Multicomponent diffusion data of major elements in basaltic melt. Calibrated using 18 diffusion couple experiments for SiO2–TiO2–Al2O3–FeO–MgO–CaO–Na2O–K2O basaltic melt compositions and data from Guo and Zhang (2018). SiO2 is here the dependent variable. This contains the invariant eigenvectors and the diagonal matrices containing the pre-exponential factors and activation energies of the eigenvalues to compute the diffusion matrix at a given temperature. It is supposed here that the eigenvectors are not temperature dependent but eigenvalues are. When computing the diffusion matrix, the order of the species is TiO2, Al2O3, FeO, MgO, CaO, Na2O, K2O. From Guo and Zhang (2020) (https://10.1016/j.chemgeo.2020.119700).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Multicomponent.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>


For this particular case, to initiate the parameters, use `SetMulticompChemicalDiffusion`:

```julia
    melt_multicomponent = Melt.Melt_multicomponent_major_Guo2020_SiO2_basaltic
    melt_multicomponent = SetMulticompChemicalDiffusion(melt_multicomponent)
```


The function `compute_D` can be used as usual to calculate the diffusion coefficients but will in this case return a [static](https://juliaarrays.github.io/StaticArrays.jl/dev/) matrix of diffusion coefficients for the major elements in the melt:

```julia
    D = compute_D(melt_multicomponent, T=1000C)
7×7 SMatrix{7, 7, Quantity{Float64, 𝐋²·⁰ 𝐓⁻¹·⁰, Unitful.FreeUnits{(m²·⁰, s⁻¹·⁰), 𝐋²·⁰ 𝐓⁻¹·⁰, nothing}}, 49} with indices SOneTo(7)×SOneTo(7):
  1.03519e-13 m²·⁰ s⁻¹·⁰   8.93006e-13 m²·⁰ s⁻¹·⁰  -1.95924e-12 m²·⁰ s⁻¹·⁰  …  -5.84922e-12 m²·⁰ s⁻¹·⁰   8.53056e-12 m²·⁰ s⁻¹·⁰   2.43836e-12 m²·⁰ s⁻¹·⁰
  4.41504e-14 m²·⁰ s⁻¹·⁰   5.07365e-15 m²·⁰ s⁻¹·⁰   1.91304e-13 m²·⁰ s⁻¹·⁰      1.16129e-12 m²·⁰ s⁻¹·⁰  -1.31307e-12 m²·⁰ s⁻¹·⁰  -4.16234e-13 m²·⁰ s⁻¹·⁰
 -2.91173e-13 m²·⁰ s⁻¹·⁰  -1.74973e-13 m²·⁰ s⁻¹·⁰   2.62782e-12 m²·⁰ s⁻¹·⁰      2.77994e-12 m²·⁰ s⁻¹·⁰  -5.25369e-12 m²·⁰ s⁻¹·⁰  -1.49885e-12 m²·⁰ s⁻¹·⁰
 -4.8815e-13 m²·⁰ s⁻¹·⁰    2.53181e-12 m²·⁰ s⁻¹·⁰  -3.84208e-12 m²·⁰ s⁻¹·⁰     -1.54711e-11 m²·⁰ s⁻¹·⁰   2.16066e-11 m²·⁰ s⁻¹·⁰   5.19873e-12 m²·⁰ s⁻¹·⁰
  6.67519e-13 m²·⁰ s⁻¹·⁰  -1.47328e-12 m²·⁰ s⁻¹·⁰  -5.84858e-13 m²·⁰ s⁻¹·⁰      8.55322e-12 m²·⁰ s⁻¹·⁰  -7.44916e-12 m²·⁰ s⁻¹·⁰  -2.12405e-12 m²·⁰ s⁻¹·⁰
 -4.89738e-13 m²·⁰ s⁻¹·⁰   2.07724e-12 m²·⁰ s⁻¹·⁰  -2.8383e-12 m²·⁰ s⁻¹·⁰   …  -1.10541e-11 m²·⁰ s⁻¹·⁰   1.78501e-11 m²·⁰ s⁻¹·⁰   4.10209e-12 m²·⁰ s⁻¹·⁰
 -2.59574e-13 m²·⁰ s⁻¹·⁰   7.36022e-13 m²·⁰ s⁻¹·⁰  -5.65377e-13 m²·⁰ s⁻¹·⁰     -2.34514e-12 m²·⁰ s⁻¹·⁰   3.67054e-12 m²·⁰ s⁻¹·⁰   3.13608e-12 m²·⁰ s⁻¹·⁰
```


Additionally, the function `compute_λ` can be used directly to compute the diagonal matrix of the eigenvalues. This can be useful in the case when the diffusion matrix needs to be diagonalized:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.compute_λ' href='#GeoParams.MaterialParameters.ChemicalDiffusion.compute_λ'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.compute_λ</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_λ(data::MeltMulticompDiffusionData; T=1K, kwargs...)
```


Computes the diagonal matrix of eigenvalues [m^2/s] from the diffusion data `data` at temperature `T` [K] from a structure of type `MeltMulticompDiffusionData`. This is useful if the system needs to be diagonalized. If `T` is provided without unit, the function assumes the unit is in Kelvin and outputs the eigenvalues without unit based on the value in m^2/s. The output is a static matrix of size `n-1` x `n-1` where `n` is the number of components.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/ChemicalDiffusion.jl#L393-L400" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.compute_λ!' href='#GeoParams.MaterialParameters.ChemicalDiffusion.compute_λ!'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.compute_λ!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_λ!(λ::AbstractArray, data::MeltMulticompDiffusionData; T = ones(size(λ))K, kwargs...)
```


In-place version of `compute_λ(data::AbstractChemicalDiffusion; T=1K, kwargs...)`. `λ` should be an array of the same size as T.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/ChemicalDiffusion.jl#L419-L423" target="_blank" rel="noreferrer">source</a></Badge>

</details>


For trace-, self- or interdiffusion parameters for melt, the following functions are implemented:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Mg_Sheng1992_basaltic' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Mg_Sheng1992_basaltic'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Mg_Sheng1992_basaltic</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Mg_Sheng1992_basaltic()
```


Diffusion data of Mg self-diffusion in basaltic melt. Calibrated with experiments with Spinel and anhydrous synthetic glass between 1250 to 1550°C at ambient pressure. From Sheng et al. (1992) (https://doi.org/10.1016/0016-7037(92)90207-Y).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Mg.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sc_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sc_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sc_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Sc_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Sc in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Sc.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sc_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sc_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sc_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Sc_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Sc in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Sc.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_V_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_V_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_V_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_V_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of V in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/V.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_V_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_V_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_V_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_V_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of V in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/V.jl#L21-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Y_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Y_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Y_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Y_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Y in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Y.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Y_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Y_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Y_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Y_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Y in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Y.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Zr_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Zr_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Zr_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Zr_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Zr in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Zr.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Zr_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Zr_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Zr_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Zr_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Zr in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Zr.jl#L21-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Hf_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Hf_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Hf_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Hf_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Hf in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Hf.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nb_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nb_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nb_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Nb_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Nb in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Nb.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nb_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nb_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nb_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Nb_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Nb in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Nb.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_La_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_La_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_La_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_La_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of La in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ni and Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/La.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_La_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_La_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_La_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_La_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of La in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni and Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/La.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ce_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ce_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ce_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Ce_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Ce in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Ce.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ce_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ce_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ce_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Ce_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Ce in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Ce.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Pr_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Pr_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Pr_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Pr_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Pr in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Pr.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Pr_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Pr_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Pr_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Pr_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Pr in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Pr.jl#L21-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nd_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nd_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nd_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Nd_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Nd in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Nd.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nd_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nd_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Nd_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Nd_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Nd in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Nd.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sm_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sm_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sm_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Sm_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Sm in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Sm.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sm_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sm_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Sm_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Sm_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Sm in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Sm.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Eu_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Eu_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Eu_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Eu_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Eu in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Eu.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Gd_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Gd_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Gd_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Gd_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Gd in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Gd.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Gd_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Gd_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Gd_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Gd_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Gd in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Gd.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Tb_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Tb_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Tb_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Tb_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Tb in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Tb.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Tb_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Tb_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Tb_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Tb_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Tb in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Tb.jl#L21-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Dy_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Dy_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Dy_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Dy_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Dy in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Dy.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Dy_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Dy_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Dy_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Dy_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Dy in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Dy.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ho_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ho_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ho_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Ho_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Ho in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Ho.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ho_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ho_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Ho_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Ho_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Ho in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Ho.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Er_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Er_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Er_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Er_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Er in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Er.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Er_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Er_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Er_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Er_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Er in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Er.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Yb_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Yb_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Yb_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Yb_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Yb in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Yb.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Yb_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Yb_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Yb_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Yb_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Yb in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Yb.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Lu_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Lu_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Lu_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Lu_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Lu in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Lu.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Lu_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Lu_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Lu_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Lu_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of Lu in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Lu.jl#L22-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Th_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Th_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_Th_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_Th_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of Th in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/Th.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_U_Holycross2018_rhyolitic_highH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_U_Holycross2018_rhyolitic_highH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_U_Holycross2018_rhyolitic_highH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_U_Holycross2018_rhyolitic_highH2O()
```


Diffusion data of U in rhyolitic melt (76.77 wt% SiO2) with 6.2 wt% of H2O. Calibrated with experiments conducted between 850-935°C at 1 GPa with Ag capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/U.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_U_Holycross2018_rhyolitic_mediumH2O' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_U_Holycross2018_rhyolitic_mediumH2O'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Melt.Melt_U_Holycross2018_rhyolitic_mediumH2O</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Melt_U_Holycross2018_rhyolitic_mediumH2O()
```


Diffusion data of U in rhyolitic melt (76.77 wt% SiO2) with 4.1 wt% of H2O. Calibrated with experiments conducted between 960-1250°C at 1 GPa with Ni capsules from synthetic glass. From Holycross and Watson (2018) (https://doi.org/10.1016/j.gca.2018.04.006).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Melt/Elements/U.jl#L21-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

