---
---

## Olivine {#Olivine}

To obtain the list of implemented diffusion parameters for olivine, use:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.chemical_diffusion_list' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.chemical_diffusion_list'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.chemical_diffusion_list</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
chemical_diffusion_list(search::String="")
```


List all available chemical diffusion data for rutile. Includes an argument to search for a specific term, i.e. an element ("La") or an author.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Olivine/Olivine.jl#L9-L14" target="_blank" rel="noreferrer">source</a></Badge>

</details>


The following diffusion parameters for olivine are implemented:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Mg_Chakraborty1994_forsterite' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Mg_Chakraborty1994_forsterite'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Mg_Chakraborty1994_forsterite</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Ol_Mg_Chakraborty1994_forsterite()
```


Trace diffusion data of Mg in synthetic pure forsterite. Calibrated with experiments conducted between 1000-1300°C at atmospheric pressure and up to 10 GPa in multianvil apparatus. From Chakraborty and al. (1994) (https://doi.org/10.1007/BF00203923).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Olivine/Elements/Mg.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Fe_Mg_Dohmen2007_perp_c' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Fe_Mg_Dohmen2007_perp_c'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Fe_Mg_Dohmen2007_perp_c</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Ol_Fe_Mg_Dohmen2007_perp_c()
```


Interdiffusion data of Fe and Mg in natural gem quality olivine for all oxygen fugacities. Calibrated with experiments conducted between 700-1200°C at atmospheric pressure with oxygen fugacity between 1e-12 to 1e-5 and parallel to c axis from Dohmen et al. (2007) and data from the literature. Note that the molar fraction of Fe is required to calculate the interdiffusion coefficient. To calculate the interdiffusion coefficient parallel to c, multiply the interdiffusion coefficient by 1e6. From Dohmen and Chakraborty (2007) (https://doi.org/10.1007/s00269-007-0158-6).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Olivine/Elements/Fe_Mg.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Fe_Mg_Dohmen2007_TaMED_perp_c' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Fe_Mg_Dohmen2007_TaMED_perp_c'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Fe_Mg_Dohmen2007_TaMED_perp_c</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Ol_Fe_Mg_Dohmen2007_TaMED_perp_c()
```


Interdiffusion data of Fe and Mg for the transition metal extrinsic (TaMED) mechanism. Calibrated with experiments conducted between 900-1200°C at atmospheric pressure with oxygen fugacity between 1e-10 to 1e-5 and parallel to c axis from Dohmen et al., 2007 and data from the literature. Note that the molar fraction of Fe and the oxygen fugacity are required to calculate the interdiffusion coefficient. To calculate the interdiffusion coefficient parallel to c, multiply the interdiffusion coefficient by 1e6. From Dohmen and Chakraborty (2007) (https://doi.org/10.1007/s00269-007-0158-6).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Olivine/Elements/Fe_Mg.jl#L45-L50" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Fe_Mg_Dohmen2007_PED_perp_c' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Fe_Mg_Dohmen2007_PED_perp_c'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Olivine.Ol_Fe_Mg_Dohmen2007_PED_perp_c</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Ol_Fe_Mg_Dohmen2007_PED_perp_c()
```


Interdiffusion data of Fe and Mg in natural gem quality olivine for the purely extrinsic (PED) mechanism. Calibrated with experiments conducted between 700-900°C at atmospheric pressure with oxygen fugacity between 1e-12 to 1e-10 and parallel to c axis from Dohmen et al. (2007) and data from the literature. Note that the molar fraction of Fe is required to calculate the interdiffusion coefficient. To calculate the interdiffusion coefficient parallel to c, multiply the interdiffusion coefficient by 1e6. From Dohmen and Chakraborty (2007) (https://doi.org/10.1007/s00269-007-0158-6).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/ChemicalDiffusion/Data/Olivine/Elements/Fe_Mg.jl#L92-L97" target="_blank" rel="noreferrer">source</a></Badge>

</details>

