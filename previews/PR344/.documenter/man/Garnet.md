---
---

## Garnet {#Garnet}

To obtain the list of implemented diffusion parameters for garnet, use:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.chemical_diffusion_list' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.chemical_diffusion_list'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.chemical_diffusion_list</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
chemical_diffusion_list(search::String="")
```


List all available chemical diffusion data for garnet. Includes an argument to search for a specific term, i.e. an element ("La") or an author.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/ChemicalDiffusion/Data/Garnet/Garnet.jl#L11-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>


The following diffusion parameters for garnet are implemented:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_Mg_Chakraborty1992' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_Mg_Chakraborty1992'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_Mg_Chakraborty1992</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Grt_Mg_Chakraborty1992()
```


Diffusion data of Mg in garnet. Calibrated between 1100-1480°C and 0.14-0.43 GPa in the C-O2 system. From Chakraborty and Ganguly (1992) (https://doi.org/10.1007/BF00296579) combined with data from Loomis et al. (1985) (https://doi.org/10.1007/BF00373040).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/ChemicalDiffusion/Data/Garnet/Elements/Mg.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_Mn_Chakraborty1992' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_Mn_Chakraborty1992'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_Mn_Chakraborty1992</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Grt_Mn_Chakraborty1992()
```


Diffusion data of Mn in garnet. Calibrated between 1100-1480°C and 0.14-0.43 GPa in the C-O2 system. From Chakraborty and Ganguly (1992) (https://doi.org/10.1007/BF00296579) combined with data from Loomis et al. (1985) (https://doi.org/10.1007/BF00373040).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/ChemicalDiffusion/Data/Garnet/Elements/Mn.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_Fe_Chakraborty1992' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_Fe_Chakraborty1992'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_Fe_Chakraborty1992</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Grt_Fe_Chakraborty1992()
```


Diffusion data of Fe in garnet. Calibrated between 1100-1480°C and 0.14-0.43 GPa in the C-O2 system. From Chakraborty and Ganguly (1992) (https://doi.org/10.1007/BF00296579) combined with data from Loomis et al. (1985) (https://doi.org/10.1007/BF00373040).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/ChemicalDiffusion/Data/Garnet/Elements/Fe.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_REE_Bloch2020_slow' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_REE_Bloch2020_slow'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_REE_Bloch2020_slow</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Grt_REE_Bloch2020_slow()
```


Diffusion data of REE in garnet for slow diffusion mechanism. Calibrated with experiments conducted between 950-1050°C at 1 atm with the QFM buffer and dry condition plus data from the literature. From Bloch et al. (2020) (https://doi.org/10.1093/petrology/egaa055) combined with data from Van Orman et al. (2002), Tirone et al. (2005) and Bloch et al. (2015).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/ChemicalDiffusion/Data/Garnet/Elements/REE.jl#L46-L51" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_REE_Bloch2020_fast' href='#GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_REE_Bloch2020_fast'><span class="jlbinding">GeoParams.MaterialParameters.ChemicalDiffusion.Garnet.Grt_REE_Bloch2020_fast</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
Grt_REE_Bloch2020_fast()
```


Diffusion data of REE in garnet for fast diffusion mechanism. Calibrated with experiments conducted between 950-1050°C at 1 atm with the QFM buffer and dry condition plus data from the literature. From Bloch et al. (2020) (https://doi.org/10.1093/petrology/egaa055) combined with data from Van Orman et al. (2002), Tirone et al. (2005) and Bloch et al. (2015).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/ChemicalDiffusion/Data/Garnet/Elements/REE.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

