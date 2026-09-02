---
---

# Shear heating {#Shear-heating}

## Methods {#Methods}

Heat caused by non-recoverable deformation can be specified in 
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Shearheating.ConstantShearheating' href='#GeoParams.MaterialParameters.Shearheating.ConstantShearheating'><span class="jlbinding">GeoParams.MaterialParameters.Shearheating.ConstantShearheating</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstantShearheating(Χ=0.0NoUnits)
```


Set the shear heating efficiency [0-1] parameter

$$Χ  = cst$$

where $\Chi$ is the shear heating efficiency [NoUnits]

Shear heating is computed as

$$H_s = \Chi \cdot \tau_{ij}(\dot{\varepsilon}_{ij} - \dot{\varepsilon}^{el}_{ij})$$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Energy/Shearheating.jl#L23-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Computational routines {#Computational-routines}

To compute, use this:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Shearheating.compute_shearheating' href='#GeoParams.MaterialParameters.Shearheating.compute_shearheating'><span class="jlbinding">GeoParams.MaterialParameters.Shearheating.compute_shearheating</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
H_s = compute_shearheating(s:<AbstractShearheating, τ, ε, ε_el)
```


Computes the shear heating source term

$$H_s = \Chi \cdot \tau_{ij} ( \dot{\varepsilon}_{ij} - \dot{\varepsilon}^{el}_{ij})$$

**Parameters**
- $\Chi$ : The efficiency of shear heating (between 0-1)
  
- $\tau_{ij}$ : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
  
- $\dot{\varepsilon}_{ij}$ : The full deviatoric strainrate tensor
  
- $\dot{\varepsilon}^{el}_{ij}$ : The full elastic deviatoric strainrate tensor
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Energy/Shearheating.jl#L91-L106" target="_blank" rel="noreferrer">source</a></Badge>



```julia
H_s = ComputeShearheating(s:<AbstractShearheating, τ, ε)
```


Computes the shear heating source term when there is no elasticity

$$H_s = \Chi \cdot \tau_{ij}  \dot{\varepsilon}_{ij}$$

**Parameters**
- $\Chi$ : The efficiency of shear heating (between 0-1)
  
- $\tau_{ij}$ : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
  
- $\dot{\varepsilon}_{ij}$ : The full deviatoric strainrate tensor
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Energy/Shearheating.jl#L109-L122" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Shearheating.compute_shearheating!' href='#GeoParams.MaterialParameters.Shearheating.compute_shearheating!'><span class="jlbinding">GeoParams.MaterialParameters.Shearheating.compute_shearheating!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_shearheating!(H_s, s:<AbstractShearheating,  τ, ε, ε_el)
```


Computes the shear heating source term in-place

$$H_s = \Chi \cdot \tau_{ij} ( \dot{\varepsilon}_{ij} - \dot{\varepsilon}^{el}_{ij})$$

**Parameters**
- $\Chi$ : The efficiency of shear heating (between 0-1)
  
- $\tau_{ij}$ : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
  
- $\dot{\varepsilon}_{ij}$ : The full deviatoric strainrate tensor
  
- $\dot{\varepsilon}^{el}_{ij}$ : The full elastic deviatoric strainrate tensor
  

_NOTE:_ The shear heating terms require the full deviatoric stress & strain rate tensors, i.e.:

$$2D: \tau_{ij} = \left(
                \begin{matrix}
                    \tau_{xx} & \tau_{xz} \\
                    \tau_{zx} & \tau_{zz}
                \end{matrix}
            \right)$$

Since $\tau_{zx}=\tau_{xz}$, most geodynamic codes only take one of the terms into account; shear heating requires all components to be used!


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Energy/Shearheating.jl#L127-L154" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_shearheating!(H_s, s:<AbstractShearheating, τ, ε)
```


Computes the shear heating source term `H_s` in-place when there is no elasticity

$$H_s = \Chi \cdot \tau_{ij}  \dot{\varepsilon}_{ij}$$

**Parameters**
- $\Chi$ : The efficiency of shear heating (between 0-1)
  
- $\tau_{ij}$ : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
  
- $\dot{\varepsilon}_{ij}$ : The full deviatoric strainrate tensor
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Energy/Shearheating.jl#L157-L170" target="_blank" rel="noreferrer">source</a></Badge>

</details>

