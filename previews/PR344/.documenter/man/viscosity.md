---
---

# Effective viscosity {#Effective-viscosity}

To compute the viscosity given stress of strainrate use:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.compute_viscosity_εII' href='#GeoParams.compute_viscosity_εII'><span class="jlbinding">GeoParams.compute_viscosity_εII</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_viscosity_εII(s::AbstractConstitutiveLaw, εII, kwargs...)
```


Compute effective viscosity given a 2nd invariant of the deviatoric strain rate tensor, extra parameters are passed as a named tuple, e.g., (;T=T)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Viscosity/Viscosity.jl#L5-L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.compute_viscosity_τII' href='#GeoParams.compute_viscosity_τII'><span class="jlbinding">GeoParams.compute_viscosity_τII</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_viscosity_τII(s::AbstractConstitutiveLaw, τII, kwargs...)
```


Compute effective viscosity given a 2nd invariant of the deviatoric stress tensor and, extra parameters are passed as a named tuple, e.g., (;T=T)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/1a8e4f8dddc2ccb8a3f977ae466234307cc4b304/src/Viscosity/Viscosity.jl#L17-L21" target="_blank" rel="noreferrer">source</a></Badge>

</details>

