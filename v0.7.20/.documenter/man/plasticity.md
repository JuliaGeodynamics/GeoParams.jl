---
---

# Plasticity {#Plasticity}

Plasticity is a non-linear rheology that is activated once stresses exceed a certain yield criteria.

## Implemented laws {#Implemented-laws}

The following plastic law are implemented:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPrager' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPrager'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPrager</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DruckerPrager(ϕ=30, Ψ=0, C=10e6Pa)
```


Sets parameters for Drucker-Prager plasticity, where the yield stress $\sigma_{y}$ is computed by

$$    \sigma_{y} = (P-P_f)\tan(ϕ) + C$$

with $\phi$ being the friction angle (in degrees), $C$ cohesion, $P$ dynamic pressure and $P_f$ the fluid pressure (both positive under compression).

_Yielding_ occurs when the second invariant of the deviatoric stress tensor, $\tau_{II}=(0.5\tau_{ij}\tau_{ij})^{0.5}$ touches the yield stress. This can be computed with the yield function $F$ and the plastic flow potential $Q$, which are respectively given by

$$    F = \tau_{II} - \cos(ϕ)C - \sin(ϕ)(P-P_f)$$

$$    Q = \tau_{II} - \sin(Ψ)(P-P_f)$$

Here, Ψ is the dilation angle, which must be zero for incompressible setups.

Plasticity is activated when $F(\tau_{II}^{trial})$ (the yield function computed with a trial stress) is &gt;0. In that case, plastic strainrate $\dot{\varepsilon}^{pl}_{ij}$ is computed by:

$$    \dot{\varepsilon}^{pl}_{ij} = \dot{\lambda} \frac{\partial Q}{\partial \sigma_{ij}}$$

where $\dot{\lambda}$ is a (scalar) that is nonzero and chosen such that the resulting stress gives $F(\tau_{II}^{final})=0$, and $\sigma_{ij}=-P + \tau_{ij}$ denotes the total stress tensor.

**Example**

```julia
julia> pl = DruckerPrager(ϕ=30, C=10e6Pa)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/Plasticity/DruckerPrager.jl#L4-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPrager_regularised' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPrager_regularised'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPrager_regularised</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DruckerPrager_regularised(ϕ=30, Ψ=0, C=10e6Pa, η_vp=1e20Pa*s)
```


Sets parameters for reularised Drucker-Prager plasticity, where the yield stress $\sigma_{y}$ is computed by

$$    \sigma_{y} = (P-P_f)\tan(ϕ) + C + 2η_vpε̇II_pl$$

with $\phi$ being the friction angle (in degrees), $C$ cohesion, $P$ dynamic pressure, $P_f$ the fluid pressure (both positive under compression), $η_vp$ the regularization viscosity and $ε̇II_pl$ the invariant of the plastic strainrate

_Yielding_ occurs when the second invariant of the deviatoric stress tensor, $\tau_{II}=(0.5\tau_{ij}\tau_{ij})^{0.5}$ touches the yield stress. This can be computed with the yield function $F$ and the plastic flow potential $Q$, which are respectively given by

$$    F = \tau_{II} - \cos(ϕ)C - \sin(ϕ)(P-P_f) - 2 \eta_{vp} \dot{\varepsilon}ε̇^{pl}_{II}$$

$$    Q = \tau_{II} - \sin(Ψ)(P-P_f)$$

Here, Ψ is the dilation angle, which must be zero for incompressible setups.

Plasticity is activated when $F(\tau_{II}^{trial})$ (the yield function computed with a trial stress) is &gt;0. In that case, plastic strainrate $\dot{\varepsilon}^{pl}_{ij}$ is computed by:

$$    \dot{\varepsilon}^{pl}_{ij} = \dot{\lambda} \frac{\partial Q}{\partial \sigma_{ij}}$$

where $\dot{\lambda}$ is a (scalar) that is nonzero and chosen such that the resulting stress gives $F(\tau_{II}^{final})=0$, and $\sigma_{ij}=-P + \tau_{ij}$ denotes the total stress tensor.

**Example**

```julia
julia> pl = DruckerPrager_regularised(ϕ=30, C=10e6Pa, η_vp=1e20Pa*s)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/Plasticity/DruckerPrager_regularised.jl#L4-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPragerCap' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPragerCap'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPragerCap</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DruckerPragerCap(ϕ=30, Ψ=0, C=10e6Pa, η_vp=1e20Pa*s, Pt=-1e5Pa)
```


Sets parameters for Drucker-Prager-Cap plasticity, as described in Popov et al. (2025), Geoscientific Model Development. The yield surface has two branches selected by the stress state: a mode-2 Drucker-Prager shear branch at high pressure and a mode-1 tensile cap near the tensile strength `pT`. The yield function $F$ and plastic flow potential $Q$ are

$$    F = \tau_{II} - kP - c \quad\text{or}\quad a\left(\sqrt{\tau_{II}^2 + (P-p_y)^2} - R_y\right)$$

$$    Q = \tau_{II} - k_q P - \text{const} \quad\text{or}\quad b\left(\sqrt{\tau_{II}^2 + (P-p_q)^2} - R_f\right)$$

where the first form is the shear branch and the second the tensile cap. As in [`DruckerPrager`](/man/plasticity#GeoParams.MaterialParameters.ConstitutiveRelationships.DruckerPrager), plastic strain rate is $\dot{\varepsilon}^{pl}_{ij} = \dot{\lambda}\,\partial Q/\partial \sigma_{ij}$, with `η_vp` providing Duvaut-Lions (Duretz-type) viscoplastic regularisation.

**Fields**
- `C::T`: The cohesion parameter.
  
- `ϕ::T`: The friction angle (in degrees).
  
- `Ψ::T`: The dilatancy angle (in degrees).
  
- `η_vp::T`: The Duvaut-Lions regularisation viscosity for the plasticity model.
  
- `pT::T`: The tensile strength (should be &lt; 0).
  

**Example**

```julia
julia> pl = DruckerPragerCap(ϕ=30, C=10e6Pa, pT=-1e5Pa)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/20a6a4c5be2e0068300cfd4bcc0c7f50d1e717e0/src/Plasticity/DruckerPragerCap.jl#L6-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Computational routines {#Computational-routines}

Usually, plasticity should be defined as part of a `CompositeRheology` structure and calculations can be done as with all other rheology computations by using `compute_τII`.
