---
---

# CreepLaws {#CreepLaws}

The following viscous creep laws are implemented:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.LinearViscous' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.LinearViscous'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.LinearViscous</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
LinearViscous(η=1e20Pa*s)
```


Defines a linear viscous creeplaw

The (isotopic) linear viscous rheology is given by

$$    \tau_{ij} = 2 \eta \dot{\varepsilon}_{ij}$$

or

$$    \dot{\varepsilon}_{ij} = \frac{\tau_{ij}}{2 \eta}$$

where $\eta_0$ is the reference viscosity [Pa*s] at reference strain rate $\dot{\varepsilon}_0$[1/s], and $n$ the power law exponent [].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/CreepLaw.jl#L104-L119" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.PowerlawViscous' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.PowerlawViscous'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.PowerlawViscous</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PowerlawViscous(η0=1e18Pa*s, n=2.0NoUnits, ε0=1e-15/s)
```


Defines a power law viscous creeplaw as:

$$        \tau_{ij}^n  = 2 \eta_0 \frac{\dot{\varepsilon}_{ij}}{\dot{\varepsilon}_0}$$

where $\eta$ is the effective viscosity [Pa*s].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/CreepLaw.jl#L341-L350" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.DislocationCreep' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.DislocationCreep'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.DislocationCreep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DislocationCreep(n = 1.0NoUnits, r = 0.0NoUnits, A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = AxialCompression )
```


Defines the flow law parameter of a dislocation creep law.

The (isotropic) dislocation creep law, as used by experimentalists, is given by

$$     \dot{\gamma} = A \sigma_\mathrm{d}^n f_\mathrm{H2O}^r \exp\left(-\frac{E+PV}{RT}\right)$$

where
- $n$ is the power law exponent
  
- $r$ is the exponent of fugacity dependence
  
- $A$ is a pre-exponential factor $[\mathrm{MPa}^{-n}s^{-1}]$ (if manually defined, $n$ and $r$ must be either pre-defined or substituted)
  
- $E$ is the activation energy $\mathrm{[kJ/mol]}$
  
- $V$ is the activation volume $\mathrm{[m^3/mol]}$
  
- $\dot{\gamma}$ is the strain rate $\mathrm{[1/s]}$
  
- $\sigma_\mathrm{d}$ is the differential stress $\mathrm{[MPa]}$ which are converted into second invariants using the `Apparatus` variable that can be
  

either `AxialCompression`, `SimpleShear` or `Invariant`. If the flow law parameters are already given as a function of second invariants, choose `Apparatus=Invariant`.

**Example**

```julia
julia> x2 = DislocationCreep(n=3)
DislocationCreep: n=3, r=0.0, A=1.5 MPa^-3 s^-1, E=476.0 kJ mol^-1, V=6.0e-6 m^3 mol^-1, Apparatus=AxialCompression
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/DislocationCreep.jl#L19-L43" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.ViscosityPartialMelt_Costa_etal_2009' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.ViscosityPartialMelt_Costa_etal_2009'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.ViscosityPartialMelt_Costa_etal_2009</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ViscosityPartialMelt_Costa_etal_2009(LinearMeltViscosity())
```


The viscosity of a partially molten rock depends on the melt viscosity, melt fraction and strainrate.

This implements a parameterisation of Costa et al. [2009].

**Reference**

Costa, A., Caricchi, L., Bagdassarov, N., 2009. A model for the rheology of particle‐bearing suspensions and partially molten rocks. Geochem Geophys Geosyst 10, 2008GC002138. https://doi.org/10.1029/2008GC002138


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/MeltViscosity.jl#L172-L183" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.LinearMeltViscosity' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.LinearMeltViscosity'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.LinearMeltViscosity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
LinearMeltViscosity(η=1e20Pa*s)
```


Defines a simple temperature-dependent melt viscosity, given by

$$    \tau_{ij} = 2 \eta \dot{\varepsilon}_{ij}$$

or

$$    \dot{\varepsilon}_{ij}  = \frac{\tau_{ij}}{2 \eta}$$

where

$$    \eta =   \eta_0*10^{A + \frac{B}{T - T0}} ;$$

and `\eta_0` is the scaling viscosity, `A` and `B` are constants, and `T_0` is a reference temperature, and `T` is the temperature [in K].

Typical parameters for basalt (default) are: `A = -9.6012`, `B = 1.3374e+04K`, `T_0 = 307.8043K` and `\eta_0 = 1Pas`.

Typical parameters for rhyolite are: `A = -8.1590`, `B = 2.4050e+04K`, `T_0 = -430.9606K` and `\eta_0 = 1Pas`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/MeltViscosity.jl#L15-L35" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.GiordanoMeltViscosity' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.GiordanoMeltViscosity'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.GiordanoMeltViscosity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GiordanoMeltViscosity(; oxd_wt = oxd_wt, η0=1Pas)
```


Defines the melt viscosity model after Giordano et al. (2008) given by

$$    \tau_{ij} = 2 \eta \dot{\varepsilon}_{ij}$$

or

$$    \dot{\varepsilon}_{ij}  = \frac{\tau_{ij}}{2 \eta}$$

where

$$    \eta = \eta_0 * 10^{AT + \frac{BT}{T - CT}}$$

**Parameters**
- `oxd_wt::NTuple{9,T}`: Melt composition as 9-element Tuple containing concentrations           in [wt%] of the following oxides ordered in the exact sequence 
  
  ```
        (SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O) 
  
        Default values are for a hydrous N-MORB melt.
  ```
  
  

`compute_εII`/`compute_τII` (and their derivatives) also accept an `mH2O` keyword (melt water content, mass fraction) that overrides the struct's own `oxd_wt[9]` for that call, recomputing only the water-dependent aggregates. `mH2O` is accepted as given: GeoParams does not check it for consistency with any `Solubility` or other closure's dissolved-water output — that is the caller's/solver's responsibility.

**Reference**
- Giordano D, Russell JK, & Dingwell DB (2008). Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science Letters, 271, 123-134. (https://dx.doi.org/10.1016/j.epsl.2008.03.038)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/MeltViscosity.jl#L333-L367" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.DiffusionCreep' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.DiffusionCreep'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.DiffusionCreep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DiffusionCreep(r = 0NoUnits, p = A = 1.5MPa/s, E = 476.0kJ/mol, V = 6e-6m^3/mol, apparatus = AxialCompression )
```


Defines the flow law parameter of a dislocation creep law.

The (isotropic) diffusion creep law, as used by experimentalists, is given by

$$     \dot{\gamma} = A \sigma_\mathrm{d} d^{\mathrm{p}} f_\mathrm{H2O}^r \exp\left(-\frac{E+PV}{RT}\right)$$

where
- $r$ is the exponent of fugacity dependence
  
- $p$ is the exponent of grain size
  
- $A$ is a pre-exponential factor $[\mathrm{MPa}^{-n}s^{-1}]$ (if manually defined, $n$ and $r$ must be either pre-defined or substituted)
  
- $E$ is the activation energy $\mathrm{[kJ/mol]}$
  
- $V$ is the activation volume $\mathrm{[m^3/mol]}$
  
- $\dot{\gamma}$ is the strain rate $\mathrm{[1/s]}$
  
- $\sigma_\mathrm{d}$ is the differential stress $\mathrm{[MPa]}$
  

The experimental parameters are converted into second invariants using the `Apparatus` variable that can be either `AxialCompression`, `SimpleShear` or `Invariant`. If the flow law parameters are already given as a function of second invariants, choose `Apparatus=Invariant`.

**Example**

```julia
julia> x2 = DiffusionCreep(Name="test")
DiffusionCreep: Name = test, n=1.0, r=0.0, p=-3.0, A=1.5 m³·⁰ MPa⁻¹·⁰ s⁻¹·⁰, E=500.0 kJ mol⁻¹·⁰, V=2.4e-5 m³·⁰ mol⁻¹·⁰, FT=1.7320508075688772, FE=1.1547005383792517)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/DiffusionCreep.jl#L23-L49" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.SetDiffusionCreep' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.SetDiffusionCreep'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.SetDiffusionCreep</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
SetDiffusionCreep["Name of Diffusion Creep"]
```


This is a dictionary with pre-defined creep laws


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/DiffusionCreep.jl#L372-L375" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Computational routines for creep laws {#Computational-routines-for-creep-laws}

Once a creep rheology is defined, we can use the following routines to perform computations within the solvers
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.compute_εII' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.compute_εII'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.compute_εII</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_εII(a::DiffusionCreep, TauII::_T; T::_T, P=one(_T), f=one(_T), d=one(_T), kwargs...)
```


Returns diffusion creep strainrate as a function of 2nd invariant of the stress tensor $\tau_{II}$

$$    \dot{ε}_{II} = A τ_{II}^n d^{p} f_{H_2O}^r \exp \left(- \frac{E + PV}{RT} \right)$$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/DiffusionCreep.jl#L140-L148" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII(a::GrainBoundarySliding, TauII::_T; T::_T, P=one(_T), f=one(_T), d=one(_T), kwargs...)
```


Returns grain boundary sliding strainrate as a function of 2nd invariant of the stress tensor $\tau_{II}$

$$    \dot{ε}_{II} = A τ_{II}^n d^{p} \exp \left(- \frac{E + PV}{RT} \right)$$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/GrainBoundarySliding.jl#L128-L136" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII(s::ConstantElasticity{_T}, τII; τII_old, dt)
```


Computes elastic strainrate given the deviatoric stress at the current (`τII`) and old timestep (`τII_old`), for a timestep `dt`:

$$     \dot{\varepsilon}^{el} = \frac{1}{2 G} \frac{D \tau_{II}}{Dt} \approx \frac{1}{2 G} \frac{\tau_{II} - \tilde{\tau}_{II}^{old}}{dt}$$

Note that we here solve the scalar equation, which is sufficient for isotropic cases. In tensor form, it would be

$$    \dot{\varepsilon}^{el}_{ij} = \frac{1}{2 G} \frac{\tau_{ij} - \tilde{\tau_{ij}}^{old}}{dt}$$

here $\tilde{{\tau_{ij}}}^{old}$ is the rotated old deviatoric stress tensor to ensure objectivity (this can be done with Jaumann derivative, or also by using the full rotational formula).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Elasticity/Elasticity.jl#L100-L114" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII(p::DruckerPrager{_T,U,U1}, λdot::_T, τII::_T,  P)
```


This computes plastic strain rate invariant for a given $λdot$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Plasticity/DruckerPrager.jl#L179-L183" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII(p::DruckerPrager_regularised{_T,U,U1}, λdot::_T, τII::_T,  P)
```


This computes plastic strain rate invariant for a given $λdot$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Plasticity/DruckerPrager_regularised.jl#L185-L189" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII(v::Parallel{T,N}, τII, args; tol=1e-6, verbose=false, n=1)
```


Computing `εII` as a function of `τII` for a Parallel elements is (usually) a nonlinear problem


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CompositeRheologies/Parallel.jl#L47-L51" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII(v::CompositeRheology{T,N}, τII, args; tol=1e-6, verbose=false, n=1)
```


Computing `εII` as a function of `τII` for a composite element is the sum of the individual contributions


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CompositeRheologies/CompositeRheology.jl#L77-L81" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII(v::AbstractPlasticity, τII::_T, args; tol=1e-6, verbose=true)
```


Performs local iterations to compute the plastic strainrate. Note that the non-plastic strainrate, ε_np, should be part of `args`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CompositeRheologies/NonlinearIterations.jl#L133-L137" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.compute_εII!' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.compute_εII!'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.compute_εII!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_εII!(EpsII::AbstractArray{_T,N}, a, TauII::AbstractArray{_T,N}; T, P, f,d,kwargs...)
```


Computes strainrate as a function of stress


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/DiffusionCreep.jl#L178-L182" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII!(EpsII::AbstractArray{_T,N}, a, TauII::AbstractArray{_T,N}; T, P, f,d,kwargs...)
```


Computes strainrate as a function of stress


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/GrainBoundarySliding.jl#L168-L172" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII!(EpsII::AbstractArray{_T,N}, s::LinearMeltViscosity, TauII::AbstractArray{_T,N}; T, kwargs...)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/MeltViscosity.jl#L71-L73" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII!(EpsII::AbstractArray{_T,N}, s::ViscosityPartialMelt_Costa_etal_2009, TauII::AbstractArray{_T,N}; T, kwargs...)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/MeltViscosity.jl#L237-L239" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII!(EpsII::AbstractArray{_T,N}, s::GiordanoMeltViscosity, TauII::AbstractArray{_T,N}; T, kwargs...)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/MeltViscosity.jl#L498-L500" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII!(EpsII::AbstractArray{_T, N}, a::HerschelBulkley, TauII::AbstractArray; T = one(_T), kwargs...)
```


In-place function for the second invariant of the strain rate for Herschel-Bulkley rheology.

`T` may be a scalar, applied to every element, or an array indexed alongside `TauII`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/HerschelBulkley.jl#L58-L64" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII!(EpsII::AbstractArray, s::LinearViscous, TauII::AbstractArray)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/CreepLaw.jl#L144-L147" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII!(EpsII::AbstractArray{_T,N}, s::ArrheniusType, TauII::AbstractArray{_T,N})
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/CreepLaw.jl#L268-L271" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_εII!(ε_el::AbstractArray{_T,N}, s::ConstantElasticity{_T}; τII::AbstractArray{_T,N}, τII_old::AbstractArray{_T,N}, dt::_T, kwargs...)
```


In-place computation of the elastic shear strainrate for given deviatoric stress invariants at the previous (`τII_old`) and new (`τII`) timestep, as well as the timestep `dt`

$$     \dot{\varepsilon}^{el} = \frac{1}{2 G} \frac{D \tau_{II}}{Dt} \approx \frac{1}{2 G} \frac{\tau_{II} - \tau_{II}^{old}}{dt}$$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Elasticity/Elasticity.jl#L153-L162" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.compute_τII' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.compute_τII'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.compute_τII</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_τII(a::DislocationCreep, EpsII; P, T, f, args...)
```


Computes the stress for a Dislocation creep law given a certain strain rate


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/DislocationCreep.jl#L208-L213" target="_blank" rel="noreferrer">source</a></Badge>



```julia
computeCreepLaw_TauII(EpsII::_T, a::DiffusionCreep; T::_T, P=zero(_T), f=one(_T), d=one(_T), kwargs...)
```


Returns diffusion creep stress as a function of 2nd invariant of the strain rate


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/DiffusionCreep.jl#L243-L247" target="_blank" rel="noreferrer">source</a></Badge>



```julia
computeCreepLaw_TauII(EpsII::_T, a::GrainBoundarySliding; T::_T, P=zero(_T), d=one(_T), kwargs...)
```


Returns grain boundary sliding stress as a function of 2nd invariant of the strain rate


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/GrainBoundarySliding.jl#L224-L228" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII(a::PeierlsCreep, EpsII; P, T, f, args...)
```


Computes the stress for a peierls creep law given a certain strain rate.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/PeierlsCreep.jl#L164-L169" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII(s::LinearMeltViscosity, EpsII; kwargs...)
```


Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/MeltViscosity.jl#L106-L110" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII(s::ViscosityPartialMelt_Costa_etal_2009, EpsII; kwargs...)
```


Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/MeltViscosity.jl#L260-L264" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII(s::GiordanoMeltViscosity, EpsII; kwargs...)
```


Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/MeltViscosity.jl#L524-L528" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII(a::HerschelBulkley, EpsII; T = one(precision(a)), kwargs...)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/HerschelBulkley.jl#L79-L82" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII(s::LinearViscous, EpsII; kwargs...)
```


Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/CreepLaw.jl#L171-L175" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII(s::ArrheniusType, EpsII; kwargs...)
```


Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/CreepLaw.jl#L294-L298" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII(s::PowerlawViscous, EpsII; kwargs...)
```


Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/CreepLaw.jl#L380-L384" target="_blank" rel="noreferrer">source</a></Badge>



```julia
τII = compute_τII(v::CompositeRheology{T,N}, εII, args; tol=1e-6, verbose=false)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CompositeRheologies/CompositeRheology.jl#L215-L218" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.ConstitutiveRelationships.compute_τII!' href='#GeoParams.MaterialParameters.ConstitutiveRelationships.compute_τII!'><span class="jlbinding">GeoParams.MaterialParameters.ConstitutiveRelationships.compute_τII!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_τII!(TauII::AbstractArray{_T,N}, a::DislocationCreep, EpsII::AbstractArray;
    P = zero(_T),
    T = one(_T),
    f = one(_T))
```


Computes the deviatoric stress invariant for a dislocation creep law.

Each of `P`, `T`, and `f` may be a scalar, applied to every element, or an array indexed alongside `EpsII`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/DislocationCreep.jl#L247-L257" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII!(TauII::AbstractArray{_T,N}, a::PeierlsCreep, EpsII::AbstractArray;
    T = one(_T), args...)
```


Computes the deviatoric stress invariant for a peierls creep law.

`T` may be a scalar, applied to every element, or an array indexed alongside `EpsII`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/PeierlsCreep.jl#L201-L209" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII!(TauII::AbstractArray{_T, N}, a::HerschelBulkley, EpsII::AbstractArray; T = one(_T), kwargs...)
```


In-place function for the second invariant of the stress for Herschel-Bulkley rheology.

`T` may be a scalar, applied to every element, or an array indexed alongside `EpsII`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/CreepLaw/HerschelBulkley.jl#L95-L101" target="_blank" rel="noreferrer">source</a></Badge>



```julia
compute_τII!(τII::AbstractArray{_T,N}, s::ConstantElasticity{_T}. ε_el::AbstractArray{_T,N}; τII_old::AbstractArray{_T,N}, dt::_T, kwargs...)
```


In-place update of the elastic stress for given deviatoric strainrate invariants and stres invariant at the old (`τII_old`) timestep, as well as the timestep `dt`

$$    \tau_{II} = 2 G dt \dot{\varepsilon}^{el} + \tau_{II}^{old}$$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/fe1e2fd766bac30f9615d0d78f80f1e017d794a0/src/Elasticity/Elasticity.jl#L177-L186" target="_blank" rel="noreferrer">source</a></Badge>

</details>

