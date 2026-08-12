---
---

# MaterialParameters {#MaterialParameters}

Material properties for a given phase can be set with `SetMaterialParams`, whereas all properties are  stored in the `MaterialParams` structure. Information about the material parameter is found at `MaterialParamsInfo`. This can be employed 
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.MaterialParamsInfo' href='#GeoParams.MaterialParameters.MaterialParamsInfo'><span class="jlbinding">GeoParams.MaterialParameters.MaterialParamsInfo</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MaterialParamsInfo
```


Structure that holds information (Equation, Comment, BibTex_Reference) about a given material parameter, which can be used to create parameter tables, documentation etc.

Usually used in combination with `param_info(the_parameter_of_interest)`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/MaterialParameters.jl#L26-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.SetMaterialParams' href='#GeoParams.MaterialParameters.SetMaterialParams'><span class="jlbinding">GeoParams.MaterialParameters.SetMaterialParams</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
SetMaterialParams(; Name::String="", Phase::Int64=1,
                    Density             =   nothing,
                    Gravity             =   nothing,
                    CreepLaws           =   nothing,
                    Elasticity          =   nothing,
                    Plasticity          =   nothing,
                    CompositeRheology   =   nothing,
                    ChemDiffusion       =   nothing,
                    Conductivity        =   nothing,
                    HeatCapacity        =   nothing,
                    RadioactiveHeat     =   nothing,
                    LatentHeat          =   nothing,
                    ShearHeat           =   nothing,
                    Permeability        =   nothing,
                    Melting             =   nothing,
                    SeismicVelocity     =   nothing,
                    Solubility          =   nothing,
                    CharDim::GeoUnits   =   nothing)
```


Sets material parameters for a given phase.

If `CharDim` is specified the input parameters are non-dimensionalized. Note that if `Density` is specified, we also set `Gravity` even if not explicitly listed

**Examples**

Define two viscous creep laws & constant density:

```julia
julia> Phase = SetMaterialParams(Name="Viscous Matrix",
                       Density   = ConstantDensity(),
                       CreepLaws = (PowerlawViscous(), LinearViscous(η=1e21Pa*s)))
Phase 1 : Viscous Matrix
        | [dimensional units]
        |
        |-- Density           : Constant density: ρ=2900 kg m⁻³
        |-- Gravity           : Gravitational acceleration: g=9.81 m s⁻²
        |-- CreepLaws         : Powerlaw viscosity: η0=1.0e18 Pa s, n=2.0, ε0=1.0e-15 s⁻¹
        |                       Linear viscosity: η=1.0e21 Pa s
```


Define two viscous creep laws & P/T dependent density and nondimensionalize

```julia
julia> CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
julia> Phase = SetMaterialParams(Name="Viscous Matrix", Phase=33,
                              Density   = PT_Density(),
                              CreepLaws = (PowerlawViscous(n=3), LinearViscous(η=1e23Pa*s)),
                              CharDim   = CharUnits_GEO)
Phase 33: Viscous Matrix
        | [non-dimensional units]
        |
        |-- Density           : P/T-dependent density: ρ0=2.9e-16, α=0.038194500000000006, β=0.01, T0=0.21454659702313156, P0=0.0
        |-- Gravity           : Gravitational acceleration: g=9.810000000000002e18
        |-- CreepLaws         : Powerlaw viscosity: η0=0.1, n=3, ε0=0.001
        |                       Linear viscosity: η=10000.0
```


You can also create an array that holds several parameters:

```julia
julia> MatParam        =   Array{MaterialParams, 1}(undef, 2);
julia> Phase           =   1;
julia> MatParam[Phase] =   SetMaterialParams(Name="Upper Crust", Phase=Phase,
                            CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                            Density  = ConstantDensity(ρ=2900kg/m^3));
julia> Phase           =   2;
julia> MatParam[Phase] =   SetMaterialParams(Name="Lower Crust", Phase=Phase,
                            CreepLaws= (PowerlawViscous(n=5), LinearViscous(η=1e21Pa*s)),
                            Density  = PT_Density(ρ0=3000kg/m^3));
julia> MatParam
2-element Vector{MaterialParams}:
 Phase 1 : Upper Crust
    | [dimensional units]
    |
    |-- Density           : Constant density: ρ=2900 kg m⁻³
    |-- Gravity           : Gravitational acceleration: g=9.81 m s⁻²
    |-- CreepLaws         : Powerlaw viscosity: η0=1.0e18 Pa s, n=2.0, ε0=1.0e-15 s⁻¹
    |                       Linear viscosity: η=1.0e23 Pa s

 Phase 2 : Lower Crust
    | [dimensional units]
    |
    |-- Density           : P/T-dependent density: ρ0=3000 kg m⁻³, α=3.0e-5 K⁻¹, β=1.0e-9 Pa⁻¹, T0=0 °C, P0=0 MPa
    |-- Gravity           : Gravitational acceleration: g=9.81 m s⁻²
    |-- CreepLaws         : Powerlaw viscosity: η0=1.0e18 Pa s, n=5, ε0=1.0e-15 s⁻¹
    |                       Linear viscosity: η=1.0e21 Pa s
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/MaterialParameters.jl#L111-L199" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.MaterialParams' href='#GeoParams.MaterialParameters.MaterialParams'><span class="jlbinding">GeoParams.MaterialParameters.MaterialParams</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MaterialParams
```


Structure that holds all material parameters for a given phase


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/MaterialParameters.jl#L64-L69" target="_blank" rel="noreferrer">source</a></Badge>

</details>

