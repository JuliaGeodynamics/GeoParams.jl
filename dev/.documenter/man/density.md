---
---

# Density {#Density}

## Methods {#Methods}

The density equation of state can be specified in a number of ways
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.ConstantDensity' href='#GeoParams.MaterialParameters.Density.ConstantDensity'><span class="jlbinding">GeoParams.MaterialParameters.Density.ConstantDensity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ConstantDensity(ρ=2900kg/m^3)
```


Set a constant density:

$$    \rho  = cst$$

where $\rho$ is the density [$kg/m^3$].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L54-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.PT_Density' href='#GeoParams.MaterialParameters.Density.PT_Density'><span class="jlbinding">GeoParams.MaterialParameters.Density.PT_Density</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PT_Density(ρ0=2900kg/m^3, α=3e-5/K, β=1e-9/Pa, T0=0C, P=0MPa)
```


Set a pressure and temperature-dependent density:

$$    \rho  = \rho_0 (1.0 - \alpha (T-T_0) + \beta  (P-P_0) )$$

where $\rho_0$ is the density [$kg/m^3$] at reference temperature $T_0$ and pressure $P_0$, $\alpha$ is the temperature dependence of density and $\beta$ the pressure dependence.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L98-L107" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.T_Density' href='#GeoParams.MaterialParameters.Density.T_Density'><span class="jlbinding">GeoParams.MaterialParameters.Density.T_Density</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
T_Density(ρ0=2900kg/m^3, α=3e-5/K, T₀=273.15K)
```


Set a temperature-dependent density:

$$    \rho  = \rho_0 (1 - \alpha * (T - T\_0) )$$

where $\rho_0$ is the density [$kg/m^3$] at reference temperature $T_0$ and $\alpha$ the temperature dependence.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L192-L200" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.Compressible_Density' href='#GeoParams.MaterialParameters.Density.Compressible_Density'><span class="jlbinding">GeoParams.MaterialParameters.Density.Compressible_Density</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Compressible_Density(ρ0=2900kg/m^3, β=1e-9/Pa, P₀=0MPa)
```


Set a pressure-dependent density:

$$    \rho  = \rho_0 \exp(β*(P - P\_0))$$

where $\rho_0$ is the density [$kg/m^3$] at reference pressure $P_0$ and $\beta$ the pressure dependence.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L148-L156" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.MeltDependent_Density' href='#GeoParams.MaterialParameters.Density.MeltDependent_Density'><span class="jlbinding">GeoParams.MaterialParameters.Density.MeltDependent_Density</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MeltDependent_Density(ρsolid=ConstantDensity(), ρmelt=ConstantDensity())
```


If we use a single phase code the average density of a partially molten rock is

$$    \rho  = \phi \rho_{\textrm{melt}} + (1-\phi) \rho_{\textrm{solid}}$$

where $\rho$ is the average density [$kg/m^3$], $\rho_{\textrm{melt}}$ the melt density, $\rho_{\textrm{solid}}$ the solid density and $\phi$ the melt fraction.

Note that any density formulation can be used for melt and solid.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L236-L246" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.BubbleFlow_Density' href='#GeoParams.MaterialParameters.Density.BubbleFlow_Density'><span class="jlbinding">GeoParams.MaterialParameters.Density.BubbleFlow_Density</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BubbleFlow_Density(ρmelt=ConstantDensity(), ρgas=ConstantDensity(), c0=0e0, a=0.0041MPa^-1/2)
```


Defines the BubbleFlow_Density as described in Slezin (2003) with a default gas solubility constant of 0.0041MPa$^{-1/2}$ used in e.g. Sparks et al. (1978)

$$    \rho = \frac{1}{\frac{c_0 - c}{\rho_g} + \frac{1-(c_0-c)}{\rho_m}}$$

with

$$c =
\begin{cases}
   aP^{1/2} & \text{for } P < \frac{c_0^2}{a^2} \\
    c_0 & \text{for } P \geq \frac{c_0^2}{a^2}
\end{cases}$$

**Arguments**
- `ρmelt`: Density of the melt
  
- `ρgas`: Density of the gas
  
- `c0`: Total volatile content
  
- `a`: Gas solubility constant (default: 4.1e-6Pa$^{-1/2}$) (after Sparks et al., 1978)
  

Possible values for a are 3.2e-6-6.4e-6Pa$^{-1/2}$ where the lower value corresponds to mafic magmas at rather large pressures (400-600MPa) and the higher value to felsic magmas at low pressures (0 to 100-200MPa) (after Slezin (2003))

**Example**

```julia
rheology = SetMaterialParams(;
                      Phase=1,
                      CreepLaws=(PowerlawViscous(), LinearViscous(; η=1e21Pa * s)),
                      Gravity=ConstantGravity(; g=9.81.0m / s^2),
                      Density= BubbleFlow_Density(ρmelt=ConstantDensity(ρ=2900kg/m^3), ρgas=ConstantDensity(ρ=1kg/m^3), c0=0.0, a=0.0041MPa^-1//2),
                      )
```


**References**
- Slezin, Yu. B. (2003), The mechanism of volcanic eruptions (a steady state approach), Journal of Volcanology and Geothermal Research, 122, 7-50, https://doi.org/10.1016/S0377-0273(02)00464-X
  
- Sparks, R. S. J.(1978), The dynamics of bubble formation and growth in magmas: A review and analysis, Journal of Volcanology and Geothermal Research, 3, 1-37, https://doi.org/10.1016/0377-0273(78)90002-1
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L279-L315" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.GasPyroclast_Density' href='#GeoParams.MaterialParameters.Density.GasPyroclast_Density'><span class="jlbinding">GeoParams.MaterialParameters.Density.GasPyroclast_Density</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GasPyroclast_Density(ρmelt=ConstantDensity(), ρgas=ConstantDensity(), δ=0e0)
```


Defines the GasPyroclast_Density as described in Slezin (2003) with a default volume fraction of free gas in the flow of 0.0 This is also used to model partly destroyed foam in the conduit.

$$    \rho = \rho_g\delta + \rho_p(1 - \delta)$$

with

$$    \rho_p = \rho_m(1 - \beta) + \rho_g\beta \approx \rho_l(1 - \beta)$$

**Arguments**
- `ρmelt`: Density of the melt
  
- `ρgas`: Density of the gas
  
- `δ`: Volume fraction of free gas in the flow
  
- `β`: Gas volume fraction enclosed within the particles
  

**Example**

```julia
rheology = SetMaterialParams(;
                      Phase=1,
                      CreepLaws=(PowerlawViscous(), LinearViscous(; η=1e21Pa * s)),
                      Gravity=ConstantGravity(; g=9.81.0m / s^2),
                      Density= GasPyroclast_Density(ρmelt=ConstantDensity(ρ=2900kg/m^3), ρgas=ConstantDensity(ρ=1kg/m^3), δ=0.0, β=0.0),
                      )
```


**References**
- Slezin, Yu. B. (2003), The mechanism of volcanic eruptions (a steady state approach), Journal of Volcanology and Geothermal Research, 122, 7-50, https://doi.org/10.1016/S0377-0273(02)00464-X
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L361-L393" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.Melt_DensityX' href='#GeoParams.MaterialParameters.Density.Melt_DensityX'><span class="jlbinding">GeoParams.MaterialParameters.Density.Melt_DensityX</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Melt_DensityX(; oxd_wt = oxd_wt)
```


Set a density depending on the oxide composition after the python script by Iacovino K & Till C (2019)

**Arguments**
- `oxd_wt::NTuple{9, T}`: Melt composition as 9-element Tuple containing concentrations            in [wt%] of the following oxides ordered in the exact sequence 
  
  ```
         (SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O) 
  
         Default values are for a hydrous N-MORB composition
  ```
  
  

The callable also accepts an `mH2O` keyword (melt water content, mass fraction) that overrides the struct's own `oxd_wt[9]` for that call, recomputing only the water-dependent aggregates. `mH2O` is accepted as given: GeoParams does not check it for consistency with any `Solubility` or other closure's dissolved-water output — that is the caller's/solver's responsibility.

**References**
- Iacovino K & Till C (2019). DensityX: A program for calculating the densities of magmatic liquids up to 1,627 °C and 30 kbar. Volcanica 2(1), p 1-10. [doi:10.30909/vol.02.01.0110](https://dx.doi.org/10.30909/vol.02.01.0110)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L603-L625" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.RedlichKwong_Density' href='#GeoParams.MaterialParameters.Density.RedlichKwong_Density'><span class="jlbinding">GeoParams.MaterialParameters.Density.RedlichKwong_Density</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
RedlichKwong_Density(; coeffs=(-112.528, 127.811, 112.04), T0=273.15K, Tref=1K, Pref=1e5Pa, ρref=1e3kg/m^3)
```


Modified Redlich–Kwong gas density of Huber et al. (2010), as used by Degruyter & Huber (2014). Fitted for H2O gas over roughly $873 < T < 1173$ K and $30 < P < 400$ MPa:

$$    \rho_g = \rho_{ref}\left(a_1\,\tau^{-0.381} + a_2\,\varpi^{-1.135}
             + a_3\,\tau^{-0.411}\varpi^{0.033}\right),
    \quad \tau = \frac{T - T_0}{T_{ref}},\ \varpi = \frac{P}{P_{ref}}$$

The fit's specific units (°C, bar) enter through the reference quantities: with `T0=273.15K`, `Tref=1K` the group $\tau$ equals $T$ in °C, and with `Pref=1e5Pa` (1 bar) $\varpi$ equals $P$ in bar. Because $\tau$ and $\varpi$ are ratios of like-dimensioned quantities, the expression is dimensionally homogeneous and nondimensionalizes cleanly: only the reference `GeoUnit`s scale, the fitted dimensionless `coeffs` do not. The parameterisation ignores gas composition (pseudo-pure H2O), matching the reference.

Returns `NaN` for `T ≤ 273.15K` or `P ≤ 0` (`τ ≤ 0` or `ω ≤ 0`), since `τ^(-0.381)`/`τ^(-0.411)` and `ω^(-1.135)` require a positive base: raising to those exponents would otherwise throw an opaque low-level `DomainError` about complex exponentiation.

**References**
- Huber C., Bachmann O. , Manga M., Two Competing Effects of Volatiles on Heat Transfer in Crystal-rich Magmas: Thermal Insulation vs Defrosting, Journal of Petrology, 847–867, https://doi.org/10.1093/petrology/egq003
  
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L428-L455" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.IdealGas_Density' href='#GeoParams.MaterialParameters.Density.IdealGas_Density'><span class="jlbinding">GeoParams.MaterialParameters.Density.IdealGas_Density</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
IdealGas_Density(; Rs=461.5J/kg/K)
```


Ideal-gas density

$$    \rho_g = \frac{P}{R_s T}$$

with specific gas constant `Rs` (default: water vapor, $R/M_w = 8.314/0.01802$). Unlike [`RedlichKwong_Density`](/man/density#GeoParams.MaterialParameters.Density.RedlichKwong_Density) this is dimensionally homogeneous, so it nondimensionalizes cleanly and serves as the analytic derivative check for the Redlich–Kwong fit.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L507-L518" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.ThreePhase_Density' href='#GeoParams.MaterialParameters.Density.ThreePhase_Density'><span class="jlbinding">GeoParams.MaterialParameters.Density.ThreePhase_Density</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ThreePhase_Density(; ρmelt=ConstantDensity(ρ=2300kg/m^3), ρsolid=ConstantDensity(ρ=2600kg/m^3), ρgas=RedlichKwong_Density(), ρ=2400kg/m^3)
```


Volume-fraction mixture density of a melt + crystal + gas magma (Degruyter & Huber 2014):

$$    \rho = \varepsilon_m \rho_m + \varepsilon_g \rho_g + \varepsilon_x \rho_x,
    \qquad \varepsilon_m = 1 - \varepsilon_x - \varepsilon_g$$

The gas and crystal volume fractions `ϕ_gas`, `ϕ_x` are inputs (from the solver's ODE state / a melting law), passed through `args`; the sub-densities receive the same `args` (e.g. `P`, `T`) so an equation-of-state gas density is evaluated in place. `ρ` is a dimensional-tracking sentinel.

`ρgas` is only evaluated when `ϕ_gas != 0`: a cell with no exsolved gas never depends on the gas closure at all, so a strict EOS like [`RedlichKwong_Density`](/man/density#GeoParams.MaterialParameters.Density.RedlichKwong_Density) (which returns `NaN` outside its calibration window) cannot poison a `ϕ_gas=0` cell regardless of `P`, `T`.

**Arguments**
- `ρmelt`: melt-phase density
  
- `ρsolid`: crystal-phase density
  
- `ρgas`: gas-phase density (e.g. [`RedlichKwong_Density`](/man/density#GeoParams.MaterialParameters.Density.RedlichKwong_Density))
  

**References**
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L545-L571" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Computational routines {#Computational-routines}

To evaluate density within a user routine, use this:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.compute_density' href='#GeoParams.MaterialParameters.Density.compute_density'><span class="jlbinding">GeoParams.MaterialParameters.Density.compute_density</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_density(P,T, s::AbstractPhaseDiagramsStruct)
```


Interpolates density as a function of `T,P` from a lookup table


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L813-L816" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Density.compute_density!' href='#GeoParams.MaterialParameters.Density.compute_density!'><span class="jlbinding">GeoParams.MaterialParameters.Density.compute_density!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_density!(rho::AbstractArray{_T, ndim}, MatParam::NTuple{N,AbstractMaterialParamsStruct}, Phases::AbstractArray{_I, ndim}; P=nothing, T=nothing) where {ndim,N,_T,_I<:Integer}
```


In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays. This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.

**Example**

```julia
julia> MatParam = (SetMaterialParams(Name="Mantle", Phase=1,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PT_Density()
                        ),
                    SetMaterialParams(Name="Crust", Phase=2,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pas)),
                        Density   = ConstantDensity(ρ=2900kg/m^3))
                  );
julia> Phases = ones(Int64,10,10);
julia> Phases[:,5:end] .= 2
julia> rho     = zeros(size(Phases))
julia> T       =  ones(size(Phases))
julia> P       =  ones(size(Phases))*10
julia> args = (P=P, T=T)
julia> compute_density!(rho, MatParam, Phases, args)
julia> rho
10×10 Matrix{Float64}:
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
```


The routine is made to minimize allocations:

```julia
julia> using BenchmarkTools
julia> @btime compute_density!($rho, $MatParam, $Phases, P=$P, T=$T)
    203.468 μs (0 allocations: 0 bytes)
```


_________________________________________________________________________________________________________

```
compute_density!(rho::AbstractArray{_T, N}, MatParam::NTuple{K,AbstractMaterialParamsStruct}, PhaseRatios::AbstractArray{_T, M}, P=nothing, T=nothing)
```


In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays. This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Density/Density.jl#L838-L885" target="_blank" rel="noreferrer">source</a></Badge>

</details>


Note that density values are usually not used in itself in the governing PDE's, but usually in combination with other parameters, such as $\rho g$ or $\rho c_p$. The non-dimensional value of $\rho$ may thus have very large or small values, but multiplied with the other values one often obtains numbers that are closer to one.
