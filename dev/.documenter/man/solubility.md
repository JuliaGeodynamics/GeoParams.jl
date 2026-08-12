---
---

# Solubility {#Solubility}

## Methods {#Methods}

Coupled H2O–CO2 volatile solubility (dissolved-content) closures can be set with:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.Liu2005_Solubility' href='#GeoParams.MaterialParameters.Solubility.Liu2005_Solubility'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.Liu2005_Solubility</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Liu2005_Solubility(; coeffs=(b1..b6, c1..c4), Pref=1MPa, Tref=1K)
```


Coupled H2O–CO2 solubility for silicic (rhyolite) melt after Liu et al. (2005), as used by Degruyter & Huber (2014). [`compute_dissolved`](/man/solubility#GeoParams.MaterialParameters.Solubility.compute_dissolved) returns the dissolved H2O and CO2 mass fractions of the melt as a function of pressure, temperature, and the CO2 mole fraction of the gas `X_co2`:

$$    m_{H_2O} = (b_1 P_w^{1/2} + b_2 P_w + b_3 P_w^{3/2})\frac{T_{ref}}{T}
             + b_4 P_w^{3/2} + P_c(b_5 P_w^{1/2} + b_6 P_w)$$

$$    m_{CO_2} = P_c(c_1 + c_2 P_w)\frac{T_{ref}}{T} + P_c(c_3 P_w^{1/2} + c_4 P_w^{3/2})$$

with the dimensionless partial pressures $P_w = P(1-X_{co2})/P_{ref}$ and $P_c = P X_{co2}/P_{ref}$. With `Pref=1MPa`, `Tref=1K` these equal the reference partial pressures in MPa and reproduce the Liu (2005) numbers; because they are ratios of like-dimensioned quantities, the closure is dimensionally homogeneous and nondimensionalizes cleanly. The fitted `coeffs` are dimensionless; `isdimensional` tracks the reference `GeoUnit`s.

**References**
- Liu, Y., Zhang, Y., Behrens, H. (2005), Solubility of H2O in rhyolitic melts at low pressures and a new empirical model for mixed H2O-CO2 solubility, JVGR 143, 219-235, https://doi.org/10.1016/j.jvolgeores.2004.09.019
  
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L40-L64" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.Mafic_Solubility' href='#GeoParams.MaterialParameters.Solubility.Mafic_Solubility'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.Mafic_Solubility</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Mafic_Solubility(; coeffs=(b1..b10, c1..c4), T0=273.15K, Tref=1K, Pref=1MPa)
```


Coupled H2O–CO2 solubility for mafic (basalt) melt. Dissolved H2O follows the mafic polynomial of the Scholz/Degruyter–Huber reference in the dimensionless groups $T_C = (T-T_0)/T_{ref}$ (numerically °C) and $P_m = P/P_{ref}$ (numerically MPa); dissolved CO2 reuses the Liu (2005) rhyolite CO2 block. Nondimensionalizes like [`Liu2005_Solubility`](/man/solubility#GeoParams.MaterialParameters.Solubility.Liu2005_Solubility). The H2O polynomial has no floor at zero and goes negative outside its calibration (high `X_co2`, low `P`); the dissolved H2O output is floored at zero rather than returned as a negative mass fraction.

**References**
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
  
- Liu, Y., Zhang, Y., Behrens, H. (2005), Solubility of H2O in rhyolitic melts at low pressures and a new empirical model for mixed H2O-CO2 solubility, JVGR 143, 219-235, https://doi.org/10.1016/j.jvolgeores.2004.09.019
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L124-L139" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.GasMixture' href='#GeoParams.MaterialParameters.Solubility.GasMixture'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.GasMixture</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GasMixture(; Cp_h2o=3880J/kg/K, Cp_co2=1200J/kg/K, M_h2o=18.02e-3kg/mol, M_co2=44.01e-3kg/mol)
```


H2O–CO2 gas-mixture properties keyed on the CO2 mole fraction of the gas `X_co2` (Degruyter & Huber 2014). The effective molar mass is

$$    m_g = M_{H_2O}(1-X_{co2}) + M_{CO_2} X_{co2}$$

and the mass-weighted specific heat is

$$    c_g = \frac{M_{H_2O} c_{H_2O}(1-X_{co2}) + M_{CO_2} c_{CO_2} X_{co2}}{m_g}$$

with the reference convention $c_g = 0$ at $X_{co2}=0$ — a discontinuity, since the formula's own limit as $X_{co2} \to 0$ is $c_{H_2O}$.

All fields are `GeoUnit`s, so the struct nondimensionalizes, but both accessors read values rather than `Quantity`s: neither return carries units.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L207-L224" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Computational routines {#Computational-routines}

To evaluate dissolved H2O and CO2 (mass fractions) within a user routine, use this:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.compute_dissolved' href='#GeoParams.MaterialParameters.Solubility.compute_dissolved'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.compute_dissolved</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_dissolved(s::AbstractSolubility, P, T, X_co2) -> (m_h2o, m_co2)
```


Dissolved H2O and CO2 mass fractions of the melt (`P` in Pa, `T` in K, `X_co2` the CO2 mole fraction of the gas). Also callable as `compute_dissolved(s; P, T, X_co2)` and `compute_dissolved(s, args::NamedTuple)`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L92-L98" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.compute_dissolved!' href='#GeoParams.MaterialParameters.Solubility.compute_dissolved!'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.compute_dissolved!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_dissolved!(m_h2o, m_co2, MatParam, Phases, args)
```


In-place dissolved H2O and CO2 over a domain. `args` is a NamedTuple of `P, T, X_co2` index-matched to the output arrays. Also accepts a single solubility or phase struct in place of `(MatParam, Phases)`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L340-L346" target="_blank" rel="noreferrer">source</a></Badge>

</details>


Partial derivatives (via `ForwardDiff`, useful for implicit/Newton solves) are available for pressure, temperature, and gas composition:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.∂dissolved_∂P' href='#GeoParams.MaterialParameters.Solubility.∂dissolved_∂P'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.∂dissolved_∂P</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
∂dissolved_∂P(s, P, T, X_co2) -> (∂m_h2o/∂P, ∂m_co2/∂P)
```


ForwardDiff partial derivatives of [`compute_dissolved`](/man/solubility#GeoParams.MaterialParameters.Solubility.compute_dissolved) with respect to pressure. Companions [`∂dissolved_∂T`](/man/solubility#GeoParams.MaterialParameters.Solubility.∂dissolved_∂T), [`∂dissolved_∂Xco2`](/man/solubility#GeoParams.MaterialParameters.Solubility.∂dissolved_∂Xco2).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L277-L282" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.∂dissolved_∂T' href='#GeoParams.MaterialParameters.Solubility.∂dissolved_∂T'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.∂dissolved_∂T</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
∂dissolved_∂T(s, P, T, X_co2) -> (∂m_h2o/∂T, ∂m_co2/∂T)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L285-L287" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.∂dissolved_∂Xco2' href='#GeoParams.MaterialParameters.Solubility.∂dissolved_∂Xco2'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.∂dissolved_∂Xco2</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
∂dissolved_∂Xco2(s, P, T, X_co2) -> (∂m_h2o/∂X_co2, ∂m_co2/∂X_co2)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L290-L292" target="_blank" rel="noreferrer">source</a></Badge>

</details>


Given a target dissolved H2O content, the gas composition `X_co2` that produces it can be found with:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.find_Xco2' href='#GeoParams.MaterialParameters.Solubility.find_Xco2'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.find_Xco2</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
find_Xco2(s::AbstractSolubility, P, T, m_h2o_target; X0=0.5, tol=1e-8, max_iter=50) -> X_co2
```


Invert [`compute_dissolved`](/man/solubility#GeoParams.MaterialParameters.Solubility.compute_dissolved) for the gas composition: solve for `X_co2 ∈ [0,1]` such that `compute_dissolved(s, P, T, X_co2)[1] == m_h2o_target`, at fixed pressure `P` and temperature `T`. Uses a safeguarded Newton iteration (via [`∂dissolved_∂Xco2`](/man/solubility#GeoParams.MaterialParameters.Solubility.∂dissolved_∂Xco2)) bracketed by bisection, so an out-of-bounds Newton step always falls back to a bisection halving instead of leaving `[0,1]`.

Throws an `ErrorException` if `m_h2o_target` is infeasible at this `P,T` (outside the achievable range between `X_co2=0` and `X_co2=1`) or if the iteration fails to converge within `max_iter` steps — this never returns a clamped or out-of-tolerance answer.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L295-L309" target="_blank" rel="noreferrer">source</a></Badge>

</details>


Gas-mixture helpers for the energy equation:
<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.compute_gas_heatcapacity' href='#GeoParams.MaterialParameters.Solubility.compute_gas_heatcapacity'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.compute_gas_heatcapacity</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
compute_gas_heatcapacity(s::GasMixture, X_co2)
```


Mass-weighted specific heat of the H2O–CO2 gas mixture; zero at `X_co2 == 0` (reference convention).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L248-L253" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeoParams.MaterialParameters.Solubility.effective_molar_mass' href='#GeoParams.MaterialParameters.Solubility.effective_molar_mass'><span class="jlbinding">GeoParams.MaterialParameters.Solubility.effective_molar_mass</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
effective_molar_mass(s::GasMixture, X_co2)
```


Effective molar mass of the H2O–CO2 gas mixture at CO2 mole fraction `X_co2`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaGeodynamics/GeoParams.jl/blob/68ac1af6736097999812363ccbc1c6ec744384f3/src/Solubility/Solubility.jl#L238-L242" target="_blank" rel="noreferrer">source</a></Badge>

</details>


Note that `ε_g` (exsolved-gas volume fraction) and `X_co2` (CO2 mole fraction of the gas) are typically ODE state variables owned by the solver, not closed algebraically at a point — the closures above are the _entries_ to that mass-balance, not a replacement for it.
