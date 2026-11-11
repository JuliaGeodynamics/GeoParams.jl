# Solubility

## Methods
Coupled H2O–CO2 volatile solubility (dissolved-content) closures can be set with:
```@docs
GeoParams.MaterialParameters.Solubility.Liu2005_Solubility
GeoParams.MaterialParameters.Solubility.Mafic_Solubility
GeoParams.MaterialParameters.Solubility.SunYao2026_Solubility
GeoParams.MaterialParameters.Solubility.GasMixture
```
[`SunYao2026_Solubility`](@ref) is the recommended law: it takes the melt composition, covers H2O and CO2 and is calibrated up to 2 GPa for H2O and 1 GPa for mixed H2O–CO2 fluids; above ~2 GPa its dissolved contents exceed 100 wt%. [`Liu2005_Solubility`](@ref) applies to rhyolite only and is useful to reproduce Degruyter & Huber (2014)-type models. [`Mafic_Solubility`](@ref) is a polynomial fit for mafic melts; its dissolved H2O does not vanish as the pressure goes to zero (2.7 wt% at 850 °C in the figure), so it is meaningful only at chamber pressures.

The figure compares the dissolved H2O in equilibrium with pure H2O gas at 850 °C.
```@example solubility_laws
using GeoParams, CairoMakie
basalt = (50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.14, 0.0)
laws = (
    "Liu2005_Solubility" => Liu2005_Solubility(),
    "Mafic_Solubility" => Mafic_Solubility(),
    "SunYao2026_Solubility, rhyolite" => SunYao2026_Solubility(),
    "SunYao2026_Solubility, basalt" => SunYao2026_Solubility(; oxd_wt = basalt),
)
P = range(1.0e5, 4.0e8, 200)
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; xlabel = "pressure P (MPa)", ylabel = "dissolved H₂O (wt%)")
for (name, law) in laws
    lines!(ax, P ./ 1.0e6, [100 * compute_dissolved(law, p, 1123.15, 0.0)[1] for p in P]; label = name)
end
axislegend(ax; position = :lt)
fig
```
## Computational routines
To evaluate dissolved H2O and CO2 (mass fractions) within a user routine, use this:
```@docs
GeoParams.MaterialParameters.Solubility.compute_dissolved
GeoParams.MaterialParameters.Solubility.compute_dissolved!
```
Partial derivatives (via `ForwardDiff`, useful for implicit/Newton solves) are available for pressure, temperature, and gas composition:
```@docs
GeoParams.MaterialParameters.Solubility.∂dissolved_∂P
GeoParams.MaterialParameters.Solubility.∂dissolved_∂T
GeoParams.MaterialParameters.Solubility.∂dissolved_∂Xco2
```
Given a target dissolved H2O content, the gas composition `X_co2` that produces it can be found with:
```@docs
GeoParams.MaterialParameters.Solubility.find_Xco2
```
Given both dissolved contents, the saturation pressure and the gas composition follow from:
```@docs
GeoParams.MaterialParameters.Solubility.compute_saturation_pressure
```
The figure applies it as a melt-inclusion barometer in rhyolite at 850 °C. The isobars are the dissolved contents in equilibrium with gases of all compositions at fixed pressure. Each of six illustrative inclusions, not measured data, is labeled with its saturation pressure and the depth ``z = P/(\rho g)`` for a crust of density ``\rho = 2600`` kg/m³. The right panel compares these pressures with those of [`Liu2005_Solubility`](@ref), which are up to 10 % higher.
```@example saturation_pressure
using GeoParams, CairoMakie
rhyolite, liu = SunYao2026_Solubility(), Liu2005_Solubility()
T = 1123.15                                  # 850 °C
ρ, g = 2600.0, 9.81                          # crust density (kg/m³), gravity (m/s²)
# illustrative melt inclusions, not data: (H2O wt%, CO2 ppm)
inclusions = [(2.5, 100.0), (3.0, 400.0), (4.0, 150.0), (4.5, 900.0), (5.5, 300.0), (6.0, 1200.0)]
fig = Figure(size = (1000, 520))
ax = Axis(fig[1, 1]; xlabel = "dissolved H₂O (wt%)", ylabel = "dissolved CO₂ (ppm)", title = "SunYao2026_Solubility, rhyolite, 850 °C")
X = range(0, 1, 201)
for (k, P) in enumerate((50.0e6, 100.0e6, 200.0e6, 300.0e6, 400.0e6))
    m = [compute_dissolved(rhyolite, P, T, x) for x in X]
    lines!(ax, 100 .* first.(m), 1.0e6 .* last.(m); color = k, colormap = :viridis, colorrange = (1, 5), label = "$(round(Int, P / 1.0e6)) MPa")
end
P_sy = [compute_saturation_pressure(rhyolite, h / 100, c / 1.0e6, T)[1] for (h, c) in inclusions]
P_liu = [compute_saturation_pressure(liu, h / 100, c / 1.0e6, T)[1] for (h, c) in inclusions]
scatter!(ax, inclusions; color = :black, label = "illustrative inclusions")
for ((h, c), P) in zip(inclusions, P_sy)
    text!(ax, h, c; text = "$(round(Int, P / 1.0e6)) MPa, $(round(P / (ρ * g) / 1.0e3; digits = 1)) km", offset = (6, 2), fontsize = 11)
end
Legend(fig[2, 1:2], ax; orientation = :horizontal, framevisible = false)
ax2 = Axis(fig[1, 2]; xlabel = "P, SunYao2026_Solubility (MPa)", ylabel = "P, Liu2005_Solubility (MPa)", title = "same inclusions, two laws", aspect = 1)
ablines!(ax2, 0, 1; color = :gray, linestyle = :dash)
scatter!(ax2, P_sy ./ 1.0e6, P_liu ./ 1.0e6; color = :black)
colsize!(fig.layout, 2, Relative(0.35))
fig
```
The equilibrium volume fraction of exsolved H2O gas, given the total water content, follows in closed form from a solubility law and the melt and gas densities of a [`ThreePhase_Density`](@ref):
```@docs
GeoParams.MaterialParameters.Solubility.compute_gas_fraction
```
The figure shows the gas fraction of crystal-free rhyolite (`ϕ = 1`) at 850 °C for the
[`SunYao2026_Solubility`](@ref) law and a [`HalbachChatterjee_Density`](@ref) gas. Gas exsolves below the
saturation pressure of each water content. The dashed line is the default critical gas fraction of [`CriticalPorosityFragmentation`](@ref).

```@example gas_fraction
using GeoParams, CairoMakie
ρ = ThreePhase_Density(; ρgas = HalbachChatterjee_Density())
solubility = SunYao2026_Solubility()
P = exp10.(range(5, log10(4.0e8), 200))
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; xscale = log10, xlabel = "pressure P (MPa)", ylabel = "gas volume fraction ϕ_gas (–)")
for m_h2o_total in (0.03, 0.05, 0.07)
    ϕ_gas = [compute_gas_fraction(solubility, ρ, (; P = p, T = 1123.15, ϕ = 1.0, m_h2o_total)) for p in P]
    lines!(ax, P ./ 1.0e6, ϕ_gas; label = "total H₂O = $(round(Int, 100 * m_h2o_total)) wt%")
end
hlines!(ax, 0.75; color = :black, linestyle = :dash, label = "ϕ_crit = 0.75")
axislegend(ax; position = :rc)
fig
```
Gas-mixture helpers for the energy equation:
```@docs
GeoParams.MaterialParameters.Solubility.compute_gas_heatcapacity
GeoParams.MaterialParameters.Solubility.effective_molar_mass
```

Only for pure H2O gas in equilibrium does the exsolved-gas volume fraction `ε_g` close algebraically at a point, through `compute_gas_fraction`. With CO2, `ε_g` and `X_co2` (CO2 mole fraction of the gas) are typically ODE state variables owned by the solver; the closures above are the *entries* to that mass balance, not a replacement for it.
