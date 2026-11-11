# Heat capacity

## Methods
Heat capacity is defined as 
```@docs
GeoParams.MaterialParameters.HeatCapacity.ConstantHeatCapacity
GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Whittington
GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Moitra
GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Moitra_smoothed
GeoParams.MaterialParameters.HeatCapacity.ThreePhase_HeatCapacity
```
## Heat capacity of magma
The glass and melt heat capacity of [`T_HeatCapacity_Moitra`](@ref GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Moitra) jumps at the glass transition ``T_g``;
[`T_HeatCapacity_Moitra_smoothed`](@ref GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Moitra_smoothed) replaces the jump by a transition of width `ΔTg`. The fit is calibrated up to 1273 K.
[`T_HeatCapacity_Whittington`](@ref GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Whittington) is shown for comparison.
```@example cp_moitra
using GeoParams, CairoMakie
T = range(300, 1400, 1101)
laws = (
    "T_HeatCapacity_Moitra" => T_HeatCapacity_Moitra(),
    "T_HeatCapacity_Moitra_smoothed (ΔTg = 20 K)" => T_HeatCapacity_Moitra_smoothed(),
    "T_HeatCapacity_Whittington" => T_HeatCapacity_Whittington(),
)
fig = Figure(size = (800, 450))
ax = Axis(fig[1, 1]; xlabel = "temperature T (K)", ylabel = "heat capacity c_p (J/kg/K)")
vspan!(ax, 1273, 1400; color = (:gray, 0.2), label = "Moitra extrapolated (> 1273 K)")
vlines!(ax, 1000; color = :black, linestyle = :dot, label = "T_g = 1000 K")
for (name, law) in laws
    lines!(ax, T, [compute_heatcapacity(law, (; T = t)) for t in T]; label = name)
end
axislegend(ax; position = :lt)
fig
```

[`ThreePhase_HeatCapacity`](@ref GeoParams.MaterialParameters.HeatCapacity.ThreePhase_HeatCapacity) mass-weights the phase heat capacities. The figure shows magma with 0–40 % crystals at
100 MPa and 850 °C, with an H2O gas of [`HalbachChatterjee_Density`](@ref). The gas has a high specific heat capacity, so
``c_p`` rises with the gas fraction, but it has little mass, so the volumetric heat capacity ``\rho c_p`` that enters
the energy equation falls almost linearly.
```@example cp_threephase
using GeoParams, CairoMakie
ρ = ThreePhase_Density(; ρgas = HalbachChatterjee_Density())
cp = ThreePhase_HeatCapacity(; ρ)
ϕ_gas = range(0, 0.6, 121)
fig = Figure(size = (900, 400))
ax1 = Axis(fig[1, 1]; xlabel = "gas fraction ϕ_gas (–)", ylabel = "c_p (J/kg/K)", title = "mass-weighted c_p")
ax2 = Axis(fig[1, 2]; xlabel = "gas fraction ϕ_gas (–)", ylabel = "ρ c_p (MJ/m³/K)", title = "volumetric heat capacity ρ c_p")
for ϕ_x in (0.0, 0.2, 0.4)
    args = [(; P = 1.0e8, T = 1123.15, ϕ_x, ϕ_gas = g) for g in ϕ_gas]
    lines!(ax1, ϕ_gas, [compute_heatcapacity(cp, a) for a in args]; label = "ϕ_x = $ϕ_x")
    lines!(ax2, ϕ_gas, [compute_density(ρ, a) * compute_heatcapacity(cp, a) / 1.0e6 for a in args]; label = "ϕ_x = $ϕ_x")
end
axislegend(ax2; position = :lb)
fig
```

## Computational routines
To compute, use this:
```@docs
GeoParams.MaterialParameters.HeatCapacity.compute_heatcapacity
GeoParams.MaterialParameters.HeatCapacity.compute_heatcapacity!
```
