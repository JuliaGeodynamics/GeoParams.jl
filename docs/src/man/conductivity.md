# Conductivity

## Methods
Thermal conductivity is defined as 
```@docs
GeoParams.MaterialParameters.Conductivity.ConstantConductivity
GeoParams.MaterialParameters.Conductivity.T_Conductivity_Whittington
GeoParams.MaterialParameters.Conductivity.T_Conductivity_Whittington_parameterised
GeoParams.MaterialParameters.Conductivity.TP_Conductivity
GeoParams.MaterialParameters.Conductivity.Set_TP_Conductivity
GeoParams.MaterialParameters.Conductivity.T_Conductivity_Moitra
GeoParams.MaterialParameters.Conductivity.T_Conductivity_Moitra_smoothed
GeoParams.MaterialParameters.Conductivity.IAPWS_H2O_Conductivity
GeoParams.MaterialParameters.Conductivity.GasMixture_Conductivity
GeoParams.MaterialParameters.Conductivity.Maxwell_Conductivity
```
## Conductivity of magma
[`T_Conductivity_Moitra`](@ref GeoParams.MaterialParameters.Conductivity.T_Conductivity_Moitra) computes ``k = \rho c_p \alpha`` from the measured diffusivity of glass and melt, so it jumps at
``T_g`` like the heat capacity; [`T_Conductivity_Moitra_smoothed`](@ref GeoParams.MaterialParameters.Conductivity.T_Conductivity_Moitra_smoothed) uses the same transition.
```@example k_moitra
using GeoParams, CairoMakie
T = range(300, 1400, 1101)
fig = Figure(size = (800, 450))
ax = Axis(fig[1, 1]; xlabel = "temperature T (K)", ylabel = "conductivity k (W/m/K)")
vspan!(ax, 1273, 1400; color = (:gray, 0.2), label = "extrapolated (> 1273 K)")
vlines!(ax, 1000; color = :black, linestyle = :dot, label = "T_g = 1000 K")
for (name, law) in ("T_Conductivity_Moitra" => T_Conductivity_Moitra(), "T_Conductivity_Moitra_smoothed (ΔTg = 20 K)" => T_Conductivity_Moitra_smoothed())
    lines!(ax, T, [compute_conductivity(law, (; T = t)) for t in T]; label = name)
end
axislegend(ax; position = :lt)
fig
```

The conductivity of the gas phase: [`IAPWS_H2O_Conductivity`](@ref GeoParams.MaterialParameters.Conductivity.IAPWS_H2O_Conductivity) is solid where IAPWS R15-11 validates the formulation and
dashed where it is extrapolated. [`GasMixture_Conductivity`](@ref GeoParams.MaterialParameters.Conductivity.GasMixture_Conductivity) mixes it with the default constant CO2 conductivity.
```@example k_gas
using GeoParams, CairoMakie
# IAPWS R15-11 validates the formulation for T ≤ 1173.15 K at P ≤ 100 MPa, T ≤ 874 K at P ≤ 250 MPa, T ≤ 573 K above
validated(P, T) = (P ≤ 1.0e8 && T ≤ 1173.15) || (P ≤ 2.5e8 && T ≤ 874.0) || T ≤ 573.0
kh2o = IAPWS_H2O_Conductivity()
P = exp10.(range(5, log10(4.0e8), 200))
fig = Figure(size = (1000, 420))
ax1 = Axis(fig[1, 1]; xscale = log10, xlabel = "pressure P (MPa)", ylabel = "k (W/m/K)", title = "IAPWS_H2O_Conductivity (dashed: extrapolated)")
for (i, Tc) in enumerate((600, 850, 1100))
    T = Tc + 273.15
    k = [compute_conductivity(kh2o, (; P = p, T)) for p in P]
    lines!(ax1, P ./ 1.0e6, k; color = Cycled(i), linestyle = :dash)
    lines!(ax1, P ./ 1.0e6, [validated(p, T) ? ki : NaN for (p, ki) in zip(P, k)]; color = Cycled(i), label = "$Tc °C")
end
axislegend(ax1; position = :lt)
ax2 = Axis(fig[1, 2]; xlabel = "CO2 mole fraction of the gas X_co2 (–)", ylabel = "k (W/m/K)", title = "GasMixture_Conductivity, 850 °C")
X = range(0, 1, 101)
for p in (1.0e7, 1.0e8, 3.0e8)
    lines!(ax2, X, [compute_conductivity(GasMixture_Conductivity(), (; P = p, T = 1123.15, X_co2 = x)) for x in X]; label = "$(round(Int, p / 1.0e6)) MPa")
end
axislegend(ax2; position = :rt)
fig
```

[`Maxwell_Conductivity`](@ref GeoParams.MaterialParameters.Conductivity.Maxwell_Conductivity) combines a matrix with gas bubbles. With its defaults, a Moitra matrix and an H2O gas, the
bubbles lower the magma conductivity roughly in proportion to the gas fraction.
```@example k_maxwell
using GeoParams, CairoMakie
k = Maxwell_Conductivity()
ϕ_gas = range(0, 0.8, 161)
fig = Figure(size = (800, 450))
ax = Axis(fig[1, 1]; xlabel = "gas fraction ϕ_gas (–)", ylabel = "conductivity k (W/m/K)", title = "Maxwell_Conductivity: Moitra matrix with H2O bubbles at 50 MPa")
for T in (900.0, 1200.0)
    lines!(ax, ϕ_gas, [compute_conductivity(k, (; P = 5.0e7, T, ϕ_gas = g)) for g in ϕ_gas]; label = "T = $(round(Int, T)) K")
end
axislegend(ax; position = :rt)
fig
```

## Computational routines
To compute, use this:
```@docs
GeoParams.MaterialParameters.Conductivity.compute_conductivity
GeoParams.MaterialParameters.Conductivity.compute_conductivity!
```
