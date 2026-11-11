# Density

## Methods
The density equation of state can be specified in a number of ways
```@docs
GeoParams.MaterialParameters.Density.ConstantDensity
GeoParams.MaterialParameters.Density.PT_Density
GeoParams.MaterialParameters.Density.T_Density
GeoParams.MaterialParameters.Density.Compressible_Density
GeoParams.MaterialParameters.Density.MeltDependent_Density
GeoParams.MaterialParameters.Density.BubbleFlow_Density
GeoParams.MaterialParameters.Density.GasPyroclast_Density
GeoParams.MaterialParameters.Density.Melt_DensityX
GeoParams.MaterialParameters.Density.RedlichKwong_Density
GeoParams.MaterialParameters.Density.IdealGas_Density
GeoParams.MaterialParameters.Density.HalbachChatterjee_Density
GeoParams.MaterialParameters.Density.ThreePhase_Density
```
## H2O gas density
The figure compares the H2O gas laws at 850 °C. [`IdealGas_Density`](@ref) is close to
[`HalbachChatterjee_Density`](@ref) below about 10 MPa, but is 10 % too low near 100 MPa and 27 % too high at
400 MPa. [`RedlichKwong_Density`](@ref) is a fit stated for 30–400 MPa; it is off by up to 35 % near 30 MPa and
diverges below it.
```@example gas_density
using GeoParams, CairoMakie
T = 1123.15
P = exp10.(range(5, log10(4.0e8), 300))
hc, rk, ig = HalbachChatterjee_Density(), RedlichKwong_Density(), IdealGas_Density()
ρ(law) = [compute_density(law, (; P = p, T)) for p in P]
ρ_hc, ρ_rk, ρ_ig = ρ(hc), ρ(rk), ρ(ig)
fit = P .≥ 3.0e7    # RedlichKwong_Density is stated for 30–400 MPa
fig = Figure(size = (1000, 420))
ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "pressure P (MPa)", ylabel = "ρ (kg/m³)", title = "H2O gas density at 850 °C")
ax2 = Axis(fig[1, 2]; xscale = log10, xlabel = "pressure P (MPa)", ylabel = "ρ / ρ_HalbachChatterjee (–)", title = "relative to HalbachChatterjee_Density")
for (ax, f) in ((ax1, identity), (ax2, r -> r ./ ρ_hc))
    lines!(ax, P ./ 1.0e6, f(ρ_hc); label = "HalbachChatterjee_Density")
    lines!(ax, P ./ 1.0e6, f(ρ_rk); color = Cycled(2), linestyle = :dash, label = "RedlichKwong_Density, outside its fit")
    lines!(ax, P[fit] ./ 1.0e6, f(ρ_rk)[fit]; color = Cycled(2), label = "RedlichKwong_Density")
    lines!(ax, P ./ 1.0e6, f(ρ_ig); color = Cycled(3), label = "IdealGas_Density")
end
ylims!(ax1, 0.1, 2.0e3)
ylims!(ax2, 0.4, 1.6)
axislegend(ax1; position = :lt)
fig
```

## Computational routines
To evaluate density within a user routine, use this:
```@docs
GeoParams.MaterialParameters.Density.compute_density
GeoParams.MaterialParameters.Density.compute_density!
```
Note that density values are usually not used in itself in the governing PDE's, but usually in combination with other parameters, such as $\rho g$ or $\rho c_p$. The non-dimensional value of $\rho$ may thus have very large or small values, but multiplied with the other values one often obtains numbers that are closer to one.
