# Magma suspensions

[`Suspension`](@ref) wraps the creep law of a melt and multiplies its viscosity by the relative viscosities
``\eta_r = \eta/\eta_\mathrm{melt}`` of suspended crystals and bubbles. The crystal and gas volume fractions of
the whole suspension are passed at runtime as `ϕ_x` and `ϕ_gas`:

```@example suspensions_core
using GeoParams, CairoMakie
law = Suspension(LinearViscous(η = 1.0e3Pa * s); crystals = Roscoe1952_Crystals(), bubbles = IshiiZuber1979_Bubbles())
compute_relative_viscosity(law; ϕ_x = 0.3, ϕ_gas = 0.1)
```

## Rigid crystals and undeformable bubbles

[`Roscoe1952_Crystals`](@ref) describes rigid particles that lock up at the maximum packing fraction
``\phi_\mathrm{max}``, and [`IshiiZuber1979_Bubbles`](@ref) describes bubbles that stay spherical. The figure compares
Roscoe's (1952) parameters for equal spheres, the Krieger–Dougherty form fitted to crystal-bearing basalt by
Keller & Suckale (2019), and Ishii–Zuber bubbles as a function of their volume fraction.

```@example suspensions_core
ϕ = range(0, 0.74, 300)
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; xlabel = "crystal or gas volume fraction ϕ (–)", ylabel = "relative viscosity ηᵣ (–)", yscale = log10)
lines!(ax, ϕ, compute_relative_viscosity.(Ref(Roscoe1952_Crystals()), ϕ); label = "Roscoe1952_Crystals(), Roscoe (1952)")
lines!(ax, ϕ, compute_relative_viscosity.(Ref(Roscoe1952_Crystals(ϕ_max = 0.58, n = 2.32)), ϕ);
    label = "Roscoe1952_Crystals(ϕ_max = 0.58, n = 2.32), Keller & Suckale (2019)")
lines!(ax, ϕ, compute_relative_viscosity.(Ref(IshiiZuber1979_Bubbles()), ϕ); label = "IshiiZuber1979_Bubbles() (bubbles)")
ylims!(ax, 1, 1.0e4)
axislegend(ax; position = :lt)
fig
```

```@docs
GeoParams.MaterialParameters.ConstitutiveRelationships.Suspension
GeoParams.MaterialParameters.ConstitutiveRelationships.Roscoe1952_Crystals
GeoParams.MaterialParameters.ConstitutiveRelationships.IshiiZuber1979_Bubbles
GeoParams.compute_relative_viscosity
```

## Particle shape, shear rate and bubble deformation

The laws of this section account for the shape of crystals and for the dependence of ``\eta_r`` on the
strain rate of the suspension.

### Elongated crystals

[`Mueller2010_Crystals`](@ref) lowers the maximum packing fraction with the particle aspect ratio ``r_p``, so elongated
crystals stiffen the melt at much lower crystal fractions than spheres. The curve for ``r_p = 10`` extrapolates
slightly beyond the fitted range ``1 \le r_p \le 9.17``.

```@example suspensions_rate
using GeoParams, CairoMakie
ϕ_x = range(0, 0.74, 300)
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; xlabel = "crystal volume fraction ϕₓ (–)", ylabel = "relative viscosity ηᵣ (–)", yscale = log10)
for r_p in (1, 5, 10)
    lines!(ax, ϕ_x, compute_relative_viscosity.(Ref(Mueller2010_Crystals(aspect_ratio = r_p)), ϕ_x); label = "Mueller2010_Crystals, rₚ = $r_p")
end
lines!(ax, ϕ_x, compute_relative_viscosity.(Ref(Roscoe1952_Crystals()), ϕ_x); color = :black, linestyle = :dash,
    label = "Roscoe1952_Crystals(), Roscoe (1952)")
ylims!(ax, 1, 1.0e4)
axislegend(ax; position = :lt)
fig
```

### Shear-thinning crystal suspensions

[`Vona2011_Crystals`](@ref) makes the Maron–Pierce law of [`Mueller2010_Crystals`](@ref) shear thinning: ``\eta_r`` decreases with
the shear rate ``\dot\gamma = 2\dot\varepsilon_{II}`` and equals the Mueller et al. (2010) value at
``\dot\gamma = 1`` s⁻¹. The figure shows ``\eta_r`` relative to that value. The law was calibrated for
``0.07 \le \phi_x \le 0.27`` (shaded) and shear rates of 0.02–2.2 s⁻¹; elsewhere it is an extrapolation.

```@example suspensions_rate
ϕ_x = range(0, 0.35, 300)
f = Vona2011_Crystals()
η_r(ϕ, γ̇) = compute_relative_viscosity(f, ϕ, γ̇)
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; xlabel = "crystal volume fraction ϕₓ (–)", ylabel = "ηᵣ(γ̇) / ηᵣ(1 s⁻¹) (–)")
vspan!(ax, 0.07, 0.27; color = (:gray, 0.2), label = "calibrated range of ϕₓ")
hlines!(ax, 1; color = :black, linestyle = :dash, label = "γ̇ = 1 s⁻¹ (Mueller2010_Crystals, rₚ = 6.5)")
for γ̇ in (0.02, 0.2, 2.0)
    lines!(ax, ϕ_x, η_r.(ϕ_x, γ̇) ./ η_r.(ϕ_x, 1.0); label = "Vona2011_Crystals(), γ̇ = $γ̇ s⁻¹")
end
xlims!(ax, 0, 0.35)
axislegend(ax; position = :lt)
fig
```

### Deformable bubbles

Bubbles stiffen the melt while surface tension keeps them spherical and weaken it once the flow stretches them.
[`Mader2013_Bubbles`](@ref) interpolates between the two limits as a function of the capillary number ``\mathrm{Ca}``
(see [`AbstractCapillaryViscosity`](@ref GeoParams.MaterialParameters.ConstitutiveRelationships.AbstractCapillaryViscosity)):
``\eta_{r0} = 1/(1-\phi)``, which equals [`IshiiZuber1979_Bubbles`](@ref), for ``\mathrm{Ca} \to 0`` (dashed), and
``\eta_{r\infty} = (1-\phi)^{5/3}`` for ``\mathrm{Ca} \to \infty`` (dotted).

```@example suspensions_rate
f = Mader2013_Bubbles()
Ca = 10 .^ range(-2, 2, 300)
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; xlabel = "capillary number Ca (–)", ylabel = "relative viscosity ηᵣ (–)", xscale = log10)
ϕ_gas = (0.1, 0.3, 0.5)
colors = Makie.wong_colors()[1:3]
for (ϕ, color) in zip(ϕ_gas, colors)
    lines!(ax, Ca, compute_relative_viscosity.(Ref(f), ϕ, Ca); color)
    hlines!(ax, compute_relative_viscosity(IshiiZuber1979_Bubbles(), ϕ); color, linestyle = :dash)
    hlines!(ax, (1 - ϕ)^(5 / 3); color, linestyle = :dot)
end
elements = [[LineElement(; color) for color in colors]; LineElement(color = :black, linestyle = :dash);
    LineElement(color = :black, linestyle = :dot)]
labels = [["Mader2013_Bubbles(), ϕ_gas = $ϕ" for ϕ in ϕ_gas]; "Ca → 0: (1 − ϕ)⁻¹ (IshiiZuber1979_Bubbles)"; "Ca → ∞: (1 − ϕ)^(5/3)"]
axislegend(ax, elements, labels; position = :rt)
fig
```

!!! warning "Stress decreasing with strain rate"
    For a Newtonian melt, the capillary number is proportional to the strain rate, so the stress of the
    suspension is proportional to ``\eta_r(\mathrm{Ca})\,\mathrm{Ca}``. With ``y = (K\,\mathrm{Ca})^2`` its slope,
    ``\eta_{r\infty} + (\eta_{r0} - \eta_{r\infty})(1 - y)/(1 + y)^2``, is smallest at ``y = 3``, where it equals
    ``\eta_{r\infty} - (\eta_{r0} - \eta_{r\infty})/8``. It is negative if ``9\,\eta_{r\infty} < \eta_{r0}``, i.e.
    ``(1 - \phi)^{8/3} < 1/9`` or ``\phi > 1 - 9^{-3/8} \approx 0.561``, where ``\phi`` is the gas fraction
    ``\phi_\mathrm{gas}/(1 - \phi_x)`` of the bubbly melt. The stress then decreases with strain rate over a range
    of strain rates, which makes the momentum equation ill-posed and the strain rate at a given stress non-unique.

The threshold can be checked numerically:

```@example suspensions_rate
stress_decreases(ϕ_gas) = any(diff(compute_relative_viscosity.(Ref(Mader2013_Bubbles()), ϕ_gas, Ca) .* Ca) .< 0)
stress_decreases(0.55), stress_decreases(0.57)
```

### Three-phase suspensions

Following Truby et al. (2015), [`Suspension`](@ref) suspends the crystals in the bubbly melt, whose gas fraction
is ``\phi_b = \phi_\mathrm{gas}/(1 - \phi_x)``. At constant ``\phi_\mathrm{gas}`` the bubbles therefore stiffen
the magma more as the crystal fraction grows. Truby et al. validated this combination for spheres
(``r_p = 1``) and undeformed bubbles up to ``\phi_x \approx 0.5`` and ``\phi_b \approx 0.3``, the range plotted here.

```@example suspensions_rate
law = Suspension(LinearViscous(); crystals = Mueller2010_Crystals(aspect_ratio = 1), bubbles = IshiiZuber1979_Bubbles())
ϕ_x = range(0, 0.5, 300)
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; xlabel = "crystal volume fraction ϕₓ (–)", ylabel = "relative viscosity ηᵣ (–)", yscale = log10)
for ϕ_gas in (0.0, 0.05, 0.1, 0.15)
    lines!(ax, ϕ_x, [compute_relative_viscosity(law; ϕ_x = x, ϕ_gas) for x in ϕ_x]; label = "ϕ_gas = $ϕ_gas")
end
axislegend(ax; position = :lt)
fig
```

```@docs
GeoParams.MaterialParameters.ConstitutiveRelationships.Mueller2010_Crystals
GeoParams.MaterialParameters.ConstitutiveRelationships.Vona2011_Crystals
GeoParams.MaterialParameters.ConstitutiveRelationships.AbstractCapillaryViscosity
GeoParams.MaterialParameters.ConstitutiveRelationships.Mader2013_Bubbles
```
