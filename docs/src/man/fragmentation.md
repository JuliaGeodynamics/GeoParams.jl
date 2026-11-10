# Magma fragmentation

Fragmentation criteria decide pointwise whether a bubbly magma has broken up into a
gas–pyroclast dispersion. Each returns a dimensionless ratio that is at least 1 where the
magma is fragmented:

```julia
julia> compute_fragmentation(CriticalPorosityFragmentation(), (; ϕ_gas = 0.8))
1.0666666666666667

julia> compute_fragmentation(StrainRateFragmentation(), (; εII = 1.0, η = 1.0e8))
1.1547005383792515
```

## Strain-rate criterion as a yield surface

For a magma with Newtonian viscosity ``\eta``, ``\tau_{II} = 2\eta\dot\varepsilon_{II}``, and the criterion of
[`StrainRateFragmentation`](@ref) becomes a von Mises yield stress that does not depend on ``\eta``,
```math
    \tau_{II} \ge \tau_y = \sqrt{3}\,k\,G_\infty .
```
The critical strain rate ``\tau_y/(2\eta)``, where the fragmentation ratio equals 1, is plotted below for the range ``G_\infty`` = 3–30 GPa of
Papale (1999) and ``k = 0.01``.

```@example fragmentation
using GeoParams, CairoMakie
η = 10 .^ range(4, 12, 100)
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; xlabel = "viscosity η (Pa s)", ylabel = "critical strain rate ε̇II (ratio = 1) (s⁻¹)",
    xscale = log10, yscale = log10, xticks = LogTicks(4:2:12))
for G in (3.0e9, 1.0e10, 3.0e10)
    f = StrainRateFragmentation(G = G * Pa)
    τ_y = sqrt(3) * f.k * G
    lines!(ax, η, τ_y ./ (2 .* η); label = "G∞ = $(G / 1e9) GPa, yield stress τII = $(round(τ_y / 1e6; digits = 1)) MPa")
end
axislegend(ax; position = :rt)
fig
```

`DruckerPrager_regularised(f; η_vp)` returns this yield surface. In a viscoelastic
`CompositeRheology` whose elastic modulus is ``G_\infty``, the plastic return mapping then caps
``\tau_{II}`` at ``\tau_y`` and returns a positive plastic multiplier ``\dot\lambda`` where the magma fragments.
The regularization allows stresses above ``\tau_y`` by about ``\eta_{vp}\dot\lambda``, so `η_vp` must be much
smaller than the magma viscosity; a value such as `η_vp = 1e20Pa*s` leaves the stress of a magma practically uncapped.

```@example fragmentation
frag = StrainRateFragmentation(G = 1.0e10Pa)
magma = CompositeRheology(
    LinearViscous(η = 1.0e8Pa * s), ConstantElasticity(G = 1.0e10Pa), DruckerPrager_regularised(frag; η_vp = 1.0e5Pa * s)
)
τ_y = sqrt(3) * frag.k * 1.0e10
args = (; dt = 1.0e3, τII_old = 0.0, P = 0.0)
for εII in (0.01, 1.0, 10.0)
    τII, λ̇ = compute_τII(magma, εII, args)
    println("εII = $εII s⁻¹: τII/τ_y = $(round(τII / τ_y; digits = 3)), λ̇ = $(round(λ̇; sigdigits = 3)) s⁻¹")
end
```

```@docs
GeoParams.Fragmentation.AbstractFragmentation
GeoParams.Fragmentation.CriticalPorosityFragmentation
GeoParams.Fragmentation.StrainRateFragmentation
GeoParams.Fragmentation.compute_fragmentation
GeoParams.Fragmentation.DruckerPrager_regularised(::GeoParams.Fragmentation.StrainRateFragmentation)
```
