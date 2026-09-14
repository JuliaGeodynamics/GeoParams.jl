```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: GeoParams.jl Docs
  text: Material parameters for geoscientific models.
  tagline: Easily extendable material parameters for the direct use in numerical simulations
  actions:
    - theme: brand
      text: Nondimensionalization
      link: /man/nondimensionalize
    - theme: alt
      text: Material Parameters 📚
      link: /man/materialparameters
    - theme: alt
      text: Constitutive Relationships 🎯
      link: /man/creeplaws
    - theme: alt
      text: API Reference 📚
      link: /man/listfunctions
  image:
    src: /logo.png
    alt: GeoParams.jl

features:
  - icon: 🚀
    title: Creep laws
    details: Effortlessly switch between linear and nonlinear creep laws.
    link: /man/creeplaws

  - icon: ⚡
    title: Chemical diffusion
    details: Calculate chemical diffusion coefficients in different minerals and melts.
    link: man/chemicaldiffusion

  - icon: 📈
    title: Plotting
    details: Various plotting routines.
    link: man/plotting

  - icon: 🧩
    title: Extensibility
    details: Provides a natural repository for contributions of various new creep laws and other rheological features for use by the larger community.
    link: /man/contributing
---
```

## GeoParams.jl

Typical geodynamic simulations involve a large number of material parameters that have units that are often inconvenient to be directly used in numerical models. This package has three main features that help with this:

- Nondimensionalization object, which can be used to transfer dimensional to non-dimensional parameters (usually better for numerical solvers).
- Material parameters object in which you can specify parameters employed in the geodynamic simulations. This object is designed to be extensible and can be passed on to the solvers, such that new creep laws or features can be readily added.
- Allocation-free computational routines for CPUs and GPUs, which replace the point-wise calculations in your solver.

We also implement some typically used creep law parameters, together with tools to plot them versus and compare our results with those of published papers (to minimize mistakes).

### Installation

GeoParams is registered, so it installs with the package manager:

```julia
julia> using Pkg; Pkg.add("GeoParams")
```

Plotting routines live in a package extension, so they become available once a
[Makie](https://docs.makie.org/stable/) backend is loaded alongside GeoParams.

### Getting started

Pick characteristic values, then use them to nondimensionalize a quantity:

```jldoctest
julia> using GeoParams

julia> CharDim = GEO_units(length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
Employing GEO units 
Characteristic values: 
         length:      1000 km
         time:        0.3169 Myr
         stress:      10 MPa
         temperature: 1000.0 °C

julia> nondimensionalize(10MPa, CharDim)
0.9999999999999998
```

Material parameters for one phase are collected in a single object that the
computational routines accept directly:

```jldoctest
julia> using GeoParams

julia> phase = SetMaterialParams(Name="Matrix", Phase=1,
                                 Density   = ConstantDensity(ρ=2900kg/m^3),
                                 CreepLaws = LinearViscous(η=1e23Pa*s));

julia> compute_density(phase.Density[1], (;))
2900.0
```

The [Nondimensionalization](man/nondimensionalize.md) and
[Material Parameters](man/materialparameters.md) pages go into more detail.

### Funding
The development of this package was supported by the European Research Council (ERC CoG #771143 MAGMA) as well as by the [GPU4GEO](https://ptsolvers.github.io/GPU4GEO/) [PASC](https://www.pasc-ch.org) project.
