# GeoParams.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliageodynamics.github.io/GeoParams.jl/dev/)
[![Build Status](https://github.com/JuliaGeodynamics/GeoParams.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeodynamics/GeoParams.jl/actions)

Typical geodynamic simulations involve a large number of material parameters that have units that are often inconvenient to be directly used in numerical models
This package has two main features that help with this:
- Create a nondimensionalization object, which can be used to transfer dimensional to non-dimensional parameters (usually better for numerical solvers)
- Create an object in which you can specify material parameters employed in the geodynamic simulations

The material parameter object is designed to be extensible and can be passed on to the solvers, such that new creep laws or features can be readily added. 
We also implement some typically used creep law parameters, together with tools to plot them versus and compare our results with those of published papers (to minimize mistakes).  

### Contents
* [1. Nondimensionalization](#1-nondimensionalization) 
* [2. Material parameters](#2-material-parameters)
* [3. Plotting and output](#3-plotting-and-output)
* [4. Installation](#4-installation)
* [5. Dependencies](#5-dependencies)
* [6. Contributing](#6-contributing)
* [7. Funding](#7-funding)

### 1. Nondimensionalization 
Typical geodynamic simulations involve dimensions on the order of 10's-1000's of kilometers, and viscosities on the order of ~1e20 Pas. If such values are directly employed in numerical solvers, they may result in roundoff errors. It is therefore common practice to nondimensionalize the input parameters by dividing them by typical values such that the result gives numbers that are closer to one.
This can be done by specifying characteristic values for `length`, `stress`, `temperature` and `viscosity`. From these `basic` units all other physical units are derived and input parameters can thus be nondimensionalized accordingly (and dimensionalized again when plotting or saving output). 

As you learned in physics, the common approach to do this is by using `SI` units. Yet, as meters and seconds are not so convenient in geodynamics, where we usually deal with lengthscales on the orders of kilometers and timescales in millions of years, we also provide the `geo` object, which allows to give input parameters in more convenient units.

```julia
julia> using GeoParams
julia> CharDim = GEO_units(length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
Employing GEO units 
Characteristic values: 
         length:      1000 km
         time:        0.3169 Myrs
         stress:      10 MPa
         temperature: 1000.0 °C
```
You can use 3 `types`:
  1. *GEO* units: Units of length in the code are expected to be in kilometers, time is in million of years (Myrs) and stresses are in MPa (1e6 Pa). 
  2. *SI* units: all values are in SI units (meters, Pascal, seconds)
  3. *NONE*: all input parameters are in nondimensional units

Once a `CharDim` structure is created, you can use the derived parameters, for example:
```julia
julia> CharDim.strainrate
1.0e-13 s⁻¹
```
You can also non-dimensionalize parameters:
```julia
julia> A    =   6.3e-2MPa^-3.05*s^-1
0.063 MPa⁻³·⁰⁵ s⁻¹
julia> A_ND =   Nondimensionalize(A, CharDim);
```
or convert them to different units:
```julia
julia> uconvert(Pa^-3.05*s^-1, A)
3.157479571851836e-20 Pa⁻³·⁰⁵ s⁻¹
```
### 2. Material parameters  
  
All geodynamic simulations require specifying material parameters, such as (nonlinear) viscous constitutive relationships or an equation of state. These parameters are usually specified per `phase`. Here, we provide a framework that simplifies doing that. Thanks to the flexibility of julia, we can actually directly embed the function that does the computations in the structure itself, which makes it straightforward to extend it and add new creep laws (which can directly be used in the solvers).  

Some examples of where this is used:
#### 2.1 Constant density, constant linear viscosity
```julia
julia> Phase = SetMaterialParams(Name="Viscous Matrix", Phase=2,
                                     Density   = ConstantDensity(),
                                     CreepLaws = LinearViscous(η=1e23Pa*s))
Phase 2 : Viscous Matrix
        | [dimensional units]
        | 
        |-- Density           : Constant density: ρ=2900 kg m⁻³ 
        |-- Gravity           : Gravitational acceleration: g=9.81 m s⁻² 
        |-- CreepLaws         : Linear viscosity: η=1.0e23 Pa s 
```
The same but with non-dimensionalization of all parameters:
```julia
julia> CharDim = GEO_units(length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas);
julia> Phase = SetMaterialParams(Name="Viscous Matrix", Phase=2, 
                                     Density   = ConstantDensity(),
                                     CreepLaws = LinearViscous(η=1e23Pa*s), CharDim=CharDim)
Phase 2 : Viscous Matrix
        | [non-dimensional units]
        | 
        |-- Density           : Constant density: ρ=2.8999999999999996e-18 
        |-- Gravity           : Gravitational acceleration: g=9.81e20 
        |-- CreepLaws         : Linear viscosity: η=999.9999999999998 
```

#### 2.2 Nonlinear creep laws
to be added


### 3. Plotting and output
- WORK IN PROGRESS

A typical geodynamic simulation involves a lot of parameters. Creating data tables for scientific publications that describe all parameters employed is usually done by hand (and no-one really likes doing that). In our experience a lot of errors happen while doing this, either because the units are mixed up (some creep laws have weird units like MPa^{-n}), or because some parameters are forgotten. To help with this, we provide number of functions that 
  1)  Simplify creating plots in the same manner as in many publications that report the laboratory experiments used to create the creep laws. In that way, they can be directly compared to the original results. You can also create publication-ready figures.
  2)  Provide tools to automatically generate data tables from the input parameters. This saves time and minimizes errors. 
#### 3.1 Plotting creep laws 
A few simple functions are provided to plot creep laws.

#### 3.2 Creating data tables
to be developed


### 4. Installation
You can install this package by specifying 
```julia
julia> ]
pkg> add https://github.com/JuliaGeodynamics/GeoParams.jl
```
and test whether it works with
```
pkg> test GeoParams
```

### 5. Dependencies
We rely on:
- [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) to deal with SI units
- [Parameters.jl](https://github.com/mauro3/Parameters.jl) to have structures that are easier to modify
- [LaTeXStrings.jl](https://github.com/stevengj/LaTeXStrings.jl) to be able to add equations to the structures that describe the employed material laws
### 6. Contributing
Help with developing this package is highly appreciated. You can contribute for example by adding new creep laws or by adding new constitutive relationships. If you invest a bit of time now, it will save others in the community a lot of time! 
The simplest way to do this is by cloning the repository, and creating a new branch for your feature. Once you are happy with what you added (and after you added a test to ensure that it will keep working with future changes), create a pull request and we will evaluate & merge it.


### 7. Funding
The development of this package was supported by the European Research Council (ERC CoG #771143).