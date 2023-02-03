<h1> <img src="docs/src/assets/logo.png" alt="GeoParams.jl" width="50"> GeoParams.jl </h1>

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliageodynamics.github.io/GeoParams.jl/dev/)
[![Build Status](https://github.com/JuliaGeodynamics/GeoParams.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeodynamics/GeoParams.jl/actions)

Typical geodynamic simulations involve a large number of material parameters and nonlinear constitutive relationships. A large part of the work in writing a new code is benchmarking and debugging the implementation of such material parameters, which involve *point-wise* calculations that are independent of the discretisation method (finite difference, finite element, finite volume).

This package has three main objectives:
- Create a nondimensionalization object, which can be used to transfer dimensional to non-dimensional parameters (helps numerical solvers)
- Create an object in which you can specify material parameters employed in the geodynamic simulations
- Provide allocation-free computational routines for GPU and CPUs, that can be integrated in solvers, which replaces all point-wise calculations (to compute material parameters or equations of state, for example).

The material parameter object is designed to be extensible and can be passed on to the solvers, such that new creep laws or features can be readily added. If you use the computational routines we provide, these new features are immediately available in all your codes (finite element, finite difference, AMR codes, etc.). 

We also implement some typically used creep law parameters, together with tools to plot them versus and compare our results with those of published papers (to minimize mistakes).  

NOTE: The package remains under development and the API is not yet fully fixed. Therefore feel free to look at it, but be aware that things may still change when you incorporate it into your codes. Comments/ideas/suggestions are highly apprecciated!

### Contents
* [1. Nondimensionalization](#1-nondimensionalization) 
* [2. Material parameters](#2-material-parameters)
* [3. Plotting and output](#3-plotting-and-output)
* [4. Computational engine](#4-computational-engine)
* [5. Installation](#5-installation)
* [6. Documentation](#6-documentation)
* [7. Dependencies](#7-dependencies)
* [8. Contributing](#8-contributing)
* [9. Funding](#9-funding)

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
julia> A_ND =   nondimensionalize(A, CharDim);
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
julia> MatParam = SetMaterialParams(Name="Viscous Matrix", Phase=2,
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
julia> MatParam = SetMaterialParams(Name="Viscous Matrix", Phase=2, 
                                     Density   = ConstantDensity(),
                                     CreepLaws = LinearViscous(η=1e23Pa*s), CharDim=CharDim)
Phase 2 : Viscous Matrix
        | [non-dimensional units]
        | 
        |-- Density           : Constant density: ρ=2.8999999999999996e-18 
        |-- Gravity           : Gravitational acceleration: g=9.81e20 
        |-- CreepLaws         : Linear viscosity: η=999.9999999999998 
```
You can define a tuple with phase information like this:

```julia
julia> CharDim      = GEO_units(length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas);
julia> MatParam     = ( SetMaterialParams(Name="Viscous Matrix", Phase=1, 
                                     Density   = ConstantDensity(),
                                     CreepLaws = LinearViscous(η=1e23Pa*s), CharDim=CharDim),
                        SetMaterialParams(Name="Viscous Sinker", Phase=2, 
                                     Density   = PT_Density(),
                                     CreepLaws = LinearViscous(η=1e21Pa*s), CharDim=CharDim)
                        );
julia> MatParam
Phase 1 : Viscous Matrix
        | [non-dimensional units]
        | 
        |-- Density           : Constant density: ρ=2.8999999999999996e-18 
        |-- Gravity           : Gravitational acceleration: g=9.81e20 
        |-- CreepLaws         : Linear viscosity: η=999.9999999999998 
Phase 2 : Viscous Sinker
        | [non-dimensional units]
        | 
        |-- Density           : P/T-dependent density: ρ0=2.8999999999999996e-18, α=0.038194500000000006, β=0.01, T0=0.21454659702313156, P0=0.0 
        |-- Gravity           : Gravitational acceleration: g=9.81e20 
        |-- CreepLaws         : Linear viscosity: η=9.999999999999998                              
```

#### 2.2 Nonlinear creep laws
You can add pre-defined non-linear creep laws as:
```julia
julia> Phase = SetMaterialParams(Name="Viscous Matrix", Phase=2, 
                                                  Density   = ConstantDensity(),
                                                  CreepLaws = (SetDislocationCreep("Wet Olivine | Hirth & Kohlstedt (2003)"),
                                                              LinearViscous(η=1e23Pa*s)) )
Phase 2 : Viscous Matrix
        | [dimensional units]
        | 
        |-- Density           : Constant density: ρ=2900.0 kg m⁻³·⁰ 
        |-- Gravity           : Gravitational acceleration: g=9.81 m s⁻²·⁰ 
        |-- CreepLaws         : DislocationCreep: Name = Wet Olivine | Hirth & Kohlstedt (2003), n=3.5, r=1.2, A=90.0, E=480.0, V=1.1e-5, Apparatus=1 
        |                       Linear viscosity: η=1.0e23 
```
Note that the dictionary `DislocationCreep_info` has all pre-defined creep laws, so for an overview type:
```julia
julia> DislocationCreep_info
Dict{String, DislocationCreep} with 2 entries:
  "Dry Olivine | Hirth & Kohlstedt (2003)" => DislocationCreep: n=3.05, r=0, A=110000.0 MPa⁻³·⁰⁵ s⁻¹, E…
  "Wet Olivine | Hirth & Kohlstedt (2003)" => DislocationCreep: n=3.5, r=1.2, A=90 MPa⁻³·⁵ s⁻¹, E=480 k…
```

### 3. Plotting and output

A typical geodynamic simulation involves a lot of parameters. Creating data tables for scientific publications that describe all parameters employed is usually done by hand (and no-one really likes doing that). In our experience a lot of errors happen while doing this, either because the units are mixed up (some creep laws have weird units like MPa^{-n}), or because some parameters are forgotten. To help with this, we provide number of functions that 
  1)  Simplify creating plots in the same manner as in many publications that report the laboratory experiments used to create the creep laws. In that way, they can be directly compared to the original results. You can also create publication-ready figures.
  2)  Provide tools to automatically generate data tables from the input parameters. This saves time and minimizes errors. 
#### 3.1 Plotting 
A few simple functions are provided to plot various parameters.
For example, in order to plot a melting parameterisation, do:
```julia
julia> using GeoParams, Plots
Adding plotting routines of GeoParams
julia> p=MeltingParam_4thOrder();
julia> PlotMeltFraction(p);
```

#### 3.2 Automatically create data tables
When writing scientific papers that describes numerical modelling results, it is usually necessary to include tables that lists all model parameters employed. Doing this is error-prone and usually not a very interesting job to do. 
That is why we provide routines that fully automatize this process:
First, we need to define a phase.
```julia
julia> MatParam = (SetMaterialParams(Name="Viscous Matrix", Phase=1, Density=ConstantDensity(),CreepLaws = SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)")),
                   SetMaterialParams(Name="Viscous Sinker", Phase=2, Density= PT_Density(),CreepLaws = LinearViscous(η=1e21Pa*s)),
                   SetMaterialParams(Name="Viscous Bottom", Phase=3, Density= PT_Density(),CreepLaws = SetDislocationCreep("Diabase | Caristan (1982)")))
```
Next, you can create a LaTeX table for the defined phase ...
```julia
julia> ParameterTable(MatParam, filename="ParameterTable", format="latex", rdigits=3)
```
or a Markdown table.
```julia
julia> ParameterTable(MatParam, filename="ParameterTable", format="markdown", rdigits=3)
```

### 4. Computational engine
Once you have implemented a parameter in `GeoParams`, we provide allocation-free computational routines which can be called within your solver in the following manner:
```julia
args= (;T=Arrays.T_K, P=Arrays.P)
compute_density!(Rho, MatParam, Phases, args)
```
Here ```Rho``` is an array with densities, `MatParam` a tuple with material parameters as explained in [section #2](#2-material-parameters), `Phases` an array with integers which indicates which phase is present at every point, and `args` contains arguments to compute density (here: temperature and pressure). This computational routine works on the CPU, but also on the GPU (in combination with [ParallelStencil.jl](https://github.com/omlins/ParallelStencil.jl)).
Using a constant density, for example, can be specified with:
```julia
MatParam = (SetMaterialParams(Name="Crust", Phase=0, 
                Density   = ConstantDensity(ρ=2900kg/m^3)),) 
```
Using a density that employs a phase diagram (which depends on pressure and temperature) can be invoked with:
```julia
MatParam = (SetMaterialParams(Name="Mantle", Phase=0, 
                Density   = PerpleX_LaMEM_Diagram("test_data/Peridotite.in"), );
```

Importantly, you *do not have to change your code* if you want to use a new density parameterisation, as implementing this in `GeoParams` is sufficient. 
This makes it easy to use identical material parameters in a range of codes and eliminates the risk for making mistakes. It also saves a substantial amount of time in developing new codes as, in our experience, much of the debugging time is devoted to fixing bugs in the material parameters. 

### 5. Installation
You can install this package by specifying 
```julia
julia> ]
pkg> add GeoParams
```
and test whether it works with
```
pkg> test GeoParams
```

### 6. Documentation
The online documentation can be accessed [here](https://juliageodynamics.github.io/GeoParams.jl/dev/) or by clicking the blue button at the top of this page.

Sometimes it is also helpful to have a look at how we call routines in the [test](https://github.com/JuliaGeodynamics/GeoParams.jl/tree/main/test) directory.

### 7. Dependencies
The key packages we rely on:
- [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) to deal with SI units
- [Parameters.jl](https://github.com/mauro3/Parameters.jl) to have structures that are easier to modify
- [LaTeXStrings.jl](https://github.com/stevengj/LaTeXStrings.jl) to be able to add equations to the structures that describe the employed material laws
### 8. Contributing
Help with developing this package is highly appreciated. You can contribute for example by adding new creep laws or by adding new constitutive relationships. If you invest a bit of time now, it will save others in the community a lot of time! 
The simplest way to do this is by cloning the repository, and creating a new branch for your feature. Once you are happy with what you added (and after you added a test to ensure that it will keep working with future changes), create a pull request and we will evaluate & merge it.


### 9. Funding
The development of this package was supported by the European Research Council (ERC CoG #771143 MAGMA).
