# GeoParams.jl
Typical geodynamic simulations involve a large number of material parameters, which are often based on laboratory experiments.
This goal of this package is to facilitate this by:
- Creating a nondimensionalization object, which can be used to transfer dimensional to non-dimensional parameters (usually better for numerical solvers)
- Creating an object in which you can specify material parameters employed in geodynamic simulations

The material parameter object is designed to be extensible and can be passed on to the solvers, such that new creep laws can be readily added. 
We also implement some typically used creep law parameters, together with tools to plot them versus and compare our results with those of published papers (to minimize mistakes).  

### 1. Nondimensionalization 
Typical geodynamic simulations involve dimensions on the order of 10's-1000's of kilometers, and viscosities on the order of ~1e20 Pas. If such values are directly employed in numerical solvers, they may result in roundoff errors. It is therefore common practice to nondimensionalize the input parameters by dividing them by typical values such that the result gives numbers that are closer to one.
This can be done by specifying characteristic values for `length`, `stress`, `temperature` and `viscosity`. From these `basic` units all other physical units are derived and input parameters can thus be nondimensionalized accordingly (and dimensionalized again when plotting or saving output). 

As you learned in physics, the common approach to do this is by using `SI` units. Yet, as meters and seconds are not so convenient in geodynamics, where we usually deal with lengthscales on the orders of kilometers and timescales in millions of years, we also provide the `geo` object, which allows to give input parameters in more convenient units.

```julia
julia> using GeoParams
julia> GeoUnit = GeoUnits(type="geo",length=1000,temperature=1000,viscosity=1e20, stress=1)
```
You can use 3 `types`:
  1. *geo* units: Units of length in the code are expected to be in kilometers, time is in million of years (Myrs) and stresses are in MPa (1e6 Pa). This is the default.
  2. *SI* units: all values are in SI units (meters, Pascal, seconds)
  3. *none*: all input parameters are in nondimensional units

Once a `GeoUnit` structure is created, you can use the derived parameters  

### 2. Material parameters  
All geodynamic simulations require specifying material parameters, such as (nonlinear) viscous constitutive relationships or an equation of state for density. These parameters are usually specified per `phase`. Here, we provide a framework that simplifies doing that. Thanks to the flexibility of julia, we can actually directly embed the function that does the computations in the structure itself, which makes it straightforward to extend it and add new creep laws (which can directly be used in the solvers).  

Some examples of this is used:
#### 2.1 Constant density, constant linear viscosity
to be added

#### 2.2 Nonlinear creep laws
to be added


### 3. Plotting and output

A typical geodynamic simulation involves a lot of parameters. Creating data tables for scientific publications that describe all parameters employed is usually done by hand (and no-one really likes doing that). In our experience a lot of errors happen while doing this, either because the units are mixed up (some creep laws have weird units like MPa^{-n}). To help with this, we provide number of functions that 
  1)  Simplify creating plots in the same manner as in many publications, such that they can be directly compared to the original results. You can also create publication-ready figures.
  2)  Provide tools to automatically generate data tables from the input parameters. This saves time and minimizes errors.
#### 3.1 Plotting creep laws 
A few simple functions are provided to plot creep laws.

#### 3.2 Creating data tables
to be developed


### 4. Installation and usage
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

### 6. Contributing
Help with developing this package is highly appreciated. You can contribute for example by adding new creep laws or by adding new constitutive relationships. If you invest a bit of time now, it will save others in the community a lot of time! This

The simplest way to do this is by cloning the repository, and creating a local 
