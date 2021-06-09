# GeoParams.jl

Typical geodynamic simulations involve a large number of material parameters that have units that are often inconvenient to be directly used in numerical models This package has two main features that help with this:

- Nondimensionalization object, which can be used to transfer dimensional to non-dimensional parameters (usually better for numerical solvers).
- Material parameters object in which you can specify  parameters employed in the geodynamic simulations. This object is designed to be extensible and can be passed on to the solvers, such that new creep laws or features can be readily added. 
  
We also implement some typically used creep law parameters, together with tools to plot them versus and compare our results with those of published papers (to minimize mistakes).
