"""
Typical geodynamic simulations involve a large number of material parameters that have units that are often inconvenient to be directly used in numerical models
This package has two main features that help with this:
- Create a nondimensionalization object, which can be used to transfer dimensional to non-dimensional parameters (usually better for numerical solvers)
- Create an object in which you can specify material parameters employed in the geodynamic simulations

The material parameter object is designed to be extensible and can be passed on to the solvers, such that new creep laws or features can be readily added. 
We also implement some typically used creep law parameters, together with tools to plot them versus and compare our results with those of published papers (to minimize mistakes). 
"""
__precompile__()
module GeoParams

using Parameters        # helps setting default parameters in structures
using Unitful

export 
        @u_str, uconvert, upreffered, unit, ustrip, NoUnits,  #  Units 
        GEO_units, SI_units, NO_units, AbstractGeoUnits, 
        Nondimensionalize, Dimensionalize, superscript, upreferred, 
        km, m, cm, Mtrs, yr, s, MPa, Pa, Pas, K, C, kg, mol
    
# note that this throws a "Method definition warning regarding superscript"; that is expected & safe 
#  as we add a nicer way to create output of superscripts. I have been unable to suppress this warning
include("Units.jl")     
using .Units

end # module
