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
using Unitful           # Units
using BibTeX            # references of creep laws

export 
        @u_str, uconvert, upreffered, unit, ustrip, NoUnits,  #  Units 
        GeoUnit, GEO_units, SI_units, NO_units, AbstractGeoUnits, 
        Nondimensionalize, Nondimensionalize!, Dimensionalize, Dimensionalize!,
        superscript, upreferred, GEO, SI, NONE, isDimensional, 
        km, m, cm, mm, Myrs, yr, s, MPa, Pa, Pas, K, C, kg, mol, J, kJ
   
#         
abstract type AbstractMaterialParam end           # structure that holds material parmeters (density, elasticity, viscosity)          
abstract type AbstractMaterialParamsStruct end    # will hold all info for a phase       
        
# note that this throws a "Method definition warning regarding superscript"; that is expected & safe 
#  as we add a nicer way to create output of superscripts. I have been unable to get rid of this warning,
#  as I am indeed redefining a method originally defined in Unitful
include("Units.jl")     
using .Units

# Define Material Parameter structure
include("MaterialParameters.jl")
using  .MaterialParameters
export MaterialParams, SetMaterialParams

# Creep laws
using  .MaterialParameters.CreepLaw
export  ComputeCreepLaw_EpsII, ComputeCreepLaw_TauII, CreepLawVariables,
        LinearViscous, PowerlawViscous, 
        DislocationCreep, SetDislocationCreep

# Density
using  .MaterialParameters.Density
export  ComputeDensity,                                # computational routines
        ConstantDensity,                        
        PT_Density

# Gravitational Acceleration
using  .MaterialParameters.GravitationalAcceleration
export  ComputeGravity,                                # computational routines
        ConstantGravity


# Add plotting routines
include("Plotting.jl")
using  .Plotting
export  PlotStressStrainrate_CreepLaw

end # module
