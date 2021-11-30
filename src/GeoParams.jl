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
using Requires          # To only add plotting routines if Plots is loaded

export 
        @u_str, uconvert, upreffered, unit, ustrip, NoUnits,  #  Units 
        GeoUnit, GEO_units, SI_units, NO_units, AbstractGeoUnits, 
        Nondimensionalize, Nondimensionalize!, Dimensionalize, Dimensionalize!,
        superscript, upreferred, GEO, SI, NONE, isDimensional, Value, Unit, 
        km, m, cm, mm, Myrs, yr, s, MPa, Pa, Pas, K, C, g, kg, mol, J, kJ, Watt, μW
   
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
        DislocationCreep, SetDislocationCreep,
	DiffusionCreep, SetDiffusionCreep

# Density
using  .MaterialParameters.Density
export  ComputeDensity,                                # computational routines
        ConstantDensity,                        
        PT_Density

# Gravitational Acceleration
using  .MaterialParameters.GravitationalAcceleration
export  ComputeGravity,                                # computational routines
        ConstantGravity

# Energy parameters: Heat Capacity, Thermal conductivity, latent heat, radioactive heat         
using .MaterialParameters.HeatCapacity
export  ComputeHeatCapacity,                          
        ConstantHeatCapacity,
        T_HeatCapacity_Whittacker

using .MaterialParameters.Conductivity
export  ComputeConductivity,                           
        ConstantConductivity,
        T_Conductivity_Whittacker,
        TP_Conductivity,
        Set_TP_Conductivity

using .MaterialParameters.LatentHeat
export  ComputeLatentHeat,                           
        ConstantLatentHeat
        
using .MaterialParameters.RadioactiveHeat        
export  ComputeRadioactiveHeat,                 
        ConstantRadioactiveHeat                  

using .MaterialParameters.Shearheating        
export  ComputeShearheating, ComputeShearheating!,               
        ConstantShearheating              

# Add melting parameterizations
include("./MeltFraction/MeltingParameterization.jl")
using .MeltingParam
export  ComputeMeltingParam, ComputeMeltingParam!,       # calculation routines
        MeltingParam_Caricchi                          

# Add plotting routines
function __init__()
        @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
                print("Adding plotting routines of GeoParam")
                @eval include("./Plotting.jl")
        end
end


end # module
