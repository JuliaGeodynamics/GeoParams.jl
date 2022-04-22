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
        GeoUnit, GeoUnits, GEO_units, SI_units, NO_units, AbstractGeoUnit,
        nondimensionalize, dimensionalize,
        superscript, upreferred, GEO, SI, NONE, isDimensional, Value, NumValue, Unit, UnitValue, isdimensional,
        km, m, cm, mm, μm, Myrs, yr, s, MPa, Pa, kbar, Pas, K, C, g, kg, mol, J, kJ, Watt, μW, Quantity

export AbstractGeoUnit1, GeoUnit1

#         
abstract type AbstractMaterialParam end                                    # structure that holds material parmeters (density, elasticity, viscosity)          
abstract type AbstractMaterialParamsStruct end                             # will hold all info for a phase       
abstract type AbstractPhaseDiagramsStruct <: AbstractMaterialParam end    # will hold all info for phase diagrams 
function PerpleX_LaMEM_Diagram end                                         # necessary as we already use this function in Units, but only define it later in PhaseDiagrams
function param_info end
export AbstractMaterialParam, AbstractMaterialParamsStruct, AbstractPhaseDiagramsStruct


# note that this throws a "Method definition warning regarding superscript"; that is expected & safe 
#  as we add a nicer way to create output of superscripts. I have been unable to get rid of this warning,
#  as I am indeed redefining a method originally defined in Unitful
include("Units.jl")
using .Units
export @unpack_units, @unpack_val
export compute_units

# Define Material Parameter structure
include("MaterialParameters.jl")
using .MaterialParameters
export MaterialParams, SetMaterialParams, No_MaterialParam, MaterialParamsInfo

# Phase Diagrams
using .MaterialParameters.PhaseDiagrams
export PhaseDiagram_LookupTable, PerpleX_LaMEM_Diagram

# Density
using .MaterialParameters.Density
export compute_density,                                # computational routines
        compute_density!,
        param_info,
        AbstractDensity,
        No_Density,
        ConstantDensity,
        PT_Density,
        Compressible_Density,
        PhaseDiagram_LookupTable,
        Read_LaMEM_Perple_X_Diagram

# Creep laws
using .MaterialParameters.CreepLaw
export computeCreepLaw_EpsII, computeCreepLaw_TauII, CreepLawVariables,
        LinearViscous, PowerlawViscous,
        DislocationCreep, SetDislocationCreep,
        DiffusionCreep,   SetDiffusionCreep

# Gravitational Acceleration
using .MaterialParameters.GravitationalAcceleration
export compute_gravity,                                # computational routines
        ConstantGravity


# Energy parameters: Heat Capacity, Thermal conductivity, latent heat, radioactive heat         
using .MaterialParameters.HeatCapacity
export compute_heatcapacity,
        compute_heatcapacity!,
        ConstantHeatCapacity,
        T_HeatCapacity_Whittington

using .MaterialParameters.Conductivity
export compute_conductivity,
        compute_conductivity!,
        ConstantConductivity,
        T_Conductivity_Whittington,
        T_Conductivity_Whittington_parameterised,
        TP_Conductivity,
        Set_TP_Conductivity

using .MaterialParameters.LatentHeat
export compute_latent_heat,compute_latent_heat!,
        ConstantLatentHeat

using .MaterialParameters.RadioactiveHeat
export compute_radioactive_heat,compute_radioactive_heat!,
        ConstantRadioactiveHeat,
        ExpDepthDependentRadioactiveHeat

using .MaterialParameters.Shearheating
export compute_shearheating!, compute_shearheating,
        ConstantShearheating

# Add zircon saturation parameterizations
include("./ZirconAge/ZirconAges.jl")
using   .ZirconAges
export  ZirconAgeData, 
        compute_zircon_age_PDF,  compute_zircons_Ttpath, 
        zircon_age_PDF, compute_zircons_convert_vecs2mat 

# Seismic velocities
using .MaterialParameters.SeismicVelocity
export compute_pwave_velocity,          compute_swave_velocity,
        compute_pwave_velocity!,        compute_swave_velocity!,
        ConstantSeismicVelocity

# Add melting parameterizations
include("./MeltFraction/MeltingParameterization.jl")
using .MeltingParam
export  compute_meltfraction,   compute_meltfraction!,       # calculation routines
        compute_dϕdT,           compute_dϕdT!,
        MeltingParam_Caricchi,  MeltingParam_4thOrder, 
        MeltingParam_5thOrder,  MeltingParam_Quadratic,
        MeltingParam_Assimilation


# Add plotting routines - only activated if the "Plots.jl" package is loaded 
function __init__()
        @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
                print("Adding plotting routines of GeoParams")
                @eval include("./Plotting.jl")
        end
end

#Set functions aliases using @use
include("aliases.jl")

# export ntuple_idx()
include("Utils.jl")
export ntuple_idx

end # module
