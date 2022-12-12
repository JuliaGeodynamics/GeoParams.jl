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
using StaticArrays
using LinearAlgebra

import Base: getindex

# overload to account for cases where this is an integer
for T in (:Real, :Symbol)
    @eval begin
        Base.getindex(val::$(T), I::Vararg{Integer,N}) where {N} = val
        Base.getindex(val::$(T), I::Integer) = val
    end
end

export @u_str,
    uconvert,
    upreffered,
    unit,
    ustrip,
    NoUnits,  #  Units 
    GeoUnit,
    GeoUnits,
    GEO_units,
    SI_units,
    NO_units,
    AbstractGeoUnit,
    nondimensionalize,
    dimensionalize,
    superscript,
    upreferred,
    GEO,
    SI,
    NONE,
    isDimensional,
    Value,
    NumValue,
    Unit,
    UnitValue,
    isdimensional,
    km,
    m,
    cm,
    mm,
    μm,
    Myrs,
    yr,
    s,
    MPa,
    Pa,
    kbar,
    Pas,
    K,
    C,
    g,
    kg,
    mol,
    J,
    kJ,
    Watt,
    μW,
    Quantity

export AbstractGeoUnit1, GeoUnit1

#         
abstract type AbstractMaterialParam end                                    # structure that holds material parmeters (density, elasticity, viscosity)          
abstract type AbstractMaterialParamsStruct end                             # will hold all info for a phase       
abstract type AbstractPhaseDiagramsStruct <: AbstractMaterialParam end    # will hold all info for phase diagrams 
abstract type AbstractConstitutiveLaw{T} <: AbstractMaterialParam end
abstract type AbstractComposite <: AbstractMaterialParam end

function PerpleX_LaMEM_Diagram end                                         # necessary as we already use this function in Units, but only define it later in PhaseDiagrams
function param_info end
export AbstractMaterialParam, AbstractMaterialParamsStruct, AbstractPhaseDiagramsStruct

include("Utils.jl")

include("TensorAlgebra/TensorAlgebra.jl")
export second_invariant, second_invariant_staggered, rotate_elastic_stress

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

# Define Table output functions
include("Tables.jl")
using .Tables
export Phase2Dict, Dict2LatexTable

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

# Constitutive relationships laws
using .MaterialParameters.ConstitutiveRelationships
export AxialCompression, SimpleShear, Invariant 

#       Calculation routines
export dεII_dτII,
    dτII_dεII,
    dεII_dτII_AD,
    dτII_dεII_AD,
    dεvol_dp,
    dp_dεvol,
    compute_εII!,
    compute_εII,
    compute_εII_AD,
    compute_τII!,
    compute_τII,
    compute_τII_AD,
    compute_εvol!,
    compute_εvol,
    compute_p!,
    compute_p,
    CorrectionFactor,
    remove_tensor_correction,
    isvolumetric,

    #       Viscous creep laws
    AbstractCreepLaw,
    LinearViscous,
    PowerlawViscous,
    ArrheniusType,
    DislocationCreep,
    SetDislocationCreep,
    DiffusionCreep,
    SetDiffusionCreep,
    DislocationCreep_info,
    DiffusionCreep_info,

    #       Elasticity
    AbstractElasticity,
    ConstantElasticity,
    SetConstantElasticity,
    effective_εII,

    #       Plasticity
    AbstractPlasticity,
    compute_yieldfunction,
    compute_yieldfunction!,
    DruckerPrager,
    DruckerPrager_regularised,
    compute_plasticpotentialDerivative,
    ∂Q∂τ,
    ∂Q∂P,∂Q∂τII,
    ∂F∂τII,∂F∂P,∂F∂λ,
    
    #       Composite rheologies
    AbstractConstitutiveLaw,
    AbstractComposite,
    computeViscosity_τII,
    computeViscosity_εII,
    computeViscosity_τII!,
    computeViscosity_εII!,
    computeViscosity_εII_AD,
    local_iterations_εII,    
    local_iterations_εII_AD,
    local_iterations_τII,
    local_iterations_τII_AD,
    computeViscosity,
    InverseCreepLaw,
    KelvinVoigt,
    CompositeRheology,
    Parallel,
    create_rheology_string, print_rheology_matrix,
    compute_εII_harmonic, compute_τII_AD,
    isplastic,isvolumetricplastic,
    compute_p_τII, 
    local_iterations_εvol, 
    compute_p_harmonic
    

# Constitutive relationships laws
include("StressComputations/StressComputations.jl")
export compute_τij, compute_p_τij, compute_τij_stagcenter!, compute_p_τij_stagcenter!, compute_τij!, compute_p_τij!

include("Rheology_Utils.jl")
export time_τII_0D, time_τII_0D!, time_p_τII_0D, time_p_τII_0D!


# Gravitational Acceleration
using .MaterialParameters.GravitationalAcceleration
export compute_gravity,                                # computational routines
    ConstantGravity

# Energy parameters: Heat Capacity, Thermal conductivity, latent heat, radioactive heat         
using .MaterialParameters.HeatCapacity
export compute_heatcapacity,
    compute_heatcapacity!, ConstantHeatCapacity, T_HeatCapacity_Whittington

using .MaterialParameters.Conductivity
export compute_conductivity,
    compute_conductivity!,
    ConstantConductivity,
    T_Conductivity_Whittington,
    T_Conductivity_Whittington_parameterised,
    TP_Conductivity,
    Set_TP_Conductivity

using .MaterialParameters.LatentHeat
export compute_latent_heat, compute_latent_heat!, ConstantLatentHeat

using .MaterialParameters.RadioactiveHeat
export compute_radioactive_heat,
    compute_radioactive_heat!, ConstantRadioactiveHeat, ExpDepthDependentRadioactiveHeat

using .MaterialParameters.Shearheating
export compute_shearheating!, compute_shearheating, ConstantShearheating

# Add TAS classification
include("./RockClassification/TASclassification.jl")
using .TASclassification
export TASclassificationData, computeTASclassification, retrieveTASrockType

# Add zircon saturation parameterizations
include("./ZirconAge/ZirconAges.jl")
using .ZirconAges
export ZirconAgeData,
    compute_zircon_age_PDF,
    compute_zircons_Ttpath,
    zircon_age_PDF,
    compute_zircons_convert_vecs2mat

# Seismic velocities
using .MaterialParameters.SeismicVelocity
export compute_pwave_velocity,
    compute_wave_velocity,
    compute_wave_velocity!,
    ConstantSeismicVelocity,
    anelastic_correction,
    melt_correction,
    porosity_correction,
    correct_wavevelocities_phasediagrams,
    melt_correction_Takei

# Add melting parameterizations
include("./MeltFraction/MeltingParameterization.jl")
using .MeltingParam
export compute_meltfraction,
    compute_meltfraction!,       # calculation routines
    compute_dϕdT,
    compute_dϕdT!,
    MeltingParam_Caricchi,
    MeltingParam_4thOrder,
    MeltingParam_5thOrder,
    MeltingParam_Quadratic,
    MeltingParam_Assimilation,
    SmoothMelting

# Add 1D Strength Envelope
include("./StrengthEnvelope/StrengthEnvelope.jl")

# Add plotting routines - only activated if the "GLMakie.jl" package is loaded 
function __init__()
    @require GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a" begin
        print("Adding plotting routines of GeoParams through GLMakie \n")
        @eval include("Plotting/Plotting.jl")
        @eval include("Plotting/StrengthEnvelope.jl")
    end
end

#Set functions aliases using @use
include("aliases.jl")
export ntuple_idx

end # module
