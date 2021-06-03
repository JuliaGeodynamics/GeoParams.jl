"""
    This creates the overall material parameters structure
"""

module MaterialParameters
using Unitful
using Parameters
using ..Units

import Base.show
using GeoParams: AbstractMaterialParam

export 
    MaterialParams
    
"""
    MaterialParams
    
Structure that holds all material parameters for a given phase

"""
 @with_kw_noshow mutable struct MaterialParams
    # 
    Name::String        =   ""                  #       Description/name of the phase
    Density             =   nothing             #       Density equation of state
    CreepLaws::Tuple    =   nothing             #       Creep laws
    Elasticity          =   nothing             #       Elastic parameters
    Plasticity          =   nothing             #       Plasticity
    Conductivity        =   nothing             #       Parameters related to the energy equation 
    HeatCapacity        =   nothing             #        
    EnergySourceTerms   =   nothing             #       Source terms in energy conservation equation
end

# Link the various definitions
include("./CreepLaw/CreepLaw.jl")


end


