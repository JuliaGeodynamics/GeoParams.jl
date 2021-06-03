"""
    This creates the overall material parameters structure
"""

module MaterialParameters
using Unitful
using Parameters
using LaTeXStrings
using GeoParams.Units

import Base.show
using GeoParams: AbstractMaterialParam

export 
    MaterialParams, Nondimensionalize!,
    CreepLaw_EpsII, CreepLaw_TauII,
    LinearViscous
    
"""
    MaterialParams
    
Structure that holds all material parameters for a given phase

"""
 @with_kw_noshow mutable struct MaterialParams
    # 
    Name::String        =   ""                  #       Description/name of the phase
    Density             =   nothing             #       Density equation of state
    CreepLaws           =   nothing             #       Creep laws
    Elasticity          =   nothing             #       Elastic parameters
    Plasticity          =   nothing             #       Plasticity
    Conductivity        =   nothing             #       Parameters related to the energy equation 
    HeatCapacity        =   nothing             #        
    EnergySourceTerms   =   nothing             #       Source terms in energy conservation equation
end


# just for testing: define a few creep laws here (will be put into a different file later )
abstract type AbstractCreepLaw <: AbstractMaterialParam end

# Linear viscous rheology ------------------------------------------------
"""
    LinearViscous(eta=1e20Pa*s)
    
Defines a linear viscous creeplaw 

The (isotopic) linear viscous rheology is given by  
```math  
    \\dot{\\varepsilon}_{ij}  = {\\tau_{ij}  \\over 2 \\eta }
```
where ``\\eta`` is the effective viscosity [Pa*s].
"""
@with_kw_noshow mutable struct LinearViscous <: AbstractCreepLaw
    equation::LaTeXString   =   L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"     
    eta::GeoUnit            =   1e20Pa*s                # viscosity
end

function CreepLaw_EpsII(TauII, s::LinearViscous)
    @unpack eta   = s
    
    return EpsII = TauII/(2.0*eta);
end

function CreepLaw_TauII(EpsII, s::LinearViscous)
    @unpack eta   = s
    
    return TauII = 2.0*eta*EpsII;
end
#-------------------------------------------------------------------------

# Powerlaw viscous rheology ----------------------------------------------
"""
    PowerlawViscous(eta=1e20Pa*s, )
    
Defines a simple power law viscous creeplaw 

The (isotopic) linear viscous rheology is given by  
```math  
    {\\tau_{ij}  = \\η_0 ( \\varepsilon_{II} \\over \\varepsilon_{II} ) 
```
where ``\\eta`` is the effective viscosity [Pa*s].
"""
@with_kw_noshow mutable struct PowerlawViscous <: AbstractCreepLaw
    η0::GeoUnit     =  1e18Pa*s            # reference viscosity 
    n::GeoUnit      =  2.0*NoUnit          # powerlaw exponent
    ε0::GeoUnit     =  1e-15*1/s;          # reference strainrate
end
#-------------------------------------------------------------------------



"""
    CreepLaw_EpsII(TauII, s:<AbstractCreepLaw)

Returns the strainrate invariant ``\\dot{\\varepsilon}_{II}`` for a given deviatoric stress 
invariant ``\\tau_{II}`` for any of the viscous creep laws implemented.

"""
CreepLaw_EpsII

"""
    CreepLaw_TauII(EpsII, s:<AbstractCreepLaw)

Returns the deviatoric stress invariant ``\\tau_{II}`` for a given strain rate  
invariant ``\\dot{\\varepsilon}_{II}`` for any of the viscous creep laws implemented.

"""
CreepLaw_TauII

end


