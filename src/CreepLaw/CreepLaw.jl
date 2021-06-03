module CreepLaw

# This implements viscous creep laws and routines to compute with them
#
# Note that various simple creep laws are defined in this file; 
# more complex ones (such as DislocationCreep) are in separate files 
# in the same directory
#
# In case you want to add new creep laws, have a look at how the ones
# here are implemented. Please add tests as well!

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam

abstract type AbstractCreepLaw <: AbstractMaterialParam end

export  CreepLaw_EpsII, CreepLaw_TauII,         # calculation routines
        LinearViscous, 
        PowerlawViscous

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

# Calculation routines for linear viscous rheologies
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



# Help info for the calculation routines
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