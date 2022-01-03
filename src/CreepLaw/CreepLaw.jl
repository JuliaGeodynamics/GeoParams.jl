module CreepLaw

# This implements viscous creep laws and routines to compute with them
#
# Note that various simple creep laws are defined in this file; 
# more complex ones (such as DislocationCreep) are in separate files 
# in the same directory
#
# In case you want to add new creep laws, have a look at how the ones
# here are implemented. Please add tests as well!

using Base: Float64
using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
using BibTeX
import Base.show

abstract type AbstractCreepLaw{T} <: AbstractMaterialParam end

export  ComputeCreepLaw_EpsII, ComputeCreepLaw_TauII,       # calculation routines
        CreepLawVariables,                                     # holds additional parameters required for calculations
        LinearViscous, 
        PowerlawViscous



# NOTE: we will likely have to remove this, in favor of multiple dispatch options
"""
 Struct that holds parameters required for creep law calculations (such as P,T, grain size etc.)

You can assign the struct as
```julia-repl
 p = CreepLawVariables(P=0.0, T=0.0, d=1.0) 
```  
where you can also pass vectors or arrays as values.

Note that, if needed, this can be extended, w/out interfering with existing calculation  
"""
@with_kw struct CreepLawVariables     
    P::GeoUnit  =   100.0MPa   # pressure
    T::GeoUnit  =   500.0C     # temperature
    d::GeoUnit  =   1.0cm      # grainsize
    f::GeoUnit  =   0.0MPa     # water-fugacity         
end

include("DislocationCreep.jl")
# Linear viscous rheology ------------------------------------------------
"""
    LinearViscous(η=1e20Pa*s)
    
Defines a linear viscous creeplaw 

The (isotopic) linear viscous rheology is given by  
```math  
    \\tau_{ij} = 2 \\eta \\dot{\\varepsilon}_{ij} 
```
or
```math  
    \\dot{\\varepsilon}_{ij}  = {\\tau_{ij}  \\over 2 \\eta }
```

where ``\\eta_0`` is the reference viscosity [Pa*s] at reference strain rate ``\\dot{\\varepsilon}_0``[1/s], and ``n`` the power law exponent []. 
"""
@with_kw_noshow struct LinearViscous{T} <: AbstractCreepLaw{T}
   # equation::LaTeXString   =   L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}"     
    η::GeoUnit{T}       =   1e20Pa*s                # viscosity
    η_val::T            =   1.0
end
LinearViscous(a...) = LinearViscous{Float64}(a...)

# Calculation routines for linear viscous rheologies
function ComputeCreepLaw_EpsII(TauII, s::LinearViscous, p::CreepLawVariables)
    @unpack η   = s
    
    return EpsII = (TauII/η)*0.5;
end

function ComputeCreepLaw_TauII(EpsII, s::LinearViscous, p::CreepLawVariables)
    @unpack η   = s
    
    return TauII = 2.0*(η*EpsII);
end

# Print info 
function show(io::IO, g::LinearViscous)  
    print(io, "Linear viscosity: η=$(g.η.val)")  
end
#-------------------------------------------------------------------------

# Powerlaw viscous rheology ----------------------------------------------
"""
    PowerlawViscous(η0=1e18Pa*s, n=2.0NoUnits, ε0=1e-15/s)
    
Defines a power law viscous creeplaw as: 

```math  
    \\tau_{ij}^n  = 2 \\eta_0 \\left( \\dot{\\varepsilon}_{ij} \\over \\dot{\\varepsilon}_0 \\right)
```
where ``\\eta`` is the effective viscosity [Pa*s].
"""
@with_kw_noshow mutable struct PowerlawViscous{T} <: AbstractCreepLaw{T}
    η0::GeoUnit{T}     =  1e18Pa*s            # reference viscosity 
    n::GeoUnit{T}      =  2.0*NoUnits         # powerlaw exponent
    ε0::GeoUnit{T}     =  1e-15*1/s;          # reference strainrate
end
PowerlawViscous(a...) = PowerlawViscous{Float64}(a...)

# Print info 
function show(io::IO, g::PowerlawViscous)  
    print(io, "Powerlaw viscosity: η0=$(g.η0.val), n=$(g.n.val), ε0=$(g.ε0.val) ")
end
#-------------------------------------------------------------------------


# Help info for the calculation routines
"""
    ComputeCreepLaw_EpsII(TauII, s:<AbstractCreepLaw, p::CreepLawVariables)

Returns the strainrate invariant ``\\dot{\\varepsilon}_{II}`` for a given deviatoric stress 
invariant ``\\tau_{II}`` for any of the viscous creep laws implemented.
``p`` is a struct that can hold various parameters (such as ``P``, ``T``) that the creep law
may need for the calculations 

```math  
    \\dot{\\varepsilon}_{II}   = f( \\tau_{II} ) 
```

"""
ComputeCreepLaw_EpsII

"""
    ComputeCreepLaw_TauII(EpsII, s:<AbstractCreepLaw, p::CreepLawVariables)

Returns the deviatoric stress invariant ``\\tau_{II}`` for a given strain rate  
invariant ``\\dot{\\varepsilon}_{II}`` for any of the viscous creep laws implemented.
``p`` is a struct that can hold various parameters (such as ``P``, ``T``) that the creep law
may need for the calculations 

```math  
    \\tau_{II}  = f( \\dot{\\varepsilon}_{II} ) 
```

"""
ComputeCreepLaw_TauII


end