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
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractCreepLaw{T} <: AbstractMaterialParam end

export computeCreepLaw_EpsII,
    computeCreepLaw_TauII,       # calculation routines
    CreepLawVariables,                                     # holds additional parameters required for calculations
    LinearViscous,
    PowerlawViscous
param_info

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

@with_kw struct CreepLawVariables{_T,U1,U2,U3,U4}
    P::GeoUnit{_T,U1} = 100.0MPa   # pressure
    T::GeoUnit{_T,U2} = 500.0C     # temperature
    d::GeoUnit{_T,U3} = 1.0cm      # grainsize
    f::GeoUnit{_T,U4} = 0.0MPa     # water-fugacity         
end
CreepLawVariables(args...) = CreepLawVariables(convert.(GeoUnit, args)...)

include("DislocationCreep.jl")
include("DiffusionCreep.jl")

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
@with_kw_noshow struct LinearViscous{T,U} <: AbstractCreepLaw{T}
    η::GeoUnit{T,U} = 1e20Pa * s                # viscosity
    η_val::T = 1.0
end
LinearViscous(args...) = LinearViscous(convert(GeoUnit, args[1]), args[2])

function param_info(s::LinearViscous) # info about the struct
    return MaterialParamsInfo(; Equation=L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}")
end

# Calculation routines for linear viscous rheologies
@inline function computeCreepLaw_EpsII(TauII::T, s::LinearViscous; kwargs...) where {T}
    η = s.η

    return EpsII = (TauII / η) * 0.5
end

@inline function dεII_dτII(TauII, s::LinearViscous; kwargs...)
    η = s.η

    return η * 0.5
end

@inline function computeCreepLaw_TauII(EpsII, s::LinearViscous; kwargs...)
    η = s.η

    return TauII = 2 * (η * EpsII)
end

@inline function dτII_dεII(EpsII, s::LinearViscous; kwargs...)
    η = s.η

    return 2 * η
end
# Print info 
function show(io::IO, g::LinearViscous)
    return print(io, "Linear viscosity: η=$(g.η.val)")
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
@with_kw_noshow struct PowerlawViscous{T,U1,U2,U3} <: AbstractCreepLaw{T}
    η0::GeoUnit{T,U1} = 1e18Pa * s            # reference viscosity 
    n::GeoUnit{T,U2} = 2.0 * NoUnits         # powerlaw exponent
    ε0::GeoUnit{T,U3} = 1e-15 * 1 / s          # reference strainrate
end
PowerlawViscous(a...) = PowerlawViscous(convert.(GeoUnit, a)...)

# Print info 
function show(io::IO, g::PowerlawViscous)
    return print(io, "Powerlaw viscosity: η0=$(g.η0.val), n=$(g.n.val), ε0=$(g.ε0.val) ")
end
#-------------------------------------------------------------------------

# add methods programatically 
for myType in (:LinearViscous, :DiffusionCreep, :DislocationCreep)
    @eval begin
        function computeCreepLaw_EpsII(TauII, a::$(myType), args)
            return computeCreepLaw_EpsII(TauII, a; args...)
        end
        dεII_dτII(TauII, a::$(myType), args) = dεII_dτII(TauII, a; args...)
        function computeCreepLaw_TauII(EpsII, a::$(myType), args)
            return computeCreepLaw_TauII(EpsII, a; args...)
        end
        if Symbol($myType) !== :ConstantElasticity
            dτII_dεII(EpsII, a::$(myType), args) = dτII_dεII(EpsII, a; args...)
        end
    end
end

# Help info for the calculation routines
"""
    computeCreepLaw_EpsII(TauII, s:<AbstractCreepLaw, p::CreepLawVariables)

Returns the strainrate invariant ``\\dot{\\varepsilon}_{II}`` for a given deviatoric stress 
invariant ``\\tau_{II}`` for any of the viscous creep laws implemented.
``p`` is a struct that can hold various parameters (such as ``P``, ``T``) that the creep law
may need for the calculations 

```math  
    \\dot{\\varepsilon}_{II}   = f( \\tau_{II} ) 
```

"""
computeCreepLaw_EpsII

"""
    computeCreepLaw_TauII(EpsII, s:<AbstractCreepLaw, p::CreepLawVariables)

Returns the deviatoric stress invariant ``\\tau_{II}`` for a given strain rate  
invariant ``\\dot{\\varepsilon}_{II}`` for any of the viscous creep laws implemented.
``p`` is a struct that can hold various parameters (such as ``P``, ``T``) that the creep law
may need for the calculations 

```math  
    \\tau_{II}  = f( \\dot{\\varepsilon}_{II} ) 
```

"""
computeCreepLaw_TauII

end
