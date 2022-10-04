# This implements viscous creep laws and routines to compute with them
#
# Note that various simple creep laws are defined in this file; 
# more complex ones (such as DislocationCreep) are in separate files 
# in the same directory
#
# In case you want to add new creep laws, have a look at how the ones
# here are implemented. Please add tests as well!

abstract type AbstractCreepLaw{T} <: AbstractConstitutiveLaw{T} end

export LinearViscous, PowerlawViscous, CorrectionFactor, AbstractCreepLaw

# This computes correction factors to go from experimental data to tensor format
function CorrectionFactor(a::AbstractCreepLaw{_T}) where {_T}
    if a.Apparatus == AxialCompression
        FT = sqrt(one(_T) * 3) # relation between differential stress recorded by apparatus and TauII
        FE = 2 / FT            # relation between gamma recorded by apparatus and EpsII
        return FT, FE

    elseif a.Apparatus == SimpleShear
        FT = one(_T) * 2      # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = FT
        return FT, FE

    elseif a.Apparatus == Invariant
        FT, FE = one(_T), one(_T)
        return FT, FE
    end
end

function CorrectionFactor(a::_T) where {_T}
    if a == AxialCompression
        FT = sqrt(one(_T) * 3) # relation between differential stress recorded by apparatus and TauII
        FE = 2 / FT            # relation between gamma recorded by apparatus and EpsII
        return FT, FE

    elseif a == SimpleShear
        FT = one(_T) * 2      # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = FT
        return FT, FE

    elseif a == Invariant
        FT, FE = one(_T), one(_T)
        return FT, FE
    end
end

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

function param_info(a::LinearViscous) # info about the struct
    return MaterialParamsInfo(; Equation=L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}")
end

# Calculation routines for linear viscous rheologies
function compute_εII(a::LinearViscous, TauII; kwargs...)
    @unpack η = a

    return (TauII / η) * 0.5
end

"""
    
    compute_εII!(EpsII::AbstractArray{_T,N}, s::LinearViscous, TauII::AbstractArray{_T,N})
"""
function compute_εII!(
    EpsII::AbstractArray{_T,N}, a::LinearViscous, TauII::AbstractArray{_T,N}; kwargs...
) where {N,_T}
    if TauII[1] isa Quantity
        @unpack_units η = a
    else
        @unpack_val η = a
    end

    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i])
    end

    return nothing
end

function dεII_dτII(a::LinearViscous, TauII; kwargs...)
    @unpack η = a

    return 0.5*(1.0/η)
end

"""
    compute_τII(s::LinearViscous, EpsII; kwargs...)

Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor 
"""
function compute_τII(a::LinearViscous, EpsII; kwargs...)
    @unpack η = a

    return 2 * (η * EpsII)
end

function compute_τII!(
    TauII::AbstractArray{_T,N}, a::LinearViscous, EpsII::AbstractArray{_T,N}; kwargs...
) where {N,_T}
    if EpsII[1] isa Quantity
        @unpack_units η = a
    else
        @unpack_val η = a
    end

    @inbounds for i in eachindex(EpsII)
        TauII[i] = compute_τII(a, EpsII[i])
    end

    return nothing
end

function dτII_dεII(a::LinearViscous, EpsII; kwargs...)
    @unpack η = a

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

#=
# add methods programatically 
#for myType in (:LinearViscous, :DiffusionCreep, :DislocationCreep, :ConstantElasticity)
for myType in (:LinearViscous, :DiffusionCreep, :DislocationCreep)

    @eval begin
        compute_εII(TauII, a::$(myType), args) = compute_εII(TauII, a; args...) 
        dεII_dτII(TauII, a::$(myType), args) = dεII_dτII(TauII, a; args...) 
        compute_τII(EpsII, a::$(myType), args) = compute_τII(EpsII, a; args...) 
        if Symbol($myType) !== :ConstantElasticity
            dτII_dεII(EpsII, a::$(myType), args) = dτII_dεII(EpsII, a; args...)
        end
    end
end
=#

# Help info for the calculation routines
"""
    compute_εII(TauII, s:<AbstractCreepLaw, p::CreepLawVariables)

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
