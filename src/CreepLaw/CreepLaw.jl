# This implements viscous creep laws and routines to compute with them
#
# Note that various simple creep laws are defined in this file;
# more complex ones (such as DislocationCreep) are in separate files
# in the same directory
#
# In case you want to add new creep laws, have a look at how the ones
# here are implemented. Please add tests as well!

# include("Data/DiffusionCreep.jl")
# include("Data/DislocationCreep.jl")
# include("Data/GrainBoundarySliding.jl")
# include("Data/NonLinearPeierlsCreep.jl")
# include("Data/PeierlsCreep.jl")

abstract type AbstractCreepLaw{T} <: AbstractConstitutiveLaw{T} end

export isvolumetric,
    LinearViscous, PowerlawViscous, CorrectionFactor, AbstractCreepLaw, ArrheniusType

isvolumetric(a::AbstractCreepLaw) = false

# This computes correction factors to go from experimental data to tensor format
function CorrectionFactor(a::AbstractCreepLaw{_T}) where {_T}
    apparatus = a.Apparatus
    if apparatus == AxialCompression
        FT = sqrt(3) # relation between differential stress recorded by apparatus and TauII
        FE = 2 / FT            # relation between gamma recorded by apparatus and EpsII
        return FT, FE

    elseif apparatus == SimpleShear
        FT = 2.0      # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = FT
        return FT, FE

    elseif apparatus == Invariant
        FT = FE = 1.0
        return FT, FE

    else
        throw(ArgumentError("Unknown apparatus"))

    end
end

function CorrectionFactor(a)
    if a == AxialCompression
        FT = sqrt(3)           # relation between differential stress recorded by apparatus and TauII
        FE = 2 / FT            # relation between gamma recorded by apparatus and EpsII
        return FT, FE

    elseif a == SimpleShear
        FT = 2.0      # it is assumed that the flow law parameters were derived as a function of differential stress, not the shear stress. Must be modidified if it is not the case
        FE = FT
        return FT, FE

    elseif a == Invariant
        FT = FE = 1.0
        return FT, FE

    else
        throw(ArgumentError("Unknown apparatus"))

    end
end

@inline rat2float(x::Rational) = float(x)
@inline rat2float(x) = x

@inline unit_power(A) = typeof(A).parameters[2].parameters[1][1].power
include("DislocationCreep.jl")
include("DiffusionCreep.jl")
include("GrainBoundarySliding.jl")
include("PeierlsCreep.jl")
include("NonLinearPeierlsCreep.jl")
include("CustomRheology.jl")
include("MeltViscosity.jl")

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
@with_kw_noshow struct LinearViscous{T, U} <: AbstractCreepLaw{T}
    η::GeoUnit{T, U} = 1.0e20Pa * s                # viscosity
    η_val::T = 1.0
end
LinearViscous(args...) = LinearViscous(convert(GeoUnit, args[1]), args[2])

function param_info(a::LinearViscous) # info about the struct
    return MaterialParamsInfo(; Equation = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}")
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
        EpsII::AbstractArray{_T, N}, a::LinearViscous, TauII::AbstractArray{_T, N}; kwargs...
    ) where {N, _T}
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

    return 0.5 * (1.0 / η)
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
        TauII::AbstractArray{_T, N}, a::LinearViscous, EpsII::AbstractArray{_T, N}; kwargs...
    ) where {N, _T}
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

# ArrheniusType temperature dependent viscosity --------------------------
"""
    ArrheniusType()
    
Defines an Arrhenius-type linear viscosity in the form of 
    \\eta = \\eta_{0} * \\exp\\left(\\frac{E_{\\eta}{T+T_{0}}-\\frac{E_{\\eta}{T_{\\eta}+T_{0}}\\right) 

The (isotropic) linear viscous rheology is given by  
```math  
    \\tau_{ij} = 2 \\eta \\dot{\\varepsilon}_{ij} 
```
or
```math  
    \\dot{\\varepsilon}_{ij}  = {\\tau_{ij}  \\over 2 \\eta }
```

where ``\\eta_0`` is the reference viscosity [Pa*s] at reference strain rate ``\\dot{\\varepsilon}_0``[1/s], and ``n`` the power law exponent []. 
"""
@with_kw_noshow struct ArrheniusType{_T, U1, U2, U3} <: AbstractCreepLaw{_T}
    η0::GeoUnit{_T, U1} = 1.0NoUnits      # Pre-exponential factor
    E_η::GeoUnit{_T, U2} = 23.03NoUnits    # Activation energy (non-dimensional)
    T_O::GeoUnit{_T, U3} = 1.0NoUnits      # Offset temperature
    T_η::GeoUnit{_T, U3} = 1.0NoUnits      # Reference temperature at which viscosity is unity
end
#ArrheniusType(args...) = ArrheniusType(args[1], args[2], args[3], args[4])

function ArrheniusType(args...)
    return ArrheniusType(
        convert(GeoUnit, args[1]),
        convert(GeoUnit, args[2]),
        convert(GeoUnit, args[3]),
        convert(GeoUnit, args[4]),
    )
end

function param_info(a::ArrheniusType) # info about the struct
    return MaterialParamsInfo(;
        Equation = L"\tau_{ij} = 2 \eta_0 exp( E_η/(T + T_O) + E_η/(T_η + T_O))  \dot{\varepsilon}_{ij}",
    )
end

# Calculation routines for linear viscous rheologies
function compute_εII(a::ArrheniusType, TauII; T = one(precision(a)), kwargs...)
    @unpack_val η0, E_η, T_O, T_η = a
    η = η0 * exp(E_η / (T + T_O) - E_η / (T_η + T_O))
    return (TauII / η) * 0.5
end

"""
    
    compute_εII!(EpsII::AbstractArray{_T,N}, s::ArrheniusType, TauII::AbstractArray{_T,N})
"""
function compute_εII!(
        EpsII::AbstractArray{_T, N},
        a::ArrheniusType,
        TauII::AbstractArray{_T, N};
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T = T[i])
    end

    return nothing
end

function dεII_dτII(a::ArrheniusType, TauII; T = one(precision(a)), kwargs...)
    @unpack_val η0, E_η, T_O, T_η = a
    η = η0 * exp(E_η / (T + T_O) - E_η / (T_η + T_O))
    return 0.5 * inv(η)
end

"""
    compute_τII(s::ArrheniusType, EpsII; kwargs...)

Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor 
"""
function compute_τII(a::ArrheniusType, EpsII; T = one(precision(a)), kwargs...)
    @unpack_val η0, E_η, T_O, T_η = a

    η = η0 * exp(E_η / (T + T_O) - E_η / (T_η + T_O))

    return 2 * (η * EpsII)
end

function compute_τII!(
        TauII::AbstractArray{_T, N},
        a::ArrheniusType,
        EpsII::AbstractArray{_T, N};
        T = ones(size(EpsII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(EpsII)
        TauII[i] = compute_τII(a, EpsII[i]; T = T[i])
    end

    return nothing
end

function dτII_dεII(a::ArrheniusType, EpsII; T = one(precision(a)), kwargs...)
    @unpack_val η0, E_η, T_O, T_η = a
    η = η0 * exp(E_η / (T + T_O) - E_η / (T_η + T_O))

    return 2 * η
end

# Print info
function show(io::IO, g::ArrheniusType)
    return print(
        io,
        "ArrheniusType: η0 = $(Value(g.η0)), E_η = $(Value(g.E_η)), T_O = $(Value(g.T_O)), T_η = $(Value(g.T_η))",
    )
end

# Powerlaw viscous rheology ----------------------------------------------
"""
    PowerlawViscous(η0=1e18Pa*s, n=2.0NoUnits, ε0=1e-15/s)
    
Defines a power law viscous creeplaw as: 

```math  
    \\tau_{ij}^n  = 2 \\eta_0 \\left( \\dot{\\varepsilon}_{ij} \\over \\dot{\\varepsilon}_0 \\right)
```
where ``\\eta`` is the effective viscosity [Pa*s].
"""
@with_kw_noshow struct PowerlawViscous{T, U1, U2, U3} <: AbstractCreepLaw{T}
    η0::GeoUnit{T, U1} = 1.0e18Pa * s       # reference viscosity
    n::GeoUnit{T, U2} = 2.0 * NoUnits     # powerlaw exponent
    ε0::GeoUnit{T, U3} = 1.0e-15 / s          # reference strainrate
end
PowerlawViscous(a...) = PowerlawViscous(convert.(GeoUnit, a)...)

# Print info
function show(io::IO, g::PowerlawViscous)
    return print(io, "Powerlaw viscosity: η0=$(g.η0.val), n=$(g.n.val), ε0=$(g.ε0.val) ")
end

# Calculation routines for linear viscous rheologies
function compute_εII(a::PowerlawViscous, TauII; kwargs...)
    @unpack_val η0, n, ε0 = a

    @pow EpsII = (TauII / η0)^(1 / n) * ε0

    return EpsII
end

function dεII_dτII(a::PowerlawViscous, TauII; kwargs...)
    @unpack_val η0, n, ε0 = a

    return @pow (TauII^((1 - n) / n)) / (n * (η0^(1 / n)))
end

"""
    compute_τII(s::PowerlawViscous, EpsII; kwargs...)

Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor 
"""
function compute_τII(a::PowerlawViscous, EpsII; kwargs...)
    @unpack_val η0, n, ε0 = a

    τII = @pow ε0 * η0 * (EpsII / ε0)^n

    return τII
end

function dτII_dεII(a::PowerlawViscous, EpsII; kwargs...)
    @unpack_val η0, n, ε0 = a

    return @pow η0 * EpsII^(n - 1) * (1 / ε0)^n
end
#-------------------------------------------------------------------------

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
