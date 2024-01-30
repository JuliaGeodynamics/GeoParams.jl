# This implements viscous creep laws for melt and partially molten rocks

export LinearMeltViscosity,
        dεII_dτII,
        dτII_dεII,
        compute_εII!,
        compute_εII,
        compute_τII!,
        compute_τII

# Linear melt viscosity ------------------------------------------------
"""
    LinearMeltViscosity(η=1e20Pa*s)
    
Defines a simple temperature-dependent melt viscosity, given by 
```math  
    \\tau_{ij} = 2 \\eta \\dot{\\varepsilon}_{ij} 
```
or
```math  
    \\dot{\\varepsilon}_{ij}  = {\\tau_{ij}  \\over 2 \\eta }
```
where
```math  
    \\eta =   \\eta_0*10^( A + {B \\over (T - T0)});  
```
and `\\eta_0` is the scaling viscosity, `A` and `B` are constants, and `T_0` is a reference temperature, and `T` is the temperature [in K].

Typical parameters for basalt (default) are: `A = -9.6012`, `B = 1.3374e+04K`, `T_0 = 307.8043K` and `\\eta_0 = 1Pas`.

Typical parameters for rhyolite are: `A = -8.1590`, `B = 2.4050e+04K`, `T_0 = -430.9606K` and `\\eta_0 = 1Pas`.
"""
@with_kw_noshow struct LinearMeltViscosity{T,U,U1,U2} <: AbstractCreepLaw{T}
    A::GeoUnit{T,U}  = -9.6012NoUnits          
    B::GeoUnit{T,U1}  = 1.3374e+04K          
    T0::GeoUnit{T,U1} = 307.8043K;          # reference T
    η0::GeoUnit{T,U2} = 1Pas;               # scaling viscosity
end
LinearMeltViscosity(args...) = LinearMeltViscosity(convert.(GeoUnit, args)...)

function param_info(a::LinearMeltViscosity) # info about the struct
    return MaterialParamsInfo(; Equation=L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}; \eta=10^{A + \frac{B}{T-T_0}}")
end

# Calculation routine
@inline function compute_εII(
                    a::LinearMeltViscosity, 
                    TauII;  
                    T=one(precision(a)),
                    kwargs...,
                    )
    
    @unpack_val A,B,T0,η0 = a
    η = η0*10^(A + (B / (T - T0)))
    return (TauII / η) * 0.5
end

@inline function compute_εII(
                            a::LinearMeltViscosity, 
                            TauII::Quantity; 
                            T=1K, 
                            args...
                            )
    @unpack_units A,B,T0,η0 = a
    η = η0*10^(A + (B / (T - T0)))
    ε = (TauII / η) * 0.5

    return ε
end

"""
    compute_εII!(EpsII::AbstractArray{_T,N}, s::LinearMeltViscosity, TauII::AbstractArray{_T,N}; T, kwargs...)
"""
@inline function compute_εII!(
    EpsII::AbstractArray{_T,N}, 
    a::LinearMeltViscosity, 
    TauII::AbstractArray{_T,N}; 
    T=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...
) where {N,_T}

    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i])
    end

    return nothing
end

@inline function dεII_dτII(
    a::LinearMeltViscosity, TauII::Quantity; T=1K, kwargs...
)
    @unpack_units A,B,T0,η0 = a
    η = η0*10^(A + (B / (T - T0)))

    return  0.5 * (1.0 / η)
end

@inline function dεII_dτII(a::LinearMeltViscosity, 
    TauII;   
    T=one(precision(a)),
    kwargs...
)
    @unpack_val A,B,T0,η0 = a
    η = η0*10^(A + (B / (T - T0)))

    return 0.5 * (1.0 / η)
end

"""
    compute_τII(s::LinearMeltViscosity, EpsII; kwargs...)

Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor 
"""
@inline function compute_τII(
    a::LinearMeltViscosity, 
    EpsII; 
    T=one(precision(a)),
    kwargs...
)

    @unpack_val A,B,T0,η0 = a
    η = η0*10^(A + (B / (T - T0)))

    return 2 * (η * EpsII)
end

@inline function compute_τII(
    a::LinearMeltViscosity, EpsII::Quantity; T=1K, kwargs...
)
    @unpack_units A,B,T0,η0 = a

    η = η0*10^(A + (B / (T - T0)))

    return 2 * (η * EpsII)
end


function compute_τII!(
    TauII::AbstractArray{_T,N}, 
    a::LinearMeltViscosity, 
    EpsII::AbstractArray{_T,N}; 
    T=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...
) where {N,_T}
   
    @inbounds for i in eachindex(EpsII)
        TauII[i] = compute_τII(a, EpsII[i]; T=T[i])
    end

    return nothing
end

@inline function dτII_dεII(
    a::LinearMeltViscosity,
    EpsII;
    T=one(precision(a)),
    kwargs...,
)
    @unpack_val A,B,T0,η0 = a
    η = η0*10^(A + (B / (T - T0)))

    return 2 * η
end

@inline function dτII_dεII(
    a::LinearMeltViscosity, EpsII::Quantity; T=1K, kwargs...
)
    @unpack_units A,B,T0,η0 = a
    η = η0*10^(A + (B / (T - T0)))

    return 2 * η
end

# Print info 
function show(io::IO, g::LinearMeltViscosity)
    return print(io, "Linear melt viscosity: η=$(UnitValue(g.η0)) * 10^($(UnitValue(g.A)) + ($(UnitValue(g.B)) / (T - $(UnitValue(g.T0)))))")
end
#-------------------------------------------------------------------------
