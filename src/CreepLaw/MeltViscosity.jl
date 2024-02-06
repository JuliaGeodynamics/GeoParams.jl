# This implements viscous creep laws for melt and partially molten rocks
using SpecialFunctions

export LinearMeltViscosity, ViscosityPartialMelt_Costa_etal_2009,
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
        EpsII[i] = compute_εII(a, TauII[i], T=T[i])
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


#-------------------------------------------------------------------------
"""
    ViscosityPartialMelt_Costa_etal_2009(LinearMeltViscosity())

The viscosity of a partially molten rock depends on the melt viscosity, melt fraction and strainrate.

This implements a parameterisation of Costa et al. [2009].

Reference
===
Costa, A., Caricchi, L., Bagdassarov, N., 2009. A model for the rheology of particle‐bearing suspensions and partially molten rocks. Geochem Geophys Geosyst 10, 2008GC002138. https://doi.org/10.1029/2008GC002138

"""
@with_kw_noshow struct ViscosityPartialMelt_Costa_etal_2009{T,U, S1<:AbstractCreepLaw} <: AbstractCreepLaw{T}
    η::S1  = LinearMeltViscosity()     # Basalt melt viscosity
    ε0::GeoUnit{T,U}  = 1.0/s                           # characteristic strainrate
end
ViscosityPartialMelt_Costa_etal_2009(args...) = ViscosityPartialMelt_Costa_etal_2009(args[1]..., convert.(GeoUnit, args[2:end])...)
ViscosityPartialMelt_Costa_etal_2009(melt_viscosity::AbstractCreepLaw, args...) = ViscosityPartialMelt_Costa_etal_2009(melt_viscosity, convert.(GeoUnit, args)...)
isDimensional(g::ViscosityPartialMelt_Costa_etal_2009) = isDimensional(g.ε0)

function param_info(a::ViscosityPartialMelt_Costa_etal_2009) # info about the struct
    return MaterialParamsInfo(; Equation=L"Costa et al. (2009) effective viscosity parameterisation")
end

# viscosity correction; where c_vf=crystal volume fraction and strain_rate is the strainrate in s^-1
function viscosity_correction(c_vf, Strain_Rate)
    strain_rate = Strain_Rate
    if strain_rate < eps(Float64)
        strain_rate = 1.e-6
    end
    phi_max = 0.066499*tanh(0.913424*log10(strain_rate) + 3.850623) + 0.591806
    delta = -6.301095*tanh(0.818496*log10(strain_rate) + 2.86) + 7.462405
    alpha = -0.000378*tanh(1.148101*log10(strain_rate) + 3.92) + 0.999572
    gamma = 3.987815*tanh(0.8908*log10(strain_rate) + 3.24) + 5.099645
    num = 1.0 + (c_vf/phi_max)^delta
    x = sqrt(π)*c_vf/(2*alpha*phi_max)*(1.0 + (c_vf/phi_max)^gamma)
    den = (1.0-alpha*erf(x))^(2.5*phi_max)
    mu_r = num/den
    return mu_r
end

# Calculation routine
@inline function compute_εII(
    a::ViscosityPartialMelt_Costa_etal_2009,
    TauII;
    ϕ = one(precision(a)),
    kwargs...,
    )

    # melt viscosity
    ε      = compute_εII(a.η, TauII, kwargs)
    η_melt = TauII/(2 * ε)

    # viscosity correction factor
    ηr     = viscosity_correction(1 - ϕ, ε/a.ε0)

    η      =  ηr * η_melt

    return (TauII / η) * 0.5
end

"""
    compute_εII!(EpsII::AbstractArray{_T,N}, s::ViscosityPartialMelt_Costa_etal_2009, TauII::AbstractArray{_T,N}; T, kwargs...)
"""
@inline function compute_εII!(
    EpsII::AbstractArray{_T,N},
    a::ViscosityPartialMelt_Costa_etal_2009,
    TauII::AbstractArray{_T,N};
    ϕ=ones(size(TauII))::AbstractArray{_T,N},
    T=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...
) where {N,_T}

    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i], T=T[i])
    end

    return nothing
end


"""
    compute_τII(s::ViscosityPartialMelt_Costa_etal_2009, EpsII; kwargs...)

Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor
"""
@inline function compute_τII(
    a::ViscosityPartialMelt_Costa_etal_2009,
    EpsII;
    ϕ=one(precision(a)),
    kwargs...
)
    # melt viscosity
    τ      = compute_τII(a.η, EpsII, kwargs)
    η_melt = τ/(2 * EpsII)

    # viscosity correction factor
    ηr     = viscosity_correction(1.0 - ϕ, EpsII/a.ε0)

    η      =  ηr * η_melt
    return 2 * (η * EpsII)
end

@inline function compute_τII(
    a::ViscosityPartialMelt_Costa_etal_2009, EpsII::Quantity; ϕ=1.0, kwargs...
)

    # melt viscosity
    τ      = compute_τII(a.η, EpsII, kwargs)
    η_melt = τ/(2 * EpsII)

    # viscosity correction factor
    ηr     = viscosity_correction(1.0 - ϕ, EpsII/a.ε0)

    η      =  ηr * η_melt
    return 2 * (η * EpsII)
end


function compute_τII!(
    TauII::AbstractArray{_T,N},
    a::ViscosityPartialMelt_Costa_etal_2009,
    EpsII::AbstractArray{_T,N};
    ϕ=ones(size(TauII))::AbstractArray{_T,N},
    T=ones(size(TauII))::AbstractArray{_T,N},
    kwargs...
) where {N,_T}

    @inbounds for i in eachindex(EpsII)
        TauII[i] = compute_τII(a, EpsII[i]; ϕ=ϕ[i], T=T[i])
    end

    return nothing
end

#use AD to compute derivatives
@inline function dεII_dτII(a::ViscosityPartialMelt_Costa_etal_2009, τII; args...)
   return ForwardDiff.derivative(τII -> compute_εII(a, τII; args...), τII)
end

@inline function dτII_dεII(a::ViscosityPartialMelt_Costa_etal_2009, εII; args...)
    return ForwardDiff.derivative(εII -> compute_τII(a, εII; args...), εII)
end

# Print info
function show(io::IO, g::ViscosityPartialMelt_Costa_etal_2009)
    return print(io, "Viscosity of partially molten rocks (after Costa et al. [2009]) using melt viscosity: η=$(g.η)")
end
#-------------------------------------------------------------------------
