# This implements viscous creep laws for melt and partially molten rocks
using SpecialFunctions

export LinearMeltViscosity,
    ViscosityPartialMelt_Costa_etal_2009,
    GiordanoMeltViscosity,
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
@with_kw_noshow struct LinearMeltViscosity{T, U, U1, U2} <: AbstractCreepLaw{T}
    A::GeoUnit{T, U} = -9.6012NoUnits
    B::GeoUnit{T, U1} = 1.3374e+4K
    T0::GeoUnit{T, U1} = 307.8043K          # reference T
    η0::GeoUnit{T, U2} = 1Pas               # scaling viscosity
end
LinearMeltViscosity(args...) = LinearMeltViscosity(convert.(GeoUnit, args)...)

isDimensional(g::LinearMeltViscosity) = isDimensional(g.η0)

function param_info(a::LinearMeltViscosity) # info about the struct
    return MaterialParamsInfo(;
        Equation = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}; \eta=10^{A + \frac{B}{T-T_0}}",
    )
end

# Calculation routine
@inline function compute_εII(a::LinearMeltViscosity, TauII; T = one(precision(a)), kwargs...)
    @unpack_val A, B, T0, η0 = a
    η = η0 * 10^(A + (B / (T - T0)))
    ε = TauII / η * 0.5
    return ε
end

@inline function compute_εII(a::LinearMeltViscosity, TauII::Quantity; T = 1K, args...)
    @unpack_units A, B, T0, η0 = a
    η = η0 * 10^(A + (B / (T - T0)))
    ε = TauII / η * 0.5
    return ε
end

"""
    compute_εII!(EpsII::AbstractArray{_T,N}, s::LinearMeltViscosity, TauII::AbstractArray{_T,N}; T, kwargs...)
"""
@inline function compute_εII!(
        EpsII::AbstractArray{_T, N},
        a::LinearMeltViscosity,
        TauII::AbstractArray{_T, N};
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T = T[i])
    end

    return nothing
end

@inline function dεII_dτII(a::LinearMeltViscosity, TauII::Quantity; T = 1K, kwargs...)
    @unpack_units A, B, T0, η0 = a
    η = η0 * 10^(A + (B / (T - T0)))

    return 0.5 * (1.0 / η)
end

@inline function dεII_dτII(a::LinearMeltViscosity, TauII; T = one(precision(a)), kwargs...)
    @unpack_val A, B, T0, η0 = a
    η = η0 * 10^(A + (B / (T - T0)))

    return 0.5 * (1.0 / η)
end

"""
    compute_τII(s::LinearMeltViscosity, EpsII; kwargs...)

Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor
"""
@inline function compute_τII(a::LinearMeltViscosity, EpsII; T = one(precision(a)), kwargs...)
    @unpack_val A, B, T0, η0 = a
    η = η0 * 10^(A + (B / (T - T0)))

    return 2 * η * EpsII
end

@inline function compute_τII(a::LinearMeltViscosity, EpsII::Quantity; T = 1K, kwargs...)
    @unpack_units A, B, T0, η0 = a

    η = η0 * 10^(A + (B / (T - T0)))

    return 2 * η * EpsII
end

function compute_τII!(
        TauII::AbstractArray{_T, N},
        a::LinearMeltViscosity,
        EpsII::AbstractArray{_T, N};
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(EpsII)
        TauII[i] = compute_τII(a, EpsII[i]; T = T[i])
    end

    return nothing
end

@inline function dτII_dεII(a::LinearMeltViscosity, EpsII; T = one(precision(a)), kwargs...)
    @unpack_val A, B, T0, η0 = a
    η = η0 * 10^(A + (B / (T - T0)))

    return 2 * η
end

@inline function dτII_dεII(a::LinearMeltViscosity, EpsII::Quantity; T = 1K, kwargs...)
    @unpack_units A, B, T0, η0 = a
    η = η0 * 10^(A + (B / (T - T0)))

    return 2 * η
end

# Print info
function show(io::IO, g::LinearMeltViscosity)
    return print(
        io,
        "Linear melt viscosity: η=$(UnitValue(g.η0)) * 10^($(UnitValue(g.A)) + ($(UnitValue(g.B)) / (T - $(UnitValue(g.T0)))))",
    )
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
@with_kw_noshow struct ViscosityPartialMelt_Costa_etal_2009{T, U, S1 <: AbstractCreepLaw} <:
    AbstractCreepLaw{T}
    η::S1 = LinearMeltViscosity()     # Basalt melt viscosity
    ε0::GeoUnit{T, U} = 1.0 / s                           # characteristic strainrate
end
function ViscosityPartialMelt_Costa_etal_2009(args...)
    return ViscosityPartialMelt_Costa_etal_2009(
        args[1]..., convert.(GeoUnit, args[2:end])...
    )
end
function ViscosityPartialMelt_Costa_etal_2009(melt_viscosity::AbstractCreepLaw, args...)
    return ViscosityPartialMelt_Costa_etal_2009(melt_viscosity, convert.(GeoUnit, args)...)
end
isDimensional(g::ViscosityPartialMelt_Costa_etal_2009) = isDimensional(g.ε0)

function param_info(a::ViscosityPartialMelt_Costa_etal_2009) # info about the struct
    return MaterialParamsInfo(;
        Equation = L"Costa et al. (2009) effective viscosity parameterisation"
    )
end

# viscosity correction; where c_vf=crystal volume fraction and strain_rate is the strainrate in s^-1
function viscosity_correction(c_vf, Strain_Rate)
    strain_rate = Strain_Rate + (Strain_Rate < eps(Float64)) * 1.0e-6
    θ = log10(strain_rate)
    ϕ_max = 0.066499 * tanh(0.913424 * θ + 3.850623) + 0.591806
    δ = -6.301095 * tanh(0.818496 * θ + 2.86) + 7.462405
    α = -0.000378 * tanh(1.148101 * θ + 3.92) + 0.999572
    γ = 3.987815 * tanh(0.8908 * θ + 3.24) + 5.099645
    num = 1 + (c_vf / ϕ_max)^δ
    x = √π * c_vf / (2 * α * ϕ_max) * (1 + (c_vf / ϕ_max)^γ)
    den = (1 - α * erf(x))^(2.5 * ϕ_max)
    mu_r = num / den
    return mu_r
end

# Calculation routine
@inline function compute_εII(
        a::ViscosityPartialMelt_Costa_etal_2009, TauII; ϕ = one(precision(a)), kwargs...
    )
    # melt viscosity
    ε = compute_εII(a.η, TauII, kwargs)
    η_melt = TauII / ε * 0.5
    # viscosity correction factor
    ηr = viscosity_correction(1 - ϕ, ε / a.ε0)
    η = ηr * η_melt
    return (TauII / η) * 0.5
end

"""
    compute_εII!(EpsII::AbstractArray{_T,N}, s::ViscosityPartialMelt_Costa_etal_2009, TauII::AbstractArray{_T,N}; T, kwargs...)
"""
@inline function compute_εII!(
        EpsII::AbstractArray{_T, N},
        a::ViscosityPartialMelt_Costa_etal_2009,
        TauII::AbstractArray{_T, N};
        ϕ = ones(size(TauII))::AbstractArray{_T, N},
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T = T[i])
    end

    return nothing
end

"""
    compute_τII(s::ViscosityPartialMelt_Costa_etal_2009, EpsII; kwargs...)

Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor
"""
@inline function compute_τII(
        a::ViscosityPartialMelt_Costa_etal_2009, EpsII; ϕ = one(precision(a)), kwargs...
    )
    # melt viscosity
    τ = compute_τII(a.η, EpsII, kwargs)
    η_melt = τ / EpsII * 0.5
    # viscosity correction factor
    ηr = viscosity_correction(1.0 - ϕ, EpsII / a.ε0)
    η = ηr * η_melt
    return 2 * η * EpsII
end

@inline function compute_τII(
        a::ViscosityPartialMelt_Costa_etal_2009, EpsII::Quantity; ϕ = 1.0, kwargs...
    )

    # melt viscosity
    τ = compute_τII(a.η, EpsII, kwargs)
    η_melt = τ / EpsII * 0.5

    # viscosity correction factor
    ηr = viscosity_correction(1.0 - ϕ, EpsII / a.ε0)

    η = ηr * η_melt
    return 2 * η * EpsII
end

function compute_τII!(
        TauII::AbstractArray{_T, N},
        a::ViscosityPartialMelt_Costa_etal_2009,
        EpsII::AbstractArray{_T, N};
        ϕ = ones(size(TauII))::AbstractArray{_T, N},
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(EpsII)
        TauII[i] = compute_τII(a, EpsII[i]; ϕ = ϕ[i], T = T[i])
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
    return print(
        io,
        "Viscosity of partially molten rocks (after Costa et al. [2009]) using melt viscosity: η=$(g.η)",
    )
end
#-------------------------------------------------------------------------
"""
    GiordanoMeltViscosity(; oxd_wt = oxd_wt, η0=1Pas)

Defines the melt viscosity model after Giordano et al. (2008) given by
```math
    \\tau_{ij} = 2 \\eta \\dot{\\varepsilon}_{ij}
```
or
```math
    \\dot{\\varepsilon}_{ij}  = {\\tau_{ij}  \\over 2 \\eta }
```
where
```math
    \\eta = \\eta_0 * 10^{(AT + BT \\over (T - CT))}
```

## Parameters
- `oxd_wt::NTuple{9,T}`: Melt composition as 9-element Tuple containing concentrations
             in [wt%] of the following oxides ordered in the exact sequence \n
             (SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O) \n
             Default values are for a hydrous N-MORB melt.

## Reference
- Giordano D, Russell JK, & Dingwell DB (2008). Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science Letters, 271, 123-134. (https://dx.doi.org/10.1016/j.epsl.2008.03.038)

"""
struct GiordanoMeltViscosity{T, T1, T2, U, U1, U2} <: AbstractCreepLaw{T}
    oxd_wt::NTuple{9, T} # SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
    MW::NTuple{9, T} # Molar weights
    bb::NTuple{10, T1} # model constants
    cc::NTuple{7, T2}  # model constants
    AT::GeoUnit{T, U}
    BT::GeoUnit{T, U1}
    CT::GeoUnit{T, U1}
    η0::GeoUnit{T, U2} # scaling viscosity

    function GiordanoMeltViscosity(;
        oxd_wt = (50.42, 1.53,  15.13,  9.81, 7.76, 11.35, 2.83,  0.14, 1.0),
        MW = (60.0843, 79.8658, 101.961276, 71.8444, 40.3044,56.0774, 61.97894, 94.1960, 18.01528), # Molar weights
        bb = (159.56, -173.34, 72.13, 75.69, -38.98, -84.08, 141.54, -2.43, -0.91, 17.62).*K,
        cc = (2.75, 15.72, 8.32, 10.20, -12.29, -99.54, 0.3).*K,
        AT = -4.55NoUnits,
        BT = 5515.26K,
        CT = 489.21K,
        η0 = 1Pas
        )

        BT, CT = calculate_BT_CT(oxd_wt, MW, bb, cc)

        ATU = convert(GeoUnit, AT)
        BTU = convert(GeoUnit, BT)
        CTU = convert(GeoUnit, CT)
        η0U = convert(GeoUnit, η0)

        bbU = ntuple(i -> convert(GeoUnit, bb[i]), 10)
        ccU = ntuple(i -> convert(GeoUnit, cc[i]), 7)

        T  = eltype(oxd_wt)
        T1 = eltype(bbU)
        T2 = eltype(ccU)
        U  = typeof(ATU).types[2]
        U1 = typeof(BTU).types[2]
        U2 = typeof(η0U).types[2]

        new{T, T1, T2, U, U1, U2}(oxd_wt, MW, bbU, ccU, ATU, BTU, CTU, η0U)
    end

    function GiordanoMeltViscosity(oxd_wt, MW, bb, cc, AT, BT, CT, η0)
        return GiordanoMeltViscosity(;
            oxd_wt = oxd_wt, MW = MW, bb = bb, cc = cc, AT = AT, BT = BT, CT = CT, η0 = η0
        )
    end

end

isDimensional(g::GiordanoMeltViscosity) = isDimensional(g.η0)

function param_info(a::GiordanoMeltViscosity) # info about the struct
    return MaterialParamsInfo(; Equation = L"\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}; \eta=10^{(AT + BT \over (T - CT))}")
end

# function calculate_BT_CT!(BT, CT, oxd_wt, MW, bb, cc)
function calculate_BT_CT(oxd_wt, MW, bb, cc)
    tmp = ntuple(i -> oxd_wt[i], Val(8))
    sum_oxd_wt  = sum(tmp)
    α = (100.0 - oxd_wt[9]) / sum_oxd_wt
    oxd_wt_norm = (
        @.(tmp * α)...,
        oxd_wt[9],
    )
    oxd_mol = oxd_wt_norm ./ MW
    oxd_mol = oxd_mol ./ sum(oxd_mol) .* 100
    # Load composition-basis matrix for multiplication against model-coefficients
    siti = oxd_mol[1] + oxd_mol[2]
    tial = oxd_mol[2] + oxd_mol[3]
    fmm  = oxd_mol[4] + oxd_mol[5]
    nak  = oxd_mol[7] + oxd_mol[8]
    b1   = siti
    b2   = oxd_mol[3]
    b3   = oxd_mol[4]
    b4   = oxd_mol[5]
    b5   = oxd_mol[6]
    b6   = oxd_mol[7] + oxd_mol[9]
    b7   = oxd_mol[9] + log(1.0 + oxd_mol[9])
    b12  =  siti * fmm
    b13  = (siti + oxd_mol[3]) * (nak + oxd_mol[9])
    b14  = oxd_mol[3] * nak

    c1   = oxd_mol[1]
    c2   = tial
    c3   = fmm
    c4   = oxd_mol[6]
    c5   = nak
    c6   = log(1 + oxd_mol[9])
    c11  = (oxd_mol[3] + fmm + oxd_mol[6]) * (nak + oxd_mol[9])
    bcf  = (b1, b2, b3, b4, b5, b6, b7, b12, b13, b14)
    ccf  = (c1, c2, c3, c4, c5, c6, c11)

    BT = sum(bb .* bcf)
    CT = sum(cc .* ccf)
    return BT, CT
end


#calculation routine
function compute_εII(a::GiordanoMeltViscosity, TauII; T = one(precision(a)), kwargs...)
    @unpack_val AT, BT, CT, η0 = a

    η =  η0 * exp10(min(12, max(-6, AT + BT /(T - CT))))

    ε = (TauII / η) * 0.5
    return ε
end

function compute_εII(a::GiordanoMeltViscosity, TauII::Quantity; T = 1K, kwargs...)
    @unpack_units AT, BT, CT, η0 = a

    η =  η0 * exp10(min(12, max(-6, AT + BT /(T - CT))))
    ε = TauII / η * 0.5
    return ε
end

"""
    compute_εII!(EpsII::AbstractArray{_T,N}, s::GiordanoMeltViscosity, TauII::AbstractArray{_T,N}; T, kwargs...)
"""
@inline function compute_εII!(
        EpsII::AbstractArray{_T, N},
        a::GiordanoMeltViscosity,
        TauII::AbstractArray{_T, N};
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T = T[i])
    end
    return nothing
end

@inline function dεII_dτII(a::GiordanoMeltViscosity, TauII::Quantity; T = 1K, kwargs...)
    @unpack_units AT, BT, CT,η0 = a

    η =  η0 * exp10(min(12, max(-6, AT + BT /(T - CT))))

    return 0.5 * (1.0 / η)
end

@inline function dεII_dτII(a::GiordanoMeltViscosity, TauII; T = one(precision(a)), kwargs...)
    @unpack_val AT, BT, CT, η0 = a

    η =  η0 * exp10(min(12, max(-6, AT + BT /(T - CT))))

    return 0.5 * (1.0 / η)
end

"""
    compute_τII(s::GiordanoMeltViscosity, EpsII; kwargs...)

Returns second invariant of the stress tensor given a 2nd invariant of strain rate tensor
"""
@inline function compute_τII(a::GiordanoMeltViscosity, EpsII; T = one(precision(a)), kwargs...)
    @unpack_val AT, BT, CT, η0 = a

    η =  η0 * exp10(min(12, max(-6, AT + BT /(T - CT))))

    return 2 * η * EpsII
end

@inline function compute_τII(a::GiordanoMeltViscosity, EpsII::Quantity; T = 1K, kwargs...)
    @unpack_units AT, BT, CT,η0 = a

    η =  η0 * exp10(min(12, max(-6, AT + BT /(T - CT))))

    return 2 * η * EpsII
end

@inline function compute_τII!(
        TauII::AbstractArray{_T, N},
        a::GiordanoMeltViscosity,
        EpsII::AbstractArray{_T, N};
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(TauII)
        TauII[i] = compute_τII(a, EpsII[i]; T = T[i])
    end
    return nothing
end

@inline function dτII_dεII(a::GiordanoMeltViscosity, EpsII; T = one(precision(a)), kwargs...)
    @unpack_val AT, BT, CT, η0 = a

    η =  η0 * exp10(min(12, max(-6, AT + BT /(T - CT))))

    return 2 * η
end

@inline function dτII_dεII(a::GiordanoMeltViscosity, EpsII::Quantity; T = 1K, kwargs...)
    @unpack_units AT, BT, CT,η0 = a

    η =  η0 * exp10(min(12, max(-6, AT + BT /(T - CT))))

    return 2 * η
end

# Print info
function show(io::IO, g::GiordanoMeltViscosity)
    return print(
        io,
        "GiordanoMeltViscosity: η= $(UnitValue(g.η0)) * 10^(AT + BT / (T - CT)) with Oxide Comp $(g.oxd_wt)",
    )
end
#-------------------------------------------------------------------------
