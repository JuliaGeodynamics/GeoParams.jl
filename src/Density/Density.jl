module Density

# This implements different methods to compute density
#
# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, Unitful, LaTeXStrings, MuladdMacro
using ..Units
using ..PhaseDiagrams
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, @extractors, add_extractor_functions, AbstractPhaseDiagramsStruct
using GeoParams: fastpow, pow_check, @pow
import ..Units: isdimensional
using ..MaterialParameters: No_MaterialParam, MaterialParamsInfo
import Base.show, GeoParams.param_info
using GeoParams: LinearInterpolator, interpolate, interpolate_field

include("../Computations.jl")

abstract type AbstractDensity{T} <: AbstractMaterialParam end
abstract type ConduitDensity{T} <: AbstractDensity{T} end

export compute_density, # calculation routines
    compute_density!, # in place calculation
    compute_density_ratio,
    param_info, # info about the parameters
    AbstractDensity,
    ConduitDensity,
    ConstantDensity, # constant
    PT_Density, # P & T dependent density
    Compressible_Density, # Compressible density
    T_Density, # T dependent density
    Vector_Density, # Vector with density
    MeltDependent_Density, # Melt dependent density
    BubbleFlow_Density, # Bubble flow density
    GasPyroclast_Density, # Gas-Pyroclast mixture density
    RedlichKwong_Density, # Modified Redlich-Kwong gas EOS (Huber et al. 2010)
    IdealGas_Density, # Ideal-gas EOS
    ThreePhase_Density, # melt + crystal + gas mixture density
    Melt_DensityX,
    get_α

# Define "empty" computational routines in case nothing is defined
function compute_density!(
        rho, s::No_MaterialParam{_T}; P::_T1 = 0.0e0, T::_T2 = 0.0e0
    ) where {_T, _T1, _T2}
    return zero(_T)
end
function compute_density(s::No_MaterialParam{_T}; P::_T1 = 0.0e0, T::_T2 = 0.0e0) where {_T, _T1, _T2}
    return zero(_T)
end

# Constant Density -------------------------------------------------------
"""
    ConstantDensity(ρ=2900kg/m^3)

Set a constant density:
```math
    \\rho  = cst
```
where ``\\rho`` is the density [``kg/m^3``].
"""
@with_kw_noshow struct ConstantDensity{_T, U} <: AbstractDensity{_T}
    ρ::GeoUnit{_T, U} = 2900.0kg / m^3 # density
end
ConstantDensity(args...) = ConstantDensity(convert.(GeoUnit, args)...)
isdimensional(s::ConstantDensity) = isdimensional(s.ρ)

@inline (ρ::ConstantDensity)(; P = 0.0e0, T = 0.0e0, args...) =
    (P isa Quantity || T isa Quantity) ? UnitValue(ρ.ρ) : ρ.ρ.val
@inline (ρ::ConstantDensity)(args) = ρ(; args...)
@inline compute_density(s::ConstantDensity{_T}, args) where {_T} = s(; args...)
@inline compute_density(s::ConstantDensity{_T}) where {_T} = s()

# This assumes that density always has a single parameter. If that is not the case, we will have to extend this (to be done)
function param_info(s::ConstantDensity) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho = cst")
end

# Calculation routines
function compute_density!(rho::AbstractArray, s::ConstantDensity; kwargs...)
    @unpack_val ρ = s
    rho[:] .= ρ
    return nothing
end

function compute_density!(rho::AbstractArray, s::ConstantDensity, args)
    return compute_density!(rho, s; args...)
end

# Print info
function show(io::IO, g::ConstantDensity)
    return print(io, "Constant density: ρ=$(UnitValue(g.ρ))")
end
#-------------------------------------------------------------------------

# Pressure & Temperature dependent density -------------------------------
"""
    PT_Density(ρ0=2900kg/m^3, α=3e-5/K, β=1e-9/Pa, T0=0C, P=0MPa)

Set a pressure and temperature-dependent density:
```math
    \\rho  = \\rho_0 (1.0 - \\alpha (T-T_0) + \\beta  (P-P_0) )
```
where ``\\rho_0`` is the density [``kg/m^3``] at reference temperature ``T_0`` and pressure ``P_0``,
``\\alpha`` is the temperature dependence of density and ``\\beta`` the pressure dependence.
"""
@with_kw_noshow struct PT_Density{_T, U1, U2, U3, U4, U5} <: AbstractDensity{_T}
    ρ0::GeoUnit{_T, U1} = 2900.0kg / m^3 # density
    α::GeoUnit{_T, U2} = 3.0e-5 / K       # T-dependence of density
    β::GeoUnit{_T, U3} = 1.0e-9 / Pa      # P-dependence of density
    T0::GeoUnit{_T, U4} = 0.0C           # Reference temperature
    P0::GeoUnit{_T, U5} = 0.0MPa         # Reference pressure
end
PT_Density(args...) = PT_Density(convert.(GeoUnit, args)...)
isdimensional(s::PT_Density) = isdimensional(s.ρ0)

function param_info(s::PT_Density) # info
    return MaterialParamsInfo(;
        Equation = L"\rho = \rho_0(1.0-\alpha (T-T_0) + \beta (P-P_0)"
    )
end

# Calculation routine
@inline function (ρ::PT_Density)(; P::Number = 0.0e0, T::Number = 0.0e0, kwargs...)
    if T isa Quantity
        @unpack_units ρ0, α, β, P0, T0 = ρ
    else
        @unpack_val   ρ0, α, β, P0, T0 = ρ
    end

    return @muladd ρ0 * (1.0 - α * (T - T0) + β * (P - P0))
end

@inline (ρ::PT_Density)(args) = ρ(; args...)
@inline compute_density(s::PT_Density, args) = s(args)

# Print info
function show(io::IO, g::PT_Density)
    return print(
        io,
        "P/T-dependent density: ρ0=$(UnitValue(g.ρ0)), α=$(UnitValue(g.α)), β=$(UnitValue(g.β)), T0=$(UnitValue(g.T0)), P0=$(UnitValue(g.P0))",
    )
end
#-------------------------------------------------------------------------

# Pressure-dependent density -------------------------------
"""
    Compressible_Density(ρ0=2900kg/m^3, β=1e-9/Pa, P₀=0MPa)

Set a pressure-dependent density:
```math
    \\rho  = \\rho_0 \\exp(β*(P - P\\_0))
```
where ``\\rho_0`` is the density [``kg/m^3``] at reference pressure ``P_0`` and ``\\beta`` the pressure dependence.
"""
@with_kw_noshow struct Compressible_Density{_T, U1, U2, U3} <: AbstractDensity{_T}
    ρ0::GeoUnit{_T, U1} = 2900.0kg / m^3 # density
    β::GeoUnit{_T, U2} = 1.0e-9 / Pa     # P-dependence of density
    P0::GeoUnit{_T, U3} = 0.0MPa         # Reference pressure
end
Compressible_Density(args...) = Compressible_Density(convert.(GeoUnit, args)...)
isdimensional(s::Compressible_Density) = isdimensional(s.ρ0)

function param_info(s::Compressible_Density) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho = \rho_0\exp(\beta*(P-P_0))")
end

function (s::Compressible_Density{_T})(; P = 0.0e0, kwargs...) where {_T}
    if P isa Quantity
        @unpack_units ρ0, β, P0 = s
    else
        @unpack_val   ρ0, β, P0 = s
    end

    return ρ0 * exp(β * (P - P0))
end

@inline (s::Compressible_Density)(args) = s(; args...)
@inline compute_density(s::Compressible_Density, args) = s(args)

# Print info
function show(io::IO, g::Compressible_Density)
    return print(
        io,
        "Compressible density: ρ0=$(UnitValue(g.ρ0)), β=$(UnitValue(g.β)), P0=$(UnitValue(g.P0))",
    )
end
#-------------------------------------------------------------------------

# Temperature-dependent density -------------------------------
"""
    T_Density(ρ0=2900kg/m^3, α=3e-5/K, T₀=273.15K)

Set a temperature-dependent density:
```math
    \\rho  = \\rho_0 (1 - \\alpha * (T - T\\_0) )
```
where ``\\rho_0`` is the density [``kg/m^3``] at reference temperature ``T_0`` and ``\\alpha`` the temperature dependence.
"""
@with_kw_noshow struct T_Density{_T, U1, U2, U3} <: AbstractDensity{_T}
    ρ0::GeoUnit{_T, U1} = 2900.0kg / m^3 # density
    α::GeoUnit{_T, U2} = 3.0e-5 / K       # T-dependence of density
    T0::GeoUnit{_T, U3} = 273.15K        # Reference temperature
end
T_Density(args...) = T_Density(convert.(GeoUnit, args)...)
isdimensional(s::T_Density) = isdimensional(s.ρ0)

function param_info(s::T_Density) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho = \rho_0*(1 - \alpha*(T-T_0))")
end

function (s::T_Density{_T})(; T = 0.0e0, kwargs...) where {_T}
    if T isa Quantity
        @unpack_units ρ0, α, T0 = s
    else
        @unpack_val   ρ0, α, T0 = s
    end

    return @muladd ρ0 * (1 - α * (T - T0))
end

@inline (s::T_Density)(args) = s(; args...)
@inline compute_density(s::T_Density, args) = s(args)

# Print info
function show(io::IO, g::T_Density)
    return print(
        io,
        "Temperature dependent density:  ρ = $(UnitValue(g.ρ0))(1 - $(UnitValue(g.α))(T-$(UnitValue(g.T0))))",
    )
end
#-------------------------------------------------------------------------

# Melt-dependent density -------------------------------------------------
"""
    MeltDependent_Density(ρsolid=ConstantDensity(), ρmelt=ConstantDensity())

If we use a single phase code the average density of a partially molten rock is
```math
    \\rho  = \\phi \\rho_{\\textrm{melt}} + (1-\\phi) \\rho_{\\textrm{solid}}
```
where ``\\rho`` is the average density [``kg/m^3``], ``\\rho_{\\textrm{melt}}`` the melt density, ``\\rho_{\\textrm{solid}} `` the solid density and ``\\phi`` the melt fraction.

Note that any density formulation can be used for melt and solid.
"""
@with_kw_noshow struct MeltDependent_Density{_T, U, S1 <: AbstractDensity, S2 <: AbstractDensity} <: AbstractDensity{_T}
    ρsolid::S1 = ConstantDensity(ρ = 2900kg / m^3) # density of the solid
    ρmelt::S2 = ConstantDensity(ρ = 2200kg / m^3)  # density of the melt
    ρ::GeoUnit{_T, U} = 2900.0kg / m^3          # to keep track on whether this struct is dimensional or not
end

MeltDependent_Density(args...) = MeltDependent_Density(args[1], args[2], convert.(GeoUnit, args[3:end])...)
isdimensional(s::MeltDependent_Density) = isdimensional(s.ρsolid)

# This assumes that density always has a single parameter. If that is not the case, we will have to extend this (to be done)
function param_info(s::MeltDependent_Density) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho =  \phi \rho_{\textrm{melt}} + (1-\phi) \\rho_{\textrm{solid}}")
end

# Calculation routines
function (rho::MeltDependent_Density{_T})(; ϕ = 0.0e0, kwargs...) where {_T}
    ρsolid = compute_density(rho.ρsolid, kwargs)
    ρmelt = compute_density(rho.ρmelt, kwargs)

    return @muladd ϕ * ρmelt + (1 - ϕ) * ρsolid
end

@inline (s::MeltDependent_Density)(args) = s(; args...)
@inline compute_density(s::MeltDependent_Density, args) = s(args)

# Print info
function show(io::IO, g::MeltDependent_Density)
    return print(io, "Melt dependent density: ρ = (1-ϕ)*ρsolid + ϕ*ρmelt; ρsolid=$(g.ρsolid); ρmelt=$(g.ρmelt)")
end
#-------------------------------------------------------------------------

# Conduit densities -------------------------------------------------
"""
    BubbleFlow_Density(ρmelt=ConstantDensity(), ρgas=ConstantDensity(), c0=0e0, a=0.0041MPa^-1/2)

Defines the BubbleFlow_Density as described in Slezin (2003) with a default gas solubility constant of 0.0041MPa``^{-1/2}`` used in e.g. Sparks et al. (1978)
```math
    \\rho = \\frac{1}{\\frac{c_0 - c}{\\rho_g} + \\frac{1-(c_0-c)}{\\rho_m}}
```
with
```math
c =
\\begin{cases}
   aP^{1/2} & \\text{for } P < \\frac{c_0^2}{a^2} \\\\
    c_0 & \\text{for } P \\geq \\frac{c_0^2}{a^2}
\\end{cases}
```
# Arguments
- `ρmelt`: Density of the melt
- `ρgas`: Density of the gas
- `c0`: Total volatile content
- `a`: Gas solubility constant (default: 4.1e-6Pa``^{-1/2}``) (after Sparks et al., 1978)

Possible values for a are 3.2e-6-6.4e-6Pa``^{-1/2}`` where the lower value corresponds to mafic magmas at rather large pressures (400-600MPa) and the higher value to felsic magmas at low pressures (0 to 100-200MPa) (after Slezin (2003))

# Example
```julia
rheology = SetMaterialParams(;
                      Phase=1,
                      CreepLaws=(PowerlawViscous(), LinearViscous(; η=1e21Pa * s)),
                      Gravity=ConstantGravity(; g=9.81.0m / s^2),
                      Density= BubbleFlow_Density(ρmelt=ConstantDensity(ρ=2900kg/m^3), ρgas=ConstantDensity(ρ=1kg/m^3), c0=0.0, a=0.0041MPa^-1//2),
                      )
```

# References
- Slezin, Yu. B. (2003), The mechanism of volcanic eruptions (a steady state approach), Journal of Volcanology and Geothermal Research, 122, 7-50, https://doi.org/10.1016/S0377-0273(02)00464-X
- Sparks, R. S. J.(1978), The dynamics of bubble formation and growth in magmas: A review and analysis, Journal of Volcanology and Geothermal Research, 3, 1-37, https://doi.org/10.1016/0377-0273(78)90002-1
"""
@with_kw_noshow struct BubbleFlow_Density{_T, U1, U2, U3, S1 <: AbstractDensity, S2 <: AbstractDensity} <: ConduitDensity{_T}
    ρmelt::S1 = ConstantDensity(ρ = 2200kg / m^3)   # density of the melt
    ρgas::S2 = ConstantDensity(ρ = 1kg / m^3)       # density of the gas
    c0::GeoUnit{_T, U1} = 0.0e0 * NoUnits         # total volatile content
    a::GeoUnit{_T, U2} = 4.1e-6Pa^(-1 // 2)         # gas solubility constant
    ρ::GeoUnit{_T, U3} = 2900.0kg / m^3          # to keep track on whether this struct is dimensional or not
end

BubbleFlow_Density(args...) = BubbleFlow_Density(args[1], args[2], convert.(GeoUnit, args[3:end])...)
isdimensional(s::BubbleFlow_Density) = isdimensional(s.ρmelt)

function param_info(s::BubbleFlow_Density) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho = 1/((c_0-c)/rho_g + 1-(c_0-c)/\rho_m)")
end

# Calculation routines
@inline function (rho::BubbleFlow_Density{_T})(; P = 0.0e0, kwargs...) where {_T}
    ρmelt = compute_density(rho.ρmelt, kwargs)
    ρgas = compute_density(rho.ρgas, kwargs)
    if P isa Quantity
        @unpack_units c0, a = rho
    else
        @unpack_val c0, a = rho
    end

    cutoff = c0^2 / a^2

    if P < cutoff
        c = a * sqrt(abs(P))
    else
        c = c0
    end

    return inv((c0 - c) / ρgas + (1 - (c0 - c)) / ρmelt)
end

@inline (s::BubbleFlow_Density)(args) = s(; args...)
@inline compute_density(s::BubbleFlow_Density, args) = s(args)

# Print info
function show(io::IO, g::BubbleFlow_Density)
    return print(io, "Bubble flow density: ρ = 1/((c0-c)/ρgas + (1-(c0-c))/ρmelt); ρmelt=$(g.ρmelt); ρgas=$(g.ρgas); c0=$(UnitValue(g.c0)); a=$(UnitValue(g.a))")
end

# Gas-Pyroclast mixture density
"""
    GasPyroclast_Density(ρmelt=ConstantDensity(), ρgas=ConstantDensity(), δ=0e0)

Defines the GasPyroclast_Density as described in Slezin (2003) with a default volume fraction of free gas in the flow of 0.0
This is also used to model partly destroyed foam in the conduit.

```math
    \\rho = \\rho_g\\delta + \\rho_p(1 - \\delta)
```
with
```math
    \\rho_p = \\rho_m(1 - \\beta) + \\rho_g\\beta \\approx \\rho_l(1 - \\beta)
```

# Arguments
- `ρmelt`: Density of the melt
- `ρgas`: Density of the gas
- `δ`: Volume fraction of free gas in the flow
- `β`: Gas volume fraction enclosed within the particles

# Example
```julia
rheology = SetMaterialParams(;
                      Phase=1,
                      CreepLaws=(PowerlawViscous(), LinearViscous(; η=1e21Pa * s)),
                      Gravity=ConstantGravity(; g=9.81.0m / s^2),
                      Density= GasPyroclast_Density(ρmelt=ConstantDensity(ρ=2900kg/m^3), ρgas=ConstantDensity(ρ=1kg/m^3), δ=0.0, β=0.0),
                      )
```

# References
- Slezin, Yu. B. (2003), The mechanism of volcanic eruptions (a steady state approach), Journal of Volcanology and Geothermal Research, 122, 7-50, https://doi.org/10.1016/S0377-0273(02)00464-X
"""
@with_kw_noshow struct GasPyroclast_Density{_T, U1, U2, U3, S1 <: AbstractDensity, S2 <: AbstractDensity} <: ConduitDensity{_T}
    ρmelt::S1 = ConstantDensity(ρ = 2200kg / m^3)   # density of the melt
    ρgas::S2 = ConstantDensity(ρ = 1kg / m^3)      # density of the gas
    δ::GeoUnit{_T, U1} = 0.0e0 * NoUnits          # volume fraction of free gas in flow
    β::GeoUnit{_T, U2} = 0.0e0 * NoUnits          # gas volume fraction enclosed within the particles
    ρ::GeoUnit{_T, U3} = 2900.0kg / m^3         # to keep track on whether this struct is dimensional or not
end

GasPyroclast_Density(args...) = GasPyroclast_Density(args[1], args[2], convert.(GeoUnit, args[3:end])...)
isdimensional(s::GasPyroclast_Density) = isdimensional(s.ρmelt)

function param_info(s::GasPyroclast_Density) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho = \rho_g\delta + \rho_p(1 - \delta)")
end


# Calculation routines
@inline function (rho::GasPyroclast_Density{_T})(; kwargs...) where {_T}
    ρmelt = compute_density(rho.ρmelt, kwargs)
    ρgas = compute_density(rho.ρgas, kwargs)
    @unpack_val δ, β = rho

    return @muladd ρgas * δ + ρmelt * (1 - β) * (1 - δ)
end

@inline (s::GasPyroclast_Density)(args) = s(; args...)
@inline compute_density(s::GasPyroclast_Density, args) = s(args)

# Print info
function show(io::IO, g::GasPyroclast_Density)
    return print(io, "Gas-Pyroclast mixture density: ρ = ρgas*δ + ρmelt*(1-β); ρmelt=$(g.ρmelt); ρgas=$(g.ρgas); δ=$(UnitValue(g.δ)); β=$(UnitValue(g.β))")
end

# Gas equation-of-state densities -----------------------------------------
"""
    RedlichKwong_Density(; coeffs=(-112.528, 127.811, 112.04), T0=273.15K, Tref=1K, Pref=1e5Pa, ρref=1e3kg/m^3)

Modified Redlich–Kwong gas density of Huber et al. (2010), as used by Degruyter
& Huber (2014). Fitted for H2O gas over roughly ``873 < T < 1173`` K and
``30 < P < 400`` MPa:
```math
    \\rho_g = \\rho_{ref}\\left(a_1\\,\\tau^{-0.381} + a_2\\,\\varpi^{-1.135}
             + a_3\\,\\tau^{-0.411}\\varpi^{0.033}\\right),
    \\quad \\tau = \\frac{T - T_0}{T_{ref}},\\ \\varpi = \\frac{P}{P_{ref}}
```
The fit's specific units (°C, bar) enter through the reference quantities: with
`T0=273.15K`, `Tref=1K` the group ``\\tau`` equals ``T`` in °C, and with
`Pref=1e5Pa` (1 bar) ``\\varpi`` equals ``P`` in bar. Because ``\\tau`` and
``\\varpi`` are ratios of like-dimensioned quantities, the expression is
dimensionally homogeneous and nondimensionalizes cleanly: only the reference
`GeoUnit`s scale, the fitted dimensionless `coeffs` do not. The parameterisation
ignores gas composition (pseudo-pure H2O), matching the reference.

Returns `NaN` for `T ≤ 273.15K` or `P ≤ 0` (`τ ≤ 0` or `ω ≤ 0`), since
`τ^(-0.381)`/`τ^(-0.411)` and `ω^(-1.135)` require a positive base: raising to
those exponents would otherwise throw an opaque low-level `DomainError` about
complex exponentiation.

# References
- Huber C., Bachmann O. , Manga M., Two Competing Effects of Volatiles on Heat Transfer in Crystal-rich Magmas: Thermal Insulation vs Defrosting, Journal of Petrology, 847–867, https://doi.org/10.1093/petrology/egq003
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
"""
struct RedlichKwong_Density{T, U1, U2, U3} <: AbstractDensity{T}
    coeffs::NTuple{3, T}      # dimensionless fit coefficients (Huber et al. 2010)
    T0::GeoUnit{T, U1}        # Kelvin -> Celsius reference offset
    Tref::GeoUnit{T, U1}      # temperature scale of the fit (1 K)
    Pref::GeoUnit{T, U2}      # pressure scale of the fit (1 bar)
    ρref::GeoUnit{T, U3}      # density scale of the fit (1e3 kg/m³)

    function RedlichKwong_Density(;
            coeffs = (-112.528, 127.811, 112.04),
            T0 = 273.15K, Tref = 1.0K, Pref = 1.0e5Pa, ρref = 1.0e3kg / m^3
        )
        T0U = convert(GeoUnit, T0)
        TrefU = convert(GeoUnit, Tref)
        PrefU = convert(GeoUnit, Pref)
        ρrefU = convert(GeoUnit, ρref)
        T = typeof(ρrefU).types[1]
        U1 = typeof(T0U).types[2]
        U2 = typeof(PrefU).types[2]
        U3 = typeof(ρrefU).types[2]
        return new{T, U1, U2, U3}(T.(coeffs), T0U, TrefU, PrefU, ρrefU)
    end
    RedlichKwong_Density(coeffs, T0, Tref, Pref, ρref) = RedlichKwong_Density(; coeffs, T0, Tref, Pref, ρref)
end
isdimensional(s::RedlichKwong_Density) = isdimensional(s.ρref)

function param_info(s::RedlichKwong_Density)
    return MaterialParamsInfo(; Equation = L"\rho_g = \rho_{ref}(a_1 \tau^{-0.381} + a_2 \varpi^{-1.135} + a_3 \tau^{-0.411} \varpi^{0.033})")
end

@inline function (rho::RedlichKwong_Density)(; P = 0.0e0, T = 0.0e0, kwargs...)
    a1, a2, a3 = rho.coeffs
    if P isa Quantity
        @unpack_units T0, Tref, Pref, ρref = rho
    else
        @unpack_val T0, Tref, Pref, ρref = rho
    end
    τ = (T - T0) / Tref            # ∝ T in °C
    ω = P / Pref                   # ∝ P in bar
    # τ^(-0.381)/ω^(-1.135) undefined otherwise. No physically sensible clamp
    # exists for a fitted gas density, so it NaN's
    τ > 0 && ω > 0 || return ρref * oftype(ustrip(ρref), NaN)
    e1, e2, e3, e4 = oftype(a1, -0.381), oftype(a1, -1.135), oftype(a1, -0.411), oftype(a1, 0.033)
    return @muladd @pow (a1 * τ^e1 + a2 * ω^e2 + a3 * τ^e3 * ω^e4) * ρref
end
@inline (s::RedlichKwong_Density)(args) = s(; args...)
@inline compute_density(s::RedlichKwong_Density, args) = s(args)

function show(io::IO, g::RedlichKwong_Density)
    return print(io, "Redlich-Kwong gas density (Huber et al. 2010): ρgas = ρref*(a1*τ^-0.381 + a2*ω^-1.135 + a3*τ^-0.411*ω^0.033)")
end

"""
    IdealGas_Density(; Rs=461.5J/kg/K)

Ideal-gas density
```math
    \\rho_g = \\frac{P}{R_s T}
```
with specific gas constant `Rs` (default: water vapor, ``R/M_w = 8.314/0.01802``).
Unlike [`RedlichKwong_Density`](@ref) this is dimensionally homogeneous, so it
nondimensionalizes cleanly and serves as the analytic derivative check for the
Redlich–Kwong fit.
"""
@with_kw_noshow struct IdealGas_Density{_T, U1} <: AbstractDensity{_T}
    Rs::GeoUnit{_T, U1} = 461.5J / kg / K   # specific gas constant
end
IdealGas_Density(args...) = IdealGas_Density(convert.(GeoUnit, args)...)
isdimensional(s::IdealGas_Density) = isdimensional(s.Rs)

function param_info(s::IdealGas_Density)
    return MaterialParamsInfo(; Equation = L"\rho_g = P/(R_s T)")
end

@inline function (rho::IdealGas_Density)(; P = 0.0e0, T = 0.0e0, kwargs...)
    if P isa Quantity
        @unpack_units Rs = rho
        return uconvert(kg / m^3, P / (Rs * T))   # P/(Rs*T) reduces to kg*Pa/J otherwise
    end
    @unpack_val Rs = rho
    return P / (Rs * T)
end
@inline (s::IdealGas_Density)(args) = s(; args...)
@inline compute_density(s::IdealGas_Density, args) = s(args)

function show(io::IO, g::IdealGas_Density)
    return print(io, "Ideal-gas density: ρgas = P/(Rs*T); Rs=$(UnitValue(g.Rs))")
end

# Three-phase (melt + crystal + gas) mixture density ----------------------
"""
    ThreePhase_Density(; ρmelt=ConstantDensity(ρ=2300kg/m^3), ρsolid=ConstantDensity(ρ=2600kg/m^3), ρgas=RedlichKwong_Density(), ρ=2400kg/m^3)

Volume-fraction mixture density of a melt + crystal + gas magma (Degruyter &
Huber 2014):
```math
    \\rho = \\varepsilon_m \\rho_m + \\varepsilon_g \\rho_g + \\varepsilon_x \\rho_x,
    \\qquad \\varepsilon_m = 1 - \\varepsilon_x - \\varepsilon_g
```
The gas and crystal volume fractions `ϕ_gas`, `ϕ_x` are inputs (from the solver's
ODE state / a melting law), passed through `args`; the sub-densities receive the
same `args` (e.g. `P`, `T`) so an equation-of-state gas density is evaluated in
place. `ρ` is a dimensional-tracking sentinel.

`ρgas` is only evaluated when `ϕ_gas != 0`: a cell with no exsolved gas never
depends on the gas closure at all, so a strict EOS like
[`RedlichKwong_Density`](@ref) (which returns `NaN` outside its calibration
window) cannot poison a `ϕ_gas=0` cell regardless of `P`, `T`.

# Arguments
- `ρmelt`: melt-phase density
- `ρsolid`: crystal-phase density
- `ρgas`: gas-phase density (e.g. [`RedlichKwong_Density`](@ref))

# References
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
"""
@with_kw_noshow struct ThreePhase_Density{_T, U, S1 <: AbstractDensity, S2 <: AbstractDensity, S3 <: AbstractDensity} <: AbstractDensity{_T}
    ρmelt::S1 = ConstantDensity(ρ = 2300kg / m^3)   # melt
    ρsolid::S2 = ConstantDensity(ρ = 2600kg / m^3)      # crystals
    ρgas::S3 = IdealGas_Density()               # gas
    ρ::GeoUnit{_T, U} = 2400.0kg / m^3              # dimensional-tracking sentinel
end
ThreePhase_Density(args...) = ThreePhase_Density(args[1], args[2], args[3], convert.(GeoUnit, args[4:end])...)
isdimensional(s::ThreePhase_Density) = isdimensional(s.ρ)

function param_info(s::ThreePhase_Density)
    return MaterialParamsInfo(; Equation = L"\rho = \varepsilon_m \rho_m + \varepsilon_g \rho_g + \varepsilon_x \rho_x")
end

@inline function (rho::ThreePhase_Density)(; ϕ_gas = 0.0e0, ϕ_x = 0.0e0, kwargs...)
    ρmelt = compute_density(rho.ρmelt, kwargs)
    ρsolid = compute_density(rho.ρsolid, kwargs)
    # Skip the gas EOS entirely when no gas is present: a zero-weighted term
    # must not depend on ρgas's domain restrictions (e.g. RedlichKwong_Density
    # returns NaN outside its calibration window regardless of ϕ_gas).
    ρgas = iszero(ϕ_gas) ? zero(ρmelt) : compute_density(rho.ρgas, kwargs)
    ϕ_m = 1 - ϕ_x - ϕ_gas
    return @muladd ϕ_m * ρmelt + ϕ_gas * ρgas + ϕ_x * ρsolid
end
@inline (s::ThreePhase_Density)(args) = s(; args...)
@inline compute_density(s::ThreePhase_Density, args) = s(args)

function show(io::IO, g::ThreePhase_Density)
    return print(io, "Three-phase density: ρ = ϕ_m*ρmelt + ϕ_gas*ρgas + ϕ_x*ρsolid; ρmelt=$(g.ρmelt); ρsolid=$(g.ρsolid); ρgas=$(g.ρgas)")
end
#-------------------------------------------------------------------------
# Melt_DensityX density depending on Oxide composition
"""
    Melt_DensityX(; oxd_wt = oxd_wt)

Set a density depending on the oxide composition after the python script by Iacovino K & Till C (2019)

## Arguments
- `oxd_wt::NTuple{9, T}`: Melt composition as 9-element Tuple containing concentrations
             in [wt%] of the following oxides ordered in the exact sequence \n
             (SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O) \n
             Default values are for a hydrous N-MORB composition

The callable also accepts an `mH2O` keyword (melt water content, mass
fraction) that overrides the struct's own `oxd_wt[9]` for that call,
recomputing only the water-dependent aggregates. `mH2O` is accepted as
given: GeoParams does not check it for consistency with any `Solubility` or
other closure's dissolved-water output — that is the caller's/solver's
responsibility.

## References
- Iacovino K & Till C (2019). DensityX: A program for calculating the densities of magmatic liquids up to 1,627 °C and 30 kbar. Volcanica 2(1), p 1-10. [doi:10.30909/vol.02.01.0110](https://dx.doi.org/10.30909/vol.02.01.0110)
"""
struct Melt_DensityX{T, T1, T2, T3, T4, T5, T6, U, U1, U2, U3, U4, U5} <: AbstractDensity{T}
    oxd_wt::NTuple{9, T} # Oxide weight percent
    MW::NTuple{9, T1}    # Molar weights g/mol
    MV::NTuple{9, T2}    # Partial molar volumes
    dVdT::NTuple{9, T3}  # Partial molar volumes
    dVdP::NTuple{9, T4}  # Partial molar volumes
    Tref::NTuple{9, T5} # Reference temperature in K
    norm_MP::NTuple{9, T6} # Normalized molar proportions
    P0::GeoUnit{T, U}    # Pressure
    ρ0::GeoUnit{T, U1}
    α::GeoUnit{T, U2}
    β::GeoUnit{T, U3}
    sum_XMW::GeoUnit{T, U4}
    sum_Vliq::GeoUnit{T, U5}

    function Melt_DensityX(;
            oxd_wt = (50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.14, 1.0),
            MW = (0.0600855, 0.07988, 0.10196, 0.07185, 0.0403, 0.05608, 0.06198, 0.0942, 0.01802) .* (kg / mol), # Molar weights kg/mol
            MV = (2.686e-5, 2.832e-5, 3.742e-5, 1.268e-5, 1.202e-5, 1.69e-5, 2.965e-5, 4.728e-5, 2.29e-5) .* (m^3 / mol), # Partial molar volumes
            dVdT = (0.0, 7.24e-9, 2.62e-9, 3.69e-9, 3.27e-9, 3.74e-9, 7.68e-9, 1.208e-8, 9.5e-9) .* (m^3 / (mol * K)),
            dVdP = (-1.89e-15, -2.31e-15, -2.26e-15, -4.5e-16, 2.7e-16, 3.4e-16, -2.4e-15, -6.75e-15, -3.2e-15) .* (m^3 / (mol * Pa)),
            Tref = (1773.15, 1773.15, 1773.15, 1723.15, 1773.15, 1773.15, 1773.15, 1773.15, 1273.15) .* K, # Reference temperature in K
            norm_MP = (0.511, 0.012, 0.09, 0.083, 0.117, 0.123, 0.028, 0.001, 0.034) .* NoUnits,
            P0 = 1.0e5Pa,
            ρ0 = 2900.0kg / m^3,
            α = 2.266e-4 / K,
            β = 6.319e-11 / Pa,
            sum_XMW = 0.0kg / mol,
            sum_Vliq = 0.0m^3 / mol
        )

        sum_XMW, norm_MP = compute_XMW_norm_MP(oxd_wt, MW)

        MWU = ntuple(i -> convert(GeoUnit, MW[i]), Val(9))
        MVU = ntuple(i -> convert(GeoUnit, MV[i]), Val(9))
        dVdTU = ntuple(i -> convert(GeoUnit, dVdT[i]), Val(9))
        dVdPU = ntuple(i -> convert(GeoUnit, dVdP[i]), Val(9))
        TrefU = ntuple(i -> convert(GeoUnit, Tref[i]), Val(9))
        norm_MPU = ntuple(i -> convert(GeoUnit, norm_MP[i]), Val(9))
        P0U = convert(GeoUnit, P0)
        ρ0U = convert(GeoUnit, ρ0)
        αU = convert(GeoUnit, α)
        βU = convert(GeoUnit, β)
        sum_XMWU = convert(GeoUnit, sum_XMW)
        sum_VliqU = convert(GeoUnit, sum_Vliq)

        T = eltype(oxd_wt)
        T1 = eltype(MWU)
        T2 = eltype(MVU)
        T3 = eltype(dVdTU)
        T4 = eltype(dVdPU)
        T5 = eltype(TrefU)
        T6 = eltype(norm_MPU)
        U = typeof(P0U).types[2]
        U1 = typeof(ρ0U).types[2]
        U2 = typeof(αU).types[2]
        U3 = typeof(βU).types[2]
        U4 = typeof(sum_XMWU).types[2]
        U5 = typeof(sum_VliqU).types[2]


        return new{T, T1, T2, T3, T4, T5, T6, U, U1, U2, U3, U4, U5}(oxd_wt, MWU, MVU, dVdTU, dVdPU, TrefU, norm_MPU, P0U, ρ0U, αU, βU, sum_XMWU, sum_VliqU)
    end

    function Melt_DensityX(oxd_wt, MW, MV, dVdT, dVdP, Tref, norm_MP, P0, ρ0, α, β, sum_XMW, sum_Vliq)
        return Melt_DensityX(;
            oxd_wt = oxd_wt, MW = MW, MV = MV, dVdT = dVdT, dVdP = dVdP, Tref = Tref, norm_MP = norm_MP, P0 = P0, ρ0 = ρ0, α = α, β = β, sum_XMW = sum_XMW, sum_Vliq = sum_Vliq,
        )
    end
end

isdimensional(g::Melt_DensityX) = isdimensional(g.ρ0)

function param_info(s::Melt_DensityX) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho from an oxide composition")
end

function compute_XMW_norm_MP(oxd_wt, MW)

    tmp = ntuple(i -> oxd_wt[i], Val(9))
    # Normalize original wt% values to 100% sum
    norm_WP = tmp ./ sum(tmp) .* 100.0

    # Divide normalized wt% values by molecular weights
    part_MP = norm_WP ./ MW
    sum_MP = sum(part_MP)

    # Convert to mol fraction
    norm_MP = part_MP ./ sum_MP

    # Calculate partial X*MW
    part_XMW = norm_MP .* MW
    sum_XMW = sum(part_XMW)

    return sum_XMW, norm_MP
end

function (s::Melt_DensityX)(; P::Number = 0.0e0, T::Number = 0.0e0, mH2O = s.oxd_wt[9] / 100, kwargs...)
    # mH2O overrides the struct's own frozen water content for this call; the
    # water-dependent aggregates (sum_XMW, norm_MP) are only recomputed when
    # it actually differs, so the default path stays exactly as fast as
    # before this kwarg existed. mH2O is accepted as given: GeoParams does
    # not check it for consistency with any Solubility/other closure's
    # dissolved water — that is the caller's/solver's responsibility.
    P0, ρ0, sum_XMW, sum_Vliq, MV, dVdT, Tref, norm_MP, dVdP = if P isa Quantity
        (; MV, dVdT, dVdP, Tref, norm_MP) = s
        @unpack_units P0, ρ0, sum_XMW, sum_Vliq = s
        norm_MPv = unpack_units(norm_MP)
        if mH2O != s.oxd_wt[9] / 100
            oxd_wt = 100 * mH2O
            s.oxd_wt[9] = oxd_wt
            sum_XMW, norm_MPv = compute_XMW_norm_MP(oxd_wt, unpack_units(s.MW))
        end
        P0, ρ0, sum_XMW, sum_Vliq, unpack_units(MV), unpack_units(dVdT), unpack_units(Tref), norm_MPv, unpack_units(dVdP)

    else
        (; MV, dVdT, dVdP, Tref, norm_MP) = s
        @unpack_val P0, ρ0, sum_XMW, sum_Vliq = s
        norm_MPv = unpack_vals(norm_MP)
        if mH2O != s.oxd_wt[9] / 100
            oxd_wt = Base.setindex(s.oxd_wt, 100 * mH2O, 9)
            sum_XMW, norm_MPv = compute_XMW_norm_MP(oxd_wt, unpack_vals(s.MW))
        end
        P0, ρ0, sum_XMW, sum_Vliq, unpack_vals(MV), unpack_vals(dVdT), unpack_vals(Tref), norm_MPv, unpack_vals(dVdP)
    end

    sum_Vliq = @muladd (MV[1] + (dVdT[1] * (T - Tref[1])) + (dVdP[1] * (P - P0))) * norm_MP[1]
    Base.@nexprs 8 i -> sum_Vliq += @muladd (MV[i + 1] + dVdT[i + 1] * (T - Tref[i + 1]) + dVdP[i + 1] * (P - P0)) * norm_MP[i + 1]

    return sum_XMW / sum_Vliq
end

@inline (s::Melt_DensityX)(args) = s(; args...)
@inline compute_density(s::Melt_DensityX, args) = s(args)

# Print info
function show(io::IO, g::Melt_DensityX)
    return print(io, "Density from oxide composition with Oxide comp.: $(g.oxd_wt)")
end

# MAGEMin DB density -------------------------------
"""
    Vector_Density(_T)

Stores a vector with density data that can be retrieved by providing an `index`
"""
struct Vector_Density{_T, V <: AbstractVector} <: AbstractDensity{_T}
    rho::V       # Density
end
Vector_Density(; rho = Vector{Float64}()) = Vector_Density{eltype(rho), typeof(rho)}(rho)

# This assumes that density always has a single parameter. If that is not the case, we will have to extend this (to be done)
function param_info(s::Vector_Density) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho from a precomputed vector")
end

# Calculation routine
"""
    compute_density(s::Vector_Density; index::Int64, kwargs...)

Pointwise calculation of density from a vector where `index` is the index of the point
"""
@inline function (s::Vector_Density)(; index::Int64, kwargs...)
    return s.rho[index]
end

@inline (s::Vector_Density)(args) = s(; args...)
@inline compute_density(s::Vector_Density, args) = s(args)

# Print info
function show(io::IO, g::Vector_Density)
    return print(io, "Density from precomputed vector with $(length(g.rho)) entries.")
end
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Phase diagrams
function param_info(s::AbstractPhaseDiagramsStruct) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho = f_{PhaseDiagram}(T,P))")
end

# """
#     compute_density(P,T, s::AbstractPhaseDiagramsStruct)
# Interpolates density as a function of `T,P` from a lookup table
# """
# @inline function compute_density(s::AbstractPhaseDiagramsStruct; P, T, kwargs...)
#     fn = s.Rho
#     return fn(T, P)
# end
# @inline compute_density(s::AbstractPhaseDiagramsStruct, args) = compute_density(s; args...)
"""
    compute_density(P,T, s::AbstractPhaseDiagramsStruct)
Interpolates density as a function of `T,P` from a lookup table
"""
@inline function compute_density(s::AbstractPhaseDiagramsStruct; P = 0.0e0, T = 0.0e0, kwargs...)
    (; coefs, T0, dT, numT, Tmax, P0, dP, numP, Pmax) = s.Rho

    rho = interpolate_field(T0, dT, numT, Tmax, P0, dP, numP, Pmax, coefs, T, P)
    return rho
end
@inline compute_density(s::AbstractPhaseDiagramsStruct, args) = compute_density(s; args...)

"""
    compute_density!(rho::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::AbstractPhaseDiagramsStruct)
In-place computation of density as a function of `T,P`, in case we are using a lookup table.
"""
# function compute_density!(rho::AbstractArray{_T}, s::AbstractPhaseDiagramsStruct; P::AbstractArray{_T}=[zero(_T)],T::AbstractArray{_T}=[zero(_T)], kwargs...) where _T end

#------------------------------------------------------------------------------------------------------------------#
# Computational routines needed for computations with the MaterialParams structure
function compute_density(s::AbstractMaterialParamsStruct, args)
    return compute_density(s.Density[1], args)
end
#-------------------------------------------------------------------------------------------------------------

"""
    compute_density!(rho::AbstractArray{_T, ndim}, MatParam::NTuple{N,AbstractMaterialParamsStruct}, Phases::AbstractArray{_I, ndim}; P=nothing, T=nothing) where {ndim,N,_T,_I<:Integer}

In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.
# Example
```julia
julia> MatParam = (SetMaterialParams(Name="Mantle", Phase=1,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PT_Density()
                        ),
                    SetMaterialParams(Name="Crust", Phase=2,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pas)),
                        Density   = ConstantDensity(ρ=2900kg/m^3))
                  );
julia> Phases = ones(Int64,10,10);
julia> Phases[:,5:end] .= 2
julia> rho     = zeros(size(Phases))
julia> T       =  ones(size(Phases))
julia> P       =  ones(size(Phases))*10
julia> args = (P=P, T=T)
julia> compute_density!(rho, MatParam, Phases, args)
julia> rho
10×10 Matrix{Float64}:
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
```
The routine is made to minimize allocations:
```julia
julia> using BenchmarkTools
julia> @btime compute_density!(\$rho, \$MatParam, \$Phases, P=\$P, T=\$T)
    203.468 μs (0 allocations: 0 bytes)
```
_________________________________________________________________________________________________________

    compute_density!(rho::AbstractArray{_T, N}, MatParam::NTuple{K,AbstractMaterialParamsStruct}, PhaseRatios::AbstractArray{_T, M}, P=nothing, T=nothing)

In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)
"""
@inline compute_density!(args::Vararg{Any, N}) where {N} = compute_param!(compute_density, args...) #Multiple dispatch to rest of routines found in Computations.jl
@inline compute_density(args::Vararg{Any, N}) where {N} = compute_param(compute_density, args...)
@inline compute_density_ratio(args::Vararg{Any, N}) where {N} = compute_param_times_frac(compute_density, args...)

# extractor methods
for type in (ConstantDensity, PT_Density, Compressible_Density, T_Density, MeltDependent_Density, Vector_Density, BubbleFlow_Density, GasPyroclast_Density, RedlichKwong_Density, IdealGas_Density, ThreePhase_Density, Melt_DensityX)
    @extractors(type, :Density)
end

import GeoParams.get_α

function get_α(rho::MeltDependent_Density; ϕ::T = 0.0, kwargs...) where {T}
    αsolid = rho.ρsolid.α.val
    αmelt = rho.ρmelt.α.val
    return @muladd ϕ * αmelt + (1 - ϕ) * αsolid
end

get_α(rho::MeltDependent_Density, args) = get_α(rho; args...)

function get_α(rho::BubbleFlow_Density; P::T = 0.0, kwargs...) where {T}
    αmelt = rho.ρmelt.α.val
    αgas = rho.ρgas.α.val
    if P isa Quantity
        @unpack_units c0, a = rho
    else
        @unpack_val c0, a = rho
    end

    cutoff = c0^2 / a^2

    if P < cutoff
        c = a * sqrt(abs(P))
    else
        c = c0
    end

    return inv((c0 - c) / αgas + (1 - (c0 - c)) / αmelt)
end

get_α(rho::BubbleFlow_Density, args) = get_α(rho; args...)

function get_α(rho::GasPyroclast_Density; kwargs...)
    αmelt = rho.ρmelt.α.val
    αgas = rho.ρgas.α.val

    @unpack_val δ, β = rho

    return @muladd αgas * δ + αmelt * (1 - β) * (1 - δ)
end

get_α(rho::GasPyroclast_Density, args) = get_α(rho; args...)

end
