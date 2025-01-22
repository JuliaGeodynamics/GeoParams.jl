module Density

# This implements different methods to compute density
#
# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, Unitful, LaTeXStrings, MuladdMacro
using ..Units
using ..PhaseDiagrams
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, @extractors, add_extractor_functions
import ..Units: isdimensional
using ..MaterialParameters: No_MaterialParam, MaterialParamsInfo
import Base.show, GeoParams.param_info

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

@inline (ρ::ConstantDensity)(; args...) = ρ.ρ.val
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
@inline function (ρ::PT_Density)(; P::Number, T::Number, kwargs...)
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
    β::GeoUnit{_T, U2} = 1.0e-9 / Pa      # P-dependence of density
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

    return @muladd ρgas * δ + ρmelt * (1 - β)
end

@inline (s::GasPyroclast_Density)(args) = s(; args...)
@inline compute_density(s::GasPyroclast_Density, args) = s(args)

# Print info
function show(io::IO, g::GasPyroclast_Density)
    return print(io, "Gas-Pyroclast mixture density: ρ = ρgas*δ + ρmelt*(1-β); ρmelt=$(g.ρmelt); ρgas=$(g.ρgas); δ=$(UnitValue(g.δ)); β=$(UnitValue(g.β))")
end
#-------------------------------------------------------------------------

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
function param_info(s::PhaseDiagram_LookupTable) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho = f_{PhaseDiagram}(T,P))")
end

"""
    compute_density(P,T, s::PhaseDiagram_LookupTable)
Interpolates density as a function of `T,P` from a lookup table
"""
@inline function compute_density(s::PhaseDiagram_LookupTable; P, T, kwargs...)
    fn = s.Rho
    return fn(T, P)
end
@inline compute_density(s::PhaseDiagram_LookupTable, args) = compute_density(s; args...)

"""
    compute_density!(rho::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::PhaseDiagram_LookupTable)
In-place computation of density as a function of `T,P`, in case we are using a lookup table.
"""
# function compute_density!(rho::AbstractArray{_T}, s::PhaseDiagram_LookupTable; P::AbstractArray{_T}=[zero(_T)],T::AbstractArray{_T}=[zero(_T)], kwargs...) where _T end

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
for type in (ConstantDensity, PT_Density, Compressible_Density, T_Density, MeltDependent_Density, Vector_Density)
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

    return @muladd αgas * δ + αmelt * (1 - β)
end

get_α(rho::GasPyroclast_Density, args) = get_α(rho; args...)

end
