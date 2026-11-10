module HeatCapacity

# This implements different methods to specify heat capacity
#
# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Adapt
using Parameters, LaTeXStrings, Unitful
using ..Units
using ..PhaseDiagrams
using ..Density: ThreePhase_Density, compute_density
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, @extractors, add_extractor_functions, AbstractPhaseDiagramsStruct
import Base.show, GeoParams.param_info
import ..Units: isdimensional
using ..MaterialParameters: MaterialParamsInfo

abstract type AbstractHeatCapacity{T} <: AbstractMaterialParam end

export compute_heatcapacity, # calculation routines
    compute_heatcapacity!, # in-place routine
    ConstantHeatCapacity, # constant
    T_HeatCapacity_Whittington, # T-dependent heat capacity
    T_HeatCapacity_Moitra, # T-dependent glass/melt heat capacity
    T_HeatCapacity_Moitra_smoothed, # same, smooth across the glass transition
    ThreePhase_HeatCapacity, # melt + crystal + gas mixture
    Latent_HeatCapacity, # Implement latent heat by modifying heat capacity
    Vector_HeatCapacity, # Vector
    param_info

include("../Computations.jl")

# Constant Heat Capacity -------------------------------------------------------
"""
    ConstantHeatCapacity(Cp=1050J/kg/K)

Set a constant heat capacity:
```math
    Cp  = cst
```
where ``Cp`` is the thermal heat capacity [``J/kg/K``].
"""
@with_kw_noshow struct ConstantHeatCapacity{T, U} <: AbstractHeatCapacity{T}
    Cp::GeoUnit{T, U} = 1050J / kg / K                # heat capacity
end
ConstantHeatCapacity(args...) = promote_construct(ConstantHeatCapacity, args...)

function param_info(s::ConstantHeatCapacity) # info about the struct
    return MaterialParamsInfo(; Equation = L"c_p = cst")
end

# Calculation routine
# Carries units when `P` or `T` does, as ConstantDensity.
@inline function compute_heatcapacity(a::ConstantHeatCapacity; P = 0, T = 0, kwargs...)
    Tc = precision_of((; P, T, kwargs...))
    @unpack_like (P isa Quantity ? P : T) Tc Cp = a
    return Cp
end

# Print info
function show(io::IO, g::ConstantHeatCapacity)
    return print(io, "Constant heat capacity: Cp=$(UnitValue(g.Cp))")
end
#-------------------------------------------------------------------------

# Temperature dependent heat capacity -------------------------------
"""
    T_HeatCapacity_Whittington()

Sets a temperature-dependent heat capacity following the parameterization of Whittington et al. (2009), Nature:
```math
    Cp = (a + b T - c/T^2)/m
```

where ``Cp`` is the heat capacity [``J/kg/K``], and ``a,b,c`` are parameters that dependent on the temperature ``T`` [``K``]:
- a = 199.50 J/mol/K    if T<= 846 K
- a = 229.32 J/mol/K    if T> 846 K
- b = 0.0857J/mol/K^2   if T<= 846 K
- b = 0.0323J/mol/K^2   if T> 846 K
- c = 5e6J/mol*K        if T<= 846 K
- c = 47.9e-6J/mol*K    if T> 846 K
- molmass =   0.22178kg/mol

Note that this is slightly different than the equation in the manuscript, as Cp is in J/kg/K (rather than ``J/mol/K`` as in eq.3/4 of the paper)
"""
@with_kw_noshow struct T_HeatCapacity_Whittington{T, U1, U2, U3, U4, U5} <:
    AbstractHeatCapacity{T}
    # Note: the resulting curve was visually compared with Fig. 2 of the paper
    a0::GeoUnit{T, U1} = 199.5J / mol / K                # prefactor for low T       (T<= 846 K)
    a1::GeoUnit{T, U1} = 229.32J / mol / K               # prefactor for high T      (T>  846 K)
    b0::GeoUnit{T, U2} = 0.0857J / mol / K^2             # linear term for low T     (T<= 846 K)
    b1::GeoUnit{T, U2} = 0.0323J / mol / K^2             # linear term for high T    (T>  846 K)
    c0::GeoUnit{T, U3} = 5.0e6J / mol * K                  # quadratic term for low T  (T<= 846 K)
    c1::GeoUnit{T, U3} = 47.9e-6J / mol * K              # quadratic term for high T (T>  846 K)
    molmass::GeoUnit{T, U4} = 0.22178kg / mol               # average molar mass
    Tcutoff::GeoUnit{T, U5} = 846K                        # cutoff temperature
end
T_HeatCapacity_Whittington(args...) = promote_construct(T_HeatCapacity_Whittington, args...)

function param_info(s::T_HeatCapacity_Whittington) # info about the struct
    return MaterialParamsInfo(; Equation = L"c_p = (a + b*T - c/T^2)/m")
end

# Calculation routine
@inline function compute_heatcapacity(a::T_HeatCapacity_Whittington; T, kwargs...)
    Tc = precision_of(T)
    @unpack_like T Tc a0, a1, b0, b1, c0, c1, molmass, Tcutoff = a

    if T ≤ Tcutoff
        return (a0 + b0 * T - c0 / T^2) / molmass
    else
        return (a1 + b1 * T - c1 / T^2) / molmass
    end
end

# Temperature-dependent heat capacity of glass and melt (Moitra et al. 2018) ----
"""
    T_HeatCapacity_Moitra()

Temperature-dependent heat capacity of a basanite glass (``T < T_g``) and melt (``T \\ge T_g``) after
Moitra et al. (2018), eqs. 3 and 4:
```math
    c_p = a_0 + b_0 T + c_0/T^2 + d_0/\\sqrt{T} \\quad (T < T_g), \\qquad
    c_p = a_1 + b_1 T \\quad (T \\ge T_g)
```
with ``T`` in K and
- a0 = 2020 J/kg/K, b0 = -0.133 J/kg/K², c0 = 4.54e6 J K/kg, d0 = -21100 J/kg/K^(1/2)
- a1 = 2150 J/kg/K, b1 = -0.5 J/kg/K²
- Tg = 1000 K

The fit is calibrated on laser-flash data between 298 and 1273 K (their Table 3) and jumps at ``T_g``; above
1273 K the melt branch is extrapolated. ``c_p`` decreases above ``T_g`` because the sample crystallizes during the
measurement (Moitra et al. 2018, p. 3627). [`T_HeatCapacity_Moitra_smoothed`](@ref GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Moitra_smoothed) replaces the jump by a smooth
transition.

# References
- Moitra, P., Sonder, I., Valentine, G. A. (2018), Effects of size and temperature-dependent thermal conductivity on the cooling of pyroclasts in air, G-cubed 19, 3623-3636, https://doi.org/10.1029/2018GC007510
"""
@with_kw_noshow struct T_HeatCapacity_Moitra{T, U1, U2, U3, U4, U5} <: AbstractHeatCapacity{T}
    a0::GeoUnit{T, U1} = 2020.0J / kg / K               # glass (T < Tg)
    b0::GeoUnit{T, U2} = -0.133J / kg / K^2
    c0::GeoUnit{T, U3} = 4.54e6J * K / kg
    d0::GeoUnit{T, U4} = -21100.0J / kg / K^(1 // 2)
    a1::GeoUnit{T, U1} = 2150.0J / kg / K               # melt (T ≥ Tg)
    b1::GeoUnit{T, U2} = -0.5J / kg / K^2
    Tg::GeoUnit{T, U5} = 1000.0K                        # glass transition temperature
end
T_HeatCapacity_Moitra(args...) = promote_construct(T_HeatCapacity_Moitra, args...)

function param_info(s::T_HeatCapacity_Moitra)
    return MaterialParamsInfo(; Equation = L"c_p = a_0 + b_0 T + c_0/T^2 + d_0/\sqrt{T}\ (T < T_g);\ c_p = a_1 + b_1 T")
end

@inline function compute_heatcapacity(a::T_HeatCapacity_Moitra; T, kwargs...)
    Tc = precision_of(T)
    @unpack_like T Tc a0, b0, c0, d0, a1, b1, Tg = a
    return T < Tg ? a0 + b0 * T + c0 / T^2 + d0 / sqrt(T) : a1 + b1 * T
end

function show(io::IO, g::T_HeatCapacity_Moitra)
    return print(io, "T-dependent glass/melt heat capacity following Moitra et al. (2018), Tg=$(UnitValue(g.Tg))")
end

"""
    T_HeatCapacity_Moitra_smoothed()

[`T_HeatCapacity_Moitra`](@ref GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Moitra) with the jump at ``T_g`` replaced by a smooth transition of width ``\\Delta T_g``:
```math
    c_p = (1 - H)\\, c_{p,glass} + H\\, c_{p,melt}, \\qquad
    H = \\frac{1}{1 + \\exp\\left(-4 (T - T_g)/\\Delta T_g\\right)}
```
where ``c_{p,glass}`` and ``c_{p,melt}`` are the glass and melt branches of `T_HeatCapacity_Moitra` evaluated at
``T``. ``H`` is the smoothed Heaviside step of [`SmoothMelting`](@ref GeoParams.MeltingParam.SmoothMelting) with ``k = 2/\\Delta T_g``, centered at
``T_g``; it rises from 0.018 at ``T_g - \\Delta T_g`` to 0.982 at ``T_g + \\Delta T_g``, i.e. from 5 % to 95 %
over about ``1.5 \\Delta T_g``. The default ``\\Delta T_g = 20`` K keeps both branches within about 15 K of ``T_g``.
"""
@with_kw_noshow struct T_HeatCapacity_Moitra_smoothed{T, U1, U2, U3, U4, U5} <: AbstractHeatCapacity{T}
    a0::GeoUnit{T, U1} = 2020.0J / kg / K               # glass
    b0::GeoUnit{T, U2} = -0.133J / kg / K^2
    c0::GeoUnit{T, U3} = 4.54e6J * K / kg
    d0::GeoUnit{T, U4} = -21100.0J / kg / K^(1 // 2)
    a1::GeoUnit{T, U1} = 2150.0J / kg / K               # melt
    b1::GeoUnit{T, U2} = -0.5J / kg / K^2
    Tg::GeoUnit{T, U5} = 1000.0K                        # glass transition temperature
    ΔTg::GeoUnit{T, U5} = 20.0K                         # width of the transition
end
T_HeatCapacity_Moitra_smoothed(args...) = promote_construct(T_HeatCapacity_Moitra_smoothed, args...)

function param_info(s::T_HeatCapacity_Moitra_smoothed)
    return MaterialParamsInfo(; Equation = L"c_p = (1 - H) c_{p,glass} + H c_{p,melt};\ H = 1/(1 + \exp(-4 (T - T_g)/\Delta T_g))")
end

@inline function compute_heatcapacity(a::T_HeatCapacity_Moitra_smoothed; T, kwargs...)
    Tc = precision_of(T)
    @unpack_like T Tc a0, b0, c0, d0, a1, b1, Tg, ΔTg = a
    H = inv(1 + exp(-4 * (T - Tg) / ΔTg))
    return (1 - H) * (a0 + b0 * T + c0 / T^2 + d0 / sqrt(T)) + H * (a1 + b1 * T)
end

function show(io::IO, g::T_HeatCapacity_Moitra_smoothed)
    return print(io, "T-dependent glass/melt heat capacity following Moitra et al. (2018), smoothed: Tg=$(UnitValue(g.Tg)), ΔTg=$(UnitValue(g.ΔTg))")
end

# Mass-weighted heat capacity of a melt + crystal + gas magma ---------------
"""
    ThreePhase_HeatCapacity(; ρ=ThreePhase_Density(), Cpmelt=T_HeatCapacity_Moitra(), Cpsolid=T_HeatCapacity_Whittington(), Cpgas=ConstantHeatCapacity(Cp=3880J/kg/K))

Specific heat capacity of a melt + crystal + gas magma, the mass-weighted mean of the phase heat capacities
(Heap et al. 2020, eq. 3, extended from two to three phases):
```math
    c_p = \\frac{\\phi_m \\rho_m c_{p,m} + \\phi_x \\rho_x c_{p,x} + \\phi_{gas} \\rho_{gas} c_{p,gas}}
               {\\phi_m \\rho_m + \\phi_x \\rho_x + \\phi_{gas} \\rho_{gas}},
    \\qquad \\phi_m = 1 - \\phi_x - \\phi_{gas}
```
so that ``\\rho c_p``, with ``\\rho`` from [`ThreePhase_Density`](@ref GeoParams.MaterialParameters.Density.ThreePhase_Density),
is the volume-weighted sum of the phase values ``\\rho_i c_{p,i}``.

The phase densities come from the sub-laws of `ρ`, so the same `ThreePhase_Density` describes the magma in the
momentum and energy equations. The gas and crystal volume fractions `ϕ_gas`, `ϕ_x` are passed at runtime, as for
`ThreePhase_Density`; all sub-laws receive the remaining arguments (e.g. `P`, `T`). The gas laws are only evaluated
when `ϕ_gas != 0`. The result is `NaN` unless ``\\phi_x \\ge 0``, ``\\phi_{gas} \\ge 0`` and ``\\phi_x + \\phi_{gas} \\le 1``.

The default gas heat capacity is that of pure H2O used by Degruyter & Huber (2014); it ignores the `X_co2` that
[`GasMixture_Conductivity`](@ref GeoParams.MaterialParameters.Conductivity.GasMixture_Conductivity) takes. The default
`Cpmelt` decreases above ``T_g`` because the sample of Moitra et al. (2018) crystallizes (see [`T_HeatCapacity_Moitra`](@ref)).

# References
- Heap, M. J., et al. (2020), The thermal properties of porous andesite, JVGR 398, 106901, https://doi.org/10.1016/j.jvolgeores.2020.106901
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
"""
struct ThreePhase_HeatCapacity{T, D <: ThreePhase_Density{T}, S1 <: AbstractHeatCapacity, S2 <: AbstractHeatCapacity, S3 <: AbstractHeatCapacity} <: AbstractHeatCapacity{T}
    ρ::D            # phase densities
    Cpmelt::S1
    Cpsolid::S2
    Cpgas::S3
end
function ThreePhase_HeatCapacity(;
        ρ = ThreePhase_Density(), Cpmelt = T_HeatCapacity_Moitra(), Cpsolid = T_HeatCapacity_Whittington(),
        Cpgas = ConstantHeatCapacity(Cp = 3880.0J / kg / K)
    )
    return ThreePhase_HeatCapacity(ρ, Cpmelt, Cpsolid, Cpgas)
end

function param_info(s::ThreePhase_HeatCapacity)
    return MaterialParamsInfo(; Equation = L"c_p = \sum_i \phi_i \rho_i c_{p,i} / \sum_i \phi_i \rho_i")
end

@inline function compute_heatcapacity(a::ThreePhase_HeatCapacity; ϕ_gas = 0, ϕ_x = 0, kwargs...)
    ρm, ρx = compute_density(a.ρ.ρmelt, kwargs), compute_density(a.ρ.ρsolid, kwargs)
    cpm, cpx = compute_heatcapacity(a.Cpmelt, kwargs), compute_heatcapacity(a.Cpsolid, kwargs)
    ϕ_m = 1 - ϕ_x - ϕ_gas
    (ϕ_x >= 0 && ϕ_gas >= 0 && ϕ_m >= 0) || return cpm * precision_of(cpm)(NaN)
    ρcp = ϕ_m * ρm * cpm + ϕ_x * ρx * cpx
    ρ = ϕ_m * ρm + ϕ_x * ρx
    iszero(ϕ_gas) && return ρcp / ρ
    ρg = compute_density(a.ρ.ρgas, kwargs)
    return (ρcp + ϕ_gas * ρg * compute_heatcapacity(a.Cpgas, kwargs)) / (ρ + ϕ_gas * ρg)
end

function show(io::IO, g::ThreePhase_HeatCapacity)
    return print(io, "Three-phase mass-weighted heat capacity: Cpmelt=$(g.Cpmelt); Cpsolid=$(g.Cpsolid); Cpgas=$(g.Cpgas)")
end

# LatentHeat by modifying heat capacity  ---------------------------------
"""
    Latent_HeatCapacity(Cp=ConstantHeatCapacity(), Q_L=400kJ/kg)

This takes the effects of latent heat into account by modifying the heat capacity in the temperature equation:

```math
\\rho C_p^{\\textrm{new}} \\frac{\\partial T}{\\partial t}  = \\frac{\\partial }{\\partial x_i} \\left( k \\frac{\\partial T}{\\partial x_i} \\right)  + H_s
```

with
```math
C_p^{\\textrm{new}}  = C_p + \\frac{\\partial \\phi}{\\partial T} Q_L
```
where ``Q_L`` is the latent heat [``J/kg``], and ``\\frac{\\partial \\phi}{\\partial T}`` is the derivative of the melt fraction with respect to temperature

"""
@with_kw_noshow struct Latent_HeatCapacity{T, U, S1 <: AbstractHeatCapacity} <: AbstractHeatCapacity{T}
    Cp::S1 = ConstantHeatCapacity()
    Q_L::GeoUnit{T, U} = 400.0e3J / kg                            # Latent heat
end
Latent_HeatCapacity(args...) = promote_construct(Latent_HeatCapacity, args...)
Latent_HeatCapacity(Cp::AbstractHeatCapacity, args...) = Latent_HeatCapacity(Cp, convert.(GeoUnit, args)...)
isdimensional(g::Latent_HeatCapacity) = isdimensional(g.Q_L)

function param_info(s::Latent_HeatCapacity) # info about the struct
    return MaterialParamsInfo(; Equation = L"Q_L = cst; c_p = Cp + \frac{\partial \phi}{\partial T}")
end

# Calculation routine
@inline function compute_heatcapacity(
        a::Latent_HeatCapacity{_T};
        dϕdT = 0,
        kwargs...
    ) where {_T}
    Tc = precision_of((; dϕdT, kwargs...))
    @unpack_val Tc Q_L = a

    Cp = compute_heatcapacity(a.Cp, kwargs)

    Cp += Q_L * dϕdT      # latent heat contribution

    return Cp
end

# Print info
function show(io::IO, g::Latent_HeatCapacity)
    return print(io, "Latent heat through modifying heat capacity: Cp = Cp + dϕdT*$(Value(g.Q_L)); Cp=$(g.Cp)")
end
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
"""
    Vector_HeatCapacity(_T)

Stores a vector with heat capacity data that can be retrieved by providing an `index`
"""
struct Vector_HeatCapacity{_T, V <: AbstractVector} <: AbstractHeatCapacity{_T}
    Cp::V       # Heat capacity
end
Adapt.@adapt_structure Vector_HeatCapacity
Vector_HeatCapacity(Cp::AbstractVector) = Vector_HeatCapacity{eltype(Cp), typeof(Cp)}(Cp)
Vector_HeatCapacity(; Cp = Vector{Float64}()) = Vector_HeatCapacity(Cp)

function param_info(s::Vector_HeatCapacity) # info about the struct
    return MaterialParamsInfo(; Equation = L"Cp from a precomputed vector")
end

"""
    compute_heatcapacity(a::Vector_HeatCapacity; index::Int64, kwargs...)

Pointwise calculation of heat capacity from a vector where `index` is the index of the point
"""
@inline function compute_heatcapacity(s::Vector_HeatCapacity; index::Int64 = 1, kwargs...)
    return s.Cp[index]
end

# Print info
function show(io::IO, g::Vector_HeatCapacity)
    return print(io, "Heat capacity from precomputed vector with $(length(g.Cp)) entries")
end

#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
"""
    compute_heatcapacity!(Cp_array::AbstractArray{_T, N},s::T_HeatCapacity_Whittington{_T}, T::_T=zero(_T), P::_T=zero(_T)) where {_T,N}

Computes T-dependent heat capacity in-place
"""
# function compute_heatcapacity!(Cp_array::AbstractArray{_T, N},s::T_HeatCapacity_Whittington{_T}; T::AbstractArray{_T, N}, kwargs...) where {_T,N} end

# add methods programmatically
for myType in (:ConstantHeatCapacity, :T_HeatCapacity_Whittington, :T_HeatCapacity_Moitra, :T_HeatCapacity_Moitra_smoothed, :ThreePhase_HeatCapacity, :Latent_HeatCapacity, :Vector_HeatCapacity)
    @eval begin
        #(s::$(myType))(args) = s(; args...)
        compute_heatcapacity(s::$(myType), args) = compute_heatcapacity(s; args...)
    end
end

# Print info
function show(io::IO, g::T_HeatCapacity_Whittington)
    return print(
        io,
        "T-dependent heat capacity following Whittington et al. (2009) for average crust. \n",
    )
end
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# Heat capacity from phase diagram
"""
    compute_heatcapacity(P,T, s::AbstractPhaseDiagramsStruct)
Interpolates heat capacity as a function of `T,P` from a lookup table
"""
@inline compute_heatcapacity(s::AbstractPhaseDiagramsStruct; P, T, kwargs...) = s.SpecificCp(T, P)
@inline compute_heatcapacity(s::AbstractPhaseDiagramsStruct, args) = compute_heatcapacity(s; args...)

#-------------------------------------------------------------------------

# Help info for the calculation routines
"""
    Cp = compute_heatcapacity(s::AbstractHeatCapacity, P, T)

Returns the heat capacity `Cp` at any temperature `T` and pressure `P` using any of the heat capacity laws implemented.

Currently available:
- ConstantHeatCapacity
- T\\_HeatCapacity_Whittington
- T\\_HeatCapacity_Moitra
- T\\_HeatCapacity_Moitra\\_smoothed
- ThreePhase\\_HeatCapacity
- Latent_HeatCapacity

# Example
```jldoctest
julia> cp = T_HeatCapacity_Whittington()
T-dependent heat capacity following Whittington et al. (2009) for average crust.

julia> compute_heatcapacity(cp, (; T = 1000.0))
1179.6374785821629
```

The in-place [`compute_heatcapacity!`](@ref) fills a preallocated array instead:

```jldoctest
julia> cp = T_HeatCapacity_Whittington();

julia> T = collect(250.0:250.0:1250.0);

julia> Cp = similar(T);

julia> compute_heatcapacity!(Cp, cp, (; T = T));

julia> Cp
5-element Vector{Float64}:
  635.4269997294616
 1002.5701145279105
 1149.2745563671706
 1179.6374785821629
 1216.0474343943067
```


"""
compute_heatcapacity()

"""
    Cp = ComputeHeatCapacity(T::Any, s::AbstractHeatCapacity)

Computes heat capacity if only temperature (and not pressure) is specified
"""
# compute_heatcapacity(s::AbstractHeatCapacity, T::AbstractArray{_T}) where _T =  compute_heatcapacity(s,similar(T), T)
# compute_heatcapacity!(Cp_array::AbstractArray{_T}, s::AbstractHeatCapacity, T::AbstractArray{_T}) where _T =  compute_heatcapacity!(Cp_array,s,similar(T), T)

"""
    Cp = ComputeHeatCapacity(s::ConstantHeatCapacity)

Returns heat capacity if we are sure that we will only employ constant heat capacity in the simulation
"""

# Computational routines needed for computations with the MaterialParams structure
function compute_heatcapacity(s::AbstractMaterialParamsStruct, args)
    return compute_heatcapacity(s.HeatCapacity[1], args)
end

"""
    compute_heatcapacity!(Cp::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat})

In-place computation of heat capacity `Cp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.
"""
compute_heatcapacity!()

compute_heatcapacity(MatParam, arg, args::Vararg{Any, N}) where {N} = compute_param(compute_heatcapacity, MatParam, arg, args...)
compute_heatcapacity!(args::Vararg{Any, N}) where {N} = compute_param!(compute_heatcapacity, args...)

# extractor methods
for type in (ConstantHeatCapacity, T_HeatCapacity_Whittington, T_HeatCapacity_Moitra, T_HeatCapacity_Moitra_smoothed, ThreePhase_HeatCapacity, Latent_HeatCapacity, Vector_HeatCapacity)
    @extractors(type, :HeatCapacity)
end

end
