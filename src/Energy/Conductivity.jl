module Conductivity

# This implements different methods to specify conductivity of rocks
#
# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using ..Density: AbstractDensity, ConstantDensity, HalbachChatterjee_Density, compute_density
using ..HeatCapacity: AbstractHeatCapacity, T_HeatCapacity_Moitra, T_HeatCapacity_Moitra_smoothed, compute_heatcapacity
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, @extractors, add_extractor_functions
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractConductivity{T} <: AbstractMaterialParam end

export compute_conductivity, # calculation routines
    compute_conductivity!,
    param_info,
    ConstantConductivity, # constant
    T_Conductivity_Whittington, # T-dependent heat capacity
    T_Conductivity_Whittington_parameterised, # parameterised T-dependent heat capacity
    TP_Conductivity, # TP dependent conductivity
    T_Conductivity_Moitra, # from T-dependent glass/melt diffusivity
    T_Conductivity_Moitra_smoothed, # same, smooth across the glass transition
    IAPWS_H2O_Conductivity, # H2O, IAPWS 2011
    GasMixture_Conductivity, # H2O-CO2 gas mixture
    Maxwell_Conductivity, # matrix with gas bubbles
    Set_TP_Conductivity                         # Routine to set pre-defined parameters

include("../Computations.jl")

# Constant Conductivity -------------------------------------------------------
"""
    ConstantConductivity(k=3.0W/m/K)

Set a constant conductivity
```math
    k  = cst
```
where ``k`` is the thermal conductivity [``W/m/K``].
"""
@with_kw_noshow struct ConstantConductivity{T, U} <: AbstractConductivity{T}
    k::GeoUnit{T, U} = 3.0Watt / m / K
end
ConstantConductivity(args...) = promote_construct(ConstantConductivity, args...)

function param_info(s::ConstantConductivity) # info about the struct
    return MaterialParamsInfo(; Equation = L"k = cst")
end

# Calculation routine; carries units when `P` or `T` does, as ConstantDensity.
function (s::ConstantConductivity)(; P = 0, T = 0, kwargs...)
    Tc = precision_of((; P, T, kwargs...))
    @unpack_like (P isa Quantity ? P : T) Tc k = s

    return k
end

compute_conductivity(s::ConstantConductivity; kwargs...) = s(; kwargs...)

function (s::ConstantConductivity)(I::Integer...)
    @unpack_val k = s

    return fill(k, I...)
end

compute_conductivity(s::ConstantConductivity, I::Integer...) = s(I...)

"""
    compute_conductivity(k_array::AbstractArray{<:AbstractFloat,N},P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, s::ConstantConductivity) where N

In-place routine to compute constant conductivity
"""
function compute_conductivity!(k_array::AbstractArray, s::ConstantConductivity; kwargs...)
    Tc = precision_of(k_array)
    @unpack_val Tc k = s
    k_array .= k
    return nothing
end

# Print info
function show(io::IO, g::ConstantConductivity)
    return print(io, "Constant conductivity: k=$(g.k.val)")
end
#-------------------------------------------------------------------------

# Temperature dependent conductivity -------------------------------
"""
    T_Conductivity_Whittington()

Sets a temperature-dependent conductivity following the parameterization of *Whittington, A.G., Hofmeister, A.M., Nabelek, P.I., 2009. Temperature-dependent thermal diffusivity of the Earth’s crust and implications for magmatism. Nature 458, 319–321. https://doi.org/10.1038/nature07818.*
Their parameterization is originally given for the thermal diffusivity, together with a parameterization for thermal conductivity, which allows us to compute
```math
    Cp = a + b T - c/T^2
```
```math
    \\kappa = d/T - e \\textrm{ if } T<=846K
```
```math
    \\kappa = f - g*T \\textrm{ if } T>846K
```
```math
    \\rho = 2700 kg/m^3
```
```math
    k = \\kappa \\rho Cp
```

where ``Cp`` is the heat capacity [``J/mol/K``], and ``a,b,c`` are parameters that dependent on the temperature `T`:
- a = 199.50 J/mol/K    if T<= 846 K
- a = 199.50 J/mol/K    if T> 846 K
- b = 0.0857J/mol/K^2   if T<= 846 K
- b = 0.0323J/mol/K^2   if T> 846 K
- c = 5e6J/mol*K        if T<= 846 K
- c = 47.9e-6J/mol*K    if T> 846 K
- d = 567.3m^2/s*K
- e = 0.062m^2/s
- f = 0.732m^2/s
- g = 0.000135m^2/s/K

This looks like:

![subet1](./assets/img/Conductivity_Whittington.png)

# Example
```julia
julia> using GLMakie, GeoParams
julia> p=T_Conductivity_Whittington();
julia> T,k,plt = PlotConductivity(p)
```


"""
@with_kw_noshow struct T_Conductivity_Whittington{T, U1, U2, U3, U4, U5, U6, U7, U8, U9} <:
    AbstractConductivity{T}
    # Note: the resulting curve of k was visually compared with Fig. 2 of the paper
    a0::GeoUnit{T, U1} = 199.5J / mol / K                # prefactor for low T       (T<= 846 K)
    a1::GeoUnit{T, U1} = 229.32J / mol / K               # prefactor for high T      (T>  846 K)
    b0::GeoUnit{T, U2} = 0.0857J / mol / K^2             # linear term for low T     (T<= 846 K)
    b1::GeoUnit{T, U2} = 0.0323J / mol / K^2             # linear term for high T    (T>  846 K)
    c0::GeoUnit{T, U3} = 5.0e6J / mol * K                  # quadratic term for low T  (T<= 846 K)
    c1::GeoUnit{T, U3} = 47.9e-6J / mol * K              # quadratic term for high T (T>  846 K)
    molmass::GeoUnit{T, U4} = 0.22178kg / mol               # average molar mass
    Tcutoff::GeoUnit{T, U5} = 846.0K                      # cutoff temperature
    rho::GeoUnit{T, U6} = 2700kg / m^3                  # Density they use for an average crust
    d::GeoUnit{T, U7} = 567.3 * 1.0e-6m^2 / s * K           # diffusivity parameterization
    e::GeoUnit{T, U8} = 0.062 * 1.0e-6m^2 / s             # diffusivity parameterization
    f::GeoUnit{T, U8} = 0.732 * 1.0e-6m^2 / s             # diffusivity parameterization
    g::GeoUnit{T, U9} = 0.000135 * 1.0e-6m^2 / s / K        # diffusivity parameterization
end
T_Conductivity_Whittington(args...) = promote_construct(T_Conductivity_Whittington, args...)

function param_info(s::T_Conductivity_Whittington) # info about the structwhere {_T,N}
    return MaterialParamsInfo(; Equation = L"k = f(T) ")
end

# Calculation routine
function (s::T_Conductivity_Whittington)(; T, kwargs...)
    Tc = precision_of(T)
    @unpack_like T Tc a0, a1, b0, b1, c0, c1, molmass, Tcutoff, rho, d, e, f, g = s

    if T ≤ Tcutoff
        return (a0 + b0 * T - c0 / T^2) / molmass * (d / T - e) * rho
    else
        return (a1 + b1 * T - c1 / T^2) / molmass * (f - g * T) * rho
    end
end

function compute_conductivity(s::T_Conductivity_Whittington; T, kwargs...)
    return s(; T = T)
end

function (s::T_Conductivity_Whittington)(T::AbstractArray; kwargs...)
    Tc = precision_of(T)
    if eltype(T) <: Quantity   # array of Quantities is not itself a Quantity
        @unpack_units Tc a0, a1, b0, b1, c0, c1, molmass, Tcutoff, rho, d, e, f, g = s
    else
        @unpack_val Tc a0, a1, b0, b1, c0, c1, molmass, Tcutoff, rho, d, e, f, g = s
    end

    inv_molmass = 1 / molmass # multiplication is considerably faster than division
    # `map` lets the output element type follow the computation, so both the unitless
    # and the dimensional (Quantity array) paths return a correctly-typed array.
    return map(T) do Ti
        if Ti <= Tcutoff
            a, b, c = a0, b0, c0
            κ = d / Ti - e
        else
            a, b, c = a1, b1, c1
            κ = f - g * Ti
        end
        Cp = (a + b * Ti - c / Ti^2) * inv_molmass # conductivity
        κ * rho * Cp       # compute conductivity from diffusivity
    end
end

(s::T_Conductivity_Whittington)(T::AbstractArray, args...) = s(T; args...)
function compute_conductivity(s::T_Conductivity_Whittington, T::AbstractArray, args...)
    return s(T; args...)
end

"""
    compute_conductivity!(k_array::AbstractArray{<:AbstractFloat,N},P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, s::T_Conductivity_Whittington) where N

In-place routine to compute temperature-dependent conductivity
"""
# function compute_conductivity!(k::AbstractArray{_T,N}, s::T_Conductivity_Whittington{_T}; T::AbstractArray{_T,N}, kwargs...) where {_T,N} end

# Print info
function show(io::IO, g::T_Conductivity_Whittington) #info about the struct
    return print(
        io,
        "T-dependent conductivity following the original Whittington et al. (2009) formulation for average crust). \n",
    )
end
#-------------------------------------------------------------------------

# Temperature dependent conductivity -------------------------------
"""
    T_Conductivity_Whittington_parameterised()

Sets a temperature-dependent conductivity that is  parameterization after *Whittington, et al.  2009*

The original parameterization involves quite a few parameters; this is a polynomial fit that is roughly valid from 0-1000 Celsius
```math
    k [W/m/K] = -2 10^{-9} (T-Ts)^3 + 6 10^{-6} (T-Ts)^2 - 0.0062 (T-Ts) + 4
```
```math
    Ts = 273.15 K
```
where `T[K]` is the temperature in Kelvin (or the nondimensional equivalent of it).

The comparison of this parameterisation vs. the original one is:
![subet1](./assets/img/Conductivity_Whittington_smooth.png)


"""
@with_kw_noshow struct T_Conductivity_Whittington_parameterised{T, U1, U2, U3, U4, U5} <:
    AbstractConductivity{T}
    # Note: the resulting curve of k was visually compared with Fig. 2 of the paper
    a::GeoUnit{T, U1} = -2.0e-9Watt / m / K^4                #
    b::GeoUnit{T, U2} = 6.0e-6Watt / m / K^3              #
    c::GeoUnit{T, U3} = -0.0062Watt / m / K^2              #
    d::GeoUnit{T, U4} = 4Watt / m / K                #
    Ts::GeoUnit{T, U5} = 273.15 * K                    # Temperature shift [C->K]
end
function T_Conductivity_Whittington_parameterised(args...)
    return T_Conductivity_Whittington_parameterised(convert.(GeoUnit, args)...)
end

function param_info(s::T_Conductivity_Whittington_parameterised) # info about the struct
    return MaterialParamsInfo(;
        Equation = L"k = a*(T[K]-Ts)^3 + b*(T[K]-Ts)^2  + c*(T[K]-Ts) + d, Ts=273.15"
    )
end

# Calculation routine
function (s::T_Conductivity_Whittington_parameterised)(;
        T = 0.0e0, kwargs...
    )
    Tc = precision_of(T)
    @unpack_like T Tc a, b, c, d, Ts = s

    T_C = T - Ts
    return a * T_C^3 + b * T_C^2 + c * T_C + d
end

function (s::T_Conductivity_Whittington_parameterised)(T::AbstractArray; kwargs...)
    if eltype(T) <: Quantity   # array of Quantities is not itself a Quantity
        @unpack_units a, b, c, d, Ts = s
    else
        @unpack_val a, b, c, d, Ts = s
    end

    # `map` lets the output element type follow the computation (unitless or Quantity).
    return map(T) do Ti
        # Note: in general, we operate with SI units or the non-dimensional equivalent of that
        T_C = Ti - Ts
        a * T_C^3 + b * T_C^2 + c * T_C + d
    end
end

(s::T_Conductivity_Whittington_parameterised)(T::AbstractArray, args) = s(T; args...)
function compute_conductivity(
        s::T_Conductivity_Whittington_parameterised, T::AbstractArray, args
    )
    return s(T, args...)
end

"""
    compute_conductivity!(k::AbstractArray{_T,N}, s::T_Conductivity_Whittington_parameterised{_T}, P::AbstractArray{_T,N}, T::AbstractArray{_T,N})

In-place routine to compute temperature-dependent conductivity
"""
# function compute_conductivity!(k::AbstractArray{_T,N}, s::T_Conductivity_Whittington_parameterised{_T}; T::AbstractArray{_T,N}, kwargs...) where {_T,N} end

# Print info
function show(io::IO, g::T_Conductivity_Whittington_parameterised) #info about the struct
    return print(
        io,
        "T-dependent conductivity parameterized from Whittington et al. (2009): k=$(NumValue(g.a))*(T-Ts)^3 + $(NumValue(g.b))*(T-Ts)^2 + $(NumValue(g.c))*(T-Ts) + $(NumValue(g.d)) with T in K and Ts=$(Value(g.Ts))\n",
    )
end
#-------------------------------------------------------------------------

# Temperature (& Pressure) dependent conductivity -------------------------------
"""
    TP_Conductivity()

Sets a temperature (and pressure)-dependent conductivity parameterization as described in Gerya, Numerical Geodynamics (2nd edition, Table 21.2).
The general for

```math
    k = \\left( a_k +  \\frac{b_k}{T + c_k} \\right) (1 + d_k P)
```

where ``k`` is the conductivity [``W/K/m``], and ``a_k,b_k,c_k,d_k`` are parameters that dependent on the temperature `T` and pressure `P`:
- ``a_k`` = 1.18Watt/K/m
- ``b_k`` = 474Watt/m
- ``c_k`` = 77K
- ``d_k`` = 0/MPa
"""
@with_kw_noshow struct TP_Conductivity{T, N, U1, U2, U3, U4} <: AbstractConductivity{T}
    Name::NTuple{N, Char} = ""            # The name is encoded as a NTuple{Char} to make it isbits
    a::GeoUnit{T, U1} = 1.18Watt / K / m  # empirical fitting term
    b::GeoUnit{T, U2} = 474.0Watt / m     # empirical fitting term
    c::GeoUnit{T, U3} = 77.0K             # empirical fitting term
    d::GeoUnit{T, U4} = 0.0 / MPa         # empirical fitting term
end

function TP_Conductivity(name::AbstractString, a, b, c, d)
    vals = convert.(GeoUnit, (a, b, c, d))
    return TP_Conductivity(Tuple(name), vals...)
end

function TP_Conductivity(name::NTuple{N, Char}, a, b, c, d) where {N}
    ga, gb, gc, gd = convert.(GeoUnit, (a, b, c, d))
    Tprom = promote_type(typeof(NumValue(ga)), typeof(NumValue(gb)), typeof(NumValue(gc)), typeof(NumValue(gd)))
    vals = (
        convert(GeoUnit{Tprom}, UnitValue(ga)),
        convert(GeoUnit{Tprom}, UnitValue(gb)),
        convert(GeoUnit{Tprom}, UnitValue(gc)),
        convert(GeoUnit{Tprom}, UnitValue(gd)),
    )
    return TP_Conductivity(name, vals...)
end

function param_info(s::TP_Conductivity)
    name = String(collect(s.Name))
    eq = L"k = \left(a_k + {b_k/{T + c_k}} \right)*(1 + d_k*P) "
    if name == ""
        return MaterialParamsInfo(; Equation = eq)
    end
    return MaterialParamsInfo(; Equation = eq, Comment = TP_Conductivity_info[name][2].Comment)
end

"""
    Set_TP_Conductivity(name::String)

Returns the pre-defined temperature(-pressure) dependent conductivity law `name`.

Available laws:
- "UpperCrust"
- "LowerCrust"
- "OceanicCrust"
- "Mantle"

# Example
```jldoctest
julia> k = Set_TP_Conductivity("Mantle")
T/P dependent conductivity: Name = Mantle, k = (0.73 + 1293.0/(T + 77.0))*(1 + 4.0e-5*P)
```
"""
Set_TP_Conductivity(name::String) = TP_Conductivity_info[name][1]

TP_Conductivity_info = Dict(
    [
        (
            "UpperCrust",
            (
                TP_Conductivity(;
                    Name = "UpperCrust", a = 0.64Watt / K / m, b = 807Watt / m, c = 77K, d = 0 / MPa
                ),
                MaterialParamsInfo(;
                    Comment = "Sediment/upper crust T-dependent conductivity, as listed in table 21.2 of Gerya et al. | Reference still to be verified!",
                ),
            ),
        )
        (
            "LowerCrust",
            (
                TP_Conductivity(;
                    Name = "LowerCrust", a = 1.18Watt / K / m, b = 474Watt / m, c = 77K, d = 0 / MPa
                ),
                MaterialParamsInfo(;
                    Comment = "Lower crust T-dependent conductivity, as listed in table 21.2 of Gerya et al. | Reference still to be verified!",
                ),
            ),
        )
        (
            "OceanicCrust",
            (
                TP_Conductivity(;
                    Name = "OceanicCrust", a = 1.18Watt / K / m, b = 474Watt / m, c = 77K, d = 0 / MPa
                ),
                MaterialParamsInfo(;
                    Comment = "Oceanic crust T-dependent conductivity, as listed in table 21.2 of Gerya et al. | Reference still to be verified!",
                ),
            ),
        )
        (
            "Mantle",
            (
                TP_Conductivity(;
                    Name = "Mantle",
                    a = 0.73Watt / K / m,
                    b = 1293Watt / m,
                    c = 77K,
                    d = 0.00004 / MPa,
                ),
                MaterialParamsInfo(;
                    Comment = "Mantle T-dependent conductivity, as listed in table 21.2 of Gerya et al. | Reference still to be verified!",
                ),
            ),
        )
    ],
)

# Calculation routine
function (s::TP_Conductivity)(; P = 0.0e0, T = 0.0e0, kwargs...)
    Tc = precision_of(P)
    T = convert_precision(Tc, T)
    @unpack_like T Tc a, b, c, d = s

    if ustrip(d) == 0
        return a + b / (T + c)
    else
        return (a + b / (T + c)) * (1 + d * P)
    end
end


function (s::TP_Conductivity)(
        P::AbstractArray, T::AbstractArray; kwargs...
    )
    Tc = precision_of(T)
    if eltype(T) <: Quantity   # array of Quantities is not itself a Quantity
        @unpack_units Tc a, b, c, d = s
    else
        @unpack_val Tc a, b, c, d = s
    end

    d_is_zero = ustrip(d) == 0
    if !d_is_zero && size(T) != size(P)
        error("Size of P and T arrays should be the same")
    end

    # `map` lets the output element type follow the computation (unitless or Quantity).
    return map(eachindex(T)) do i
        base = a + b / (T[i] + c)
        d_is_zero ? base : base * (1 + d * P[i])
    end
end


#(s::TP_Conductivity)(P::AbstractArray, T::AbstractArray) = s(P, T)
compute_conductivity(s::TP_Conductivity, P::AbstractArray, T::AbstractArray) = s(P, T)

# Conductivity from the thermal diffusivity of glass and melt (Moitra et al. 2018) ----
"""
    T_Conductivity_Moitra(; ρ=ConstantDensity(ρ=3000kg/m^3), Cp=T_HeatCapacity_Moitra())

Thermal conductivity ``k = \\rho c_p \\alpha`` (Moitra et al. 2018, eq. 2) from the temperature-dependent thermal
diffusivity of a basanite glass (``T < T_g``) and melt (``T \\ge T_g``), their eqs. 5 and 6:
```math
    \\alpha = a_0 T^{-0.41} + b_0 T \\quad (T < T_g), \\qquad \\alpha = a_1 + b_1 T \\quad (T \\ge T_g)
```
with ``T`` in K and
- a0 = 5.98 mm²/s K^0.41, b0 = 9.23e-5 mm²/s/K
- a1 = -0.06 mm²/s, b1 = 4.57e-4 mm²/s/K
- Tg = 1000 K

The density `ρ` and heat capacity `Cp` are laws evaluated with the same arguments, so ``k`` follows the material of
the model while keeping the measured diffusivity. The defaults are the sample density and heat capacity of Moitra et
al. (2018), with which ``k`` reproduces their Table 3 (298-1273 K). ``k`` jumps at ``T_g``;
[`T_Conductivity_Moitra_smoothed`](@ref GeoParams.MaterialParameters.Conductivity.T_Conductivity_Moitra_smoothed)
replaces the jump by a smooth transition. Above 1273 K the melt branches are extrapolated. Moitra et al. (2018,
p. 3627) attribute the rise of ``\\alpha`` above ``T_g`` to crystallization of the sample, so the melt branch is not
that of a crystal-free melt: at 1473 K it gives ``k \\approx 2.6`` W/m/K, about twice the 1-1.5 W/m/K typical of
basaltic melt.

# References
- Moitra, P., Sonder, I., Valentine, G. A. (2018), Effects of size and temperature-dependent thermal conductivity on the cooling of pyroclasts in air, G-cubed 19, 3623-3636, https://doi.org/10.1029/2018GC007510
"""
@with_kw_noshow struct T_Conductivity_Moitra{T, U1, U2, U3, U4, S1 <: AbstractDensity, S2 <: AbstractHeatCapacity} <: AbstractConductivity{T}
    ρ::S1 = ConstantDensity(ρ = 3000.0kg / m^3)
    Cp::S2 = T_HeatCapacity_Moitra()
    a0::GeoUnit{T, U1} = 5.98e-6m^2 / s * K^(41 // 100)   # glass (T < Tg)
    b0::GeoUnit{T, U2} = 9.23e-11m^2 / s / K
    a1::GeoUnit{T, U3} = -6.0e-8m^2 / s                    # melt (T ≥ Tg)
    b1::GeoUnit{T, U2} = 4.57e-10m^2 / s / K
    Tg::GeoUnit{T, U4} = 1000.0K                           # glass transition temperature
end
T_Conductivity_Moitra(args...) = promote_construct(T_Conductivity_Moitra, args...)

function param_info(s::T_Conductivity_Moitra)
    return MaterialParamsInfo(; Equation = L"k = \rho c_p \alpha;\ \alpha = a_0 T^{-0.41} + b_0 T\ (T < T_g);\ \alpha = a_1 + b_1 T")
end

function (s::T_Conductivity_Moitra)(; T, kwargs...)
    Tc = precision_of(T)
    @unpack_like T Tc a0, b0, a1, b1, Tg = s
    α = T < Tg ? a0 * T^(-41 // 100) + b0 * T : a1 + b1 * T
    return compute_density(s.ρ, (; T, kwargs...)) * compute_heatcapacity(s.Cp, (; T, kwargs...)) * α
end

function show(io::IO, g::T_Conductivity_Moitra)
    return print(io, "Conductivity from the glass/melt diffusivity of Moitra et al. (2018): ρ=$(g.ρ); Cp=$(g.Cp)")
end

"""
    T_Conductivity_Moitra_smoothed(; ρ=ConstantDensity(ρ=3000kg/m^3), Cp=T_HeatCapacity_Moitra_smoothed())

[`T_Conductivity_Moitra`](@ref GeoParams.MaterialParameters.Conductivity.T_Conductivity_Moitra) with the jump of the
diffusivity at ``T_g`` replaced by a smooth transition of width ``\\Delta T_g``:
```math
    k = \\rho c_p \\alpha, \\qquad \\alpha = (1 - H)\\, \\alpha_{glass} + H\\, \\alpha_{melt}, \\qquad
    H = \\frac{1}{1 + \\exp\\left(-4 (T - T_g)/\\Delta T_g\\right)}
```
with the glass and melt diffusivities of `T_Conductivity_Moitra` and the same ``H`` as
[`T_HeatCapacity_Moitra_smoothed`](@ref GeoParams.MaterialParameters.HeatCapacity.T_HeatCapacity_Moitra_smoothed),
the default `Cp`. The default ``\\Delta T_g`` is 20 K.
"""
@with_kw_noshow struct T_Conductivity_Moitra_smoothed{T, U1, U2, U3, U4, S1 <: AbstractDensity, S2 <: AbstractHeatCapacity} <: AbstractConductivity{T}
    ρ::S1 = ConstantDensity(ρ = 3000.0kg / m^3)
    Cp::S2 = T_HeatCapacity_Moitra_smoothed()
    a0::GeoUnit{T, U1} = 5.98e-6m^2 / s * K^(41 // 100)   # glass
    b0::GeoUnit{T, U2} = 9.23e-11m^2 / s / K
    a1::GeoUnit{T, U3} = -6.0e-8m^2 / s                    # melt
    b1::GeoUnit{T, U2} = 4.57e-10m^2 / s / K
    Tg::GeoUnit{T, U4} = 1000.0K                           # glass transition temperature
    ΔTg::GeoUnit{T, U4} = 20.0K                            # width of the transition
end
T_Conductivity_Moitra_smoothed(args...) = promote_construct(T_Conductivity_Moitra_smoothed, args...)

function param_info(s::T_Conductivity_Moitra_smoothed)
    return MaterialParamsInfo(; Equation = L"k = \rho c_p \alpha;\ \alpha = (1 - H) \alpha_{glass} + H \alpha_{melt};\ H = 1/(1 + \exp(-4 (T - T_g)/\Delta T_g))")
end

function (s::T_Conductivity_Moitra_smoothed)(; T, kwargs...)
    Tc = precision_of(T)
    @unpack_like T Tc a0, b0, a1, b1, Tg, ΔTg = s
    H = inv(1 + exp(-4 * (T - Tg) / ΔTg))
    α = (1 - H) * (a0 * T^(-41 // 100) + b0 * T) + H * (a1 + b1 * T)
    return compute_density(s.ρ, (; T, kwargs...)) * compute_heatcapacity(s.Cp, (; T, kwargs...)) * α
end

function show(io::IO, g::T_Conductivity_Moitra_smoothed)
    return print(io, "Conductivity from the glass/melt diffusivity of Moitra et al. (2018), smoothed: ΔTg=$(UnitValue(g.ΔTg)); ρ=$(g.ρ); Cp=$(g.Cp)")
end

# Conductivity of H2O, IAPWS 2011 (Huber et al. 2012) without critical enhancement ----
"""
    IAPWS_H2O_Conductivity(; ρ=HalbachChatterjee_Density(), L0, L1)

Thermal conductivity of H2O from the IAPWS 2011 formulation (Huber et al. 2012; IAPWS R15-11, eqs. 15-17) without
the critical-enhancement term ``\\bar\\lambda_2``:
```math
    \\lambda = \\lambda^* \\bar\\lambda_0(\\bar T)\\, \\bar\\lambda_1(\\bar T, \\bar\\rho), \\qquad
    \\bar\\lambda_0 = \\frac{\\sqrt{\\bar T}}{\\sum_{k=0}^{4} L_k / \\bar T^k}, \\qquad
    \\bar\\lambda_1 = \\exp\\left[\\bar\\rho \\sum_{i=0}^{4} \\left(\\frac{1}{\\bar T} - 1\\right)^i
        \\sum_{j=0}^{5} L_{ij} (\\bar\\rho - 1)^j\\right]
```
with ``\\bar T = T/T^*``, ``\\bar\\rho = \\rho/\\rho^*``, ``T^* = 647.096`` K, ``\\rho^* = 322`` kg/m³ and
``\\lambda^* = 10^{-3}`` W/m/K (eqs. 1, 3, 4). The dimensionless coefficients `L0` (``L_k``, Table 1) and `L1`
(``L_{ij}``, Table 2, row ``i``, column ``j``) default to the values of IAPWS R15-11.
The density is the law `ρ` evaluated with the same arguments (`P`, `T`).

IAPWS R15-11 validates the formulation (its eq. 14) up to 1173.15 K for ``P \\le`` 100 MPa, up to 874 K for
100-250 MPa, and up to 573 K above 250 MPa; most of the range 600-1200 °C, 0.1-400 MPa is therefore an extrapolation
of the standard. Within that range, omitting ``\\bar\\lambda_2`` changes ``\\lambda`` by at most 1.99 % (600 °C,
80 MPa) and not at all above about 697 °C, and with the default density ``\\lambda`` is within 1.5 % of the full
(extrapolated) IAPWS formulation. The law is not valid near the critical point (374 °C, 22 MPa), where
``\\bar\\lambda_2`` diverges. The reference quantities `Tref`, `ρref` and `kref`
are ``T^*``, ``\\rho^*`` and ``\\lambda^*``; they are not material parameters.

# References
- Huber, M. L., Perkins, R. A., Friend, D. G., et al. (2012), New International Formulation for the Thermal Conductivity of H2O, J. Phys. Chem. Ref. Data 41, 033102, https://doi.org/10.1063/1.4738955
- IAPWS R15-11 (2011), Release on the IAPWS Formulation 2011 for the Thermal Conductivity of Ordinary Water Substance, https://iapws.org/technical-guidance/release/ThCond
"""
@with_kw_noshow struct IAPWS_H2O_Conductivity{T, U1, U2, U3, S <: AbstractDensity} <: AbstractConductivity{T}
    ρ::S = HalbachChatterjee_Density()
    Tref::GeoUnit{T, U1} = 647.096K
    ρref::GeoUnit{T, U2} = 322.0kg / m^3
    kref::GeoUnit{T, U3} = 1.0e-3Watt / m / K
    L0::NTuple{5, T} = (2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4)
    L1::NTuple{5, NTuple{6, T}} = (
        (1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258),
        (2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245),
        (2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816),
        (-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0.0, 0.0),
        (-2.720337, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842),
    )
end
IAPWS_H2O_Conductivity(args...) = promote_construct(IAPWS_H2O_Conductivity, args...)

function param_info(s::IAPWS_H2O_Conductivity)
    return MaterialParamsInfo(; Equation = L"\lambda = \lambda^* \bar\lambda_0(\bar T) \bar\lambda_1(\bar T, \bar\rho)")
end

function (s::IAPWS_H2O_Conductivity)(; P = 0, T = 0, kwargs...)
    Tc = precision_of((; P, T))
    P, T = convert_precision(Tc, P), convert_precision(Tc, T)
    @unpack_like P Tc Tref, ρref, kref = s
    t = NoUnits(T / Tref)
    d = NoUnits(compute_density(s.ρ, (; P, T, kwargs...)) / ρref)
    λ0 = sqrt(t) / evalpoly(1 / t, convert_precision(Tc, s.L0))
    λ1 = exp(d * evalpoly(1 / t - 1, map(L -> evalpoly(d - 1, L), convert_precision(Tc, s.L1))))
    return kref * λ0 * λ1
end

function show(io::IO, g::IAPWS_H2O_Conductivity)
    return print(io, "H2O conductivity, IAPWS 2011 without critical enhancement: ρ=$(g.ρ)")
end

# Conductivity of an H2O-CO2 gas mixture (Udoetok 2013) ----------------------
"""
    GasMixture_Conductivity(; kh2o=IAPWS_H2O_Conductivity(), kco2=ConstantConductivity(k=0.071W/m/K))

Thermal conductivity of a binary H2O-CO2 gas mixture after Udoetok (2013), eq. 16, the mean of the parallel (eq. 14)
and series (eq. 9) arrangements of the two gases:
```math
    k = \\frac{1}{2}\\left(X_{H_2O} k_{H_2O} + X_{CO_2} k_{CO_2}
        + \\frac{k_{H_2O} k_{CO_2}}{X_{H_2O} k_{CO_2} + X_{CO_2} k_{H_2O}}\\right),
    \\qquad X_{H_2O} = 1 - X_{CO_2}
```
The CO2 mole fraction of the gas `X_co2` is passed at runtime, as for the H2O-CO2 solubility laws; `kh2o` and
`kco2` are conductivity laws evaluated with the remaining arguments (e.g. `T`). The result is `NaN` outside
``0 \\le X_{CO_2} \\le 1``. The default `kh2o` depends on `P` and `T`; the default `kco2` = 0.071 W/m/K is the
conductivity of CO2 at 1000 K and low pressure. Over 600-1200 °C and 0.1-400 MPa the conductivity of CO2 ranges from
0.06 to 0.15 W/m/K, so the default underestimates it by up to a factor of about 2 at high pressure.

# References
- Udoetok, E. S. (2013), Thermal conductivity of binary mixtures of gases, Frontiers in Heat and Mass Transfer 4, 023008, https://doi.org/10.5098/hmt.v4.2.3008
"""
struct GasMixture_Conductivity{T, S1 <: AbstractConductivity{T}, S2 <: AbstractConductivity} <: AbstractConductivity{T}
    kh2o::S1
    kco2::S2
end
function GasMixture_Conductivity(; kh2o = IAPWS_H2O_Conductivity(), kco2 = ConstantConductivity(k = 0.071Watt / m / K))
    return GasMixture_Conductivity(kh2o, kco2)
end

function param_info(s::GasMixture_Conductivity)
    return MaterialParamsInfo(; Equation = L"k = (X_w k_w + X_c k_c + k_w k_c/(X_w k_c + X_c k_w))/2")
end

function (s::GasMixture_Conductivity)(; X_co2 = 0, kwargs...)
    kw, kc = compute_conductivity(s.kh2o, kwargs), compute_conductivity(s.kco2, kwargs)
    0 <= X_co2 <= 1 || return kw * precision_of(kw)(NaN)
    Xw = 1 - X_co2
    return (Xw * kw + X_co2 * kc + kw * kc / (Xw * kc + X_co2 * kw)) / 2
end

function show(io::IO, g::GasMixture_Conductivity)
    return print(io, "H2O-CO2 gas-mixture conductivity following Udoetok (2013): kh2o=$(g.kh2o); kco2=$(g.kco2)")
end

# Conductivity of a magma with spherical gas bubbles (Maxwell; Heap et al. 2020) ----
"""
    Maxwell_Conductivity(; kmatrix=T_Conductivity_Moitra(), kgas=GasMixture_Conductivity())

Effective thermal conductivity of a matrix containing non-interacting spherical gas bubbles, the Maxwell equation
used by Heap et al. (2020), eq. 1, with ``\\beta = 3(1-m)/(2+m)`` for spherical pores:
```math
    \\frac{k}{k_0} = \\frac{(1-\\phi_{gas})(1-m) + m\\beta\\phi_{gas}}{(1-\\phi_{gas})(1-m) + \\beta\\phi_{gas}}
                 = \\frac{(1-\\phi_{gas})(2k_0 + k_{gas}) + 3\\phi_{gas} k_{gas}}{(1-\\phi_{gas})(2k_0 + k_{gas}) + 3\\phi_{gas} k_0},
    \\qquad m = \\frac{k_{gas}}{k_0}
```
where ``k_0`` is the conductivity of the bubble-free magma (`kmatrix`: melt, including any crystals) and ``k_{gas}``
that of the gas (`kgas`). The second form, used here, stays finite at ``m = 1``. The gas volume fraction of the
whole magma `ϕ_gas` is passed at runtime; both sub-laws receive the remaining arguments (e.g. `P`, `T`, `X_co2`; the default `kgas` needs `P`
and `T`). `kgas` is only evaluated when `ϕ_gas != 0`. The result is `NaN` outside ``0 \\le \\phi_{gas} \\le 1``.

The Maxwell equation assumes dilute, non-interacting bubbles. For ``k_{gas} < k_0`` it is the Hashin–Shtrikman
upper bound, so it overestimates ``k`` of foams with connected or deformed bubbles; Heap et al. (2020) test it up to
``\\phi_{gas} \\approx 0.6``. The default `kmatrix` inherits the melt-branch limitation of
[`T_Conductivity_Moitra`](@ref GeoParams.MaterialParameters.Conductivity.T_Conductivity_Moitra) above ``T_g``.

# References
- Heap, M. J., et al. (2020), The thermal properties of porous andesite, JVGR 398, 106901, https://doi.org/10.1016/j.jvolgeores.2020.106901
"""
struct Maxwell_Conductivity{T, S1 <: AbstractConductivity{T}, S2 <: AbstractConductivity} <: AbstractConductivity{T}
    kmatrix::S1
    kgas::S2
end
Maxwell_Conductivity(; kmatrix = T_Conductivity_Moitra(), kgas = GasMixture_Conductivity()) = Maxwell_Conductivity(kmatrix, kgas)

function param_info(s::Maxwell_Conductivity)
    return MaterialParamsInfo(; Equation = L"k = k_0 \frac{(1-\phi)(2k_0 + k_g) + 3\phi k_g}{(1-\phi)(2k_0 + k_g) + 3\phi k_0}")
end

function (s::Maxwell_Conductivity)(; ϕ_gas = 0, kwargs...)
    k0 = compute_conductivity(s.kmatrix, kwargs)
    0 <= ϕ_gas <= 1 || return k0 * precision_of(k0)(NaN)
    iszero(ϕ_gas) && return k0
    kg = compute_conductivity(s.kgas, kwargs)
    a = (1 - ϕ_gas) * (2 * k0 + kg)
    return k0 * (a + 3 * ϕ_gas * kg) / (a + 3 * ϕ_gas * k0)
end

function show(io::IO, g::Maxwell_Conductivity)
    return print(io, "Maxwell conductivity of a bubbly magma: kmatrix=$(g.kmatrix); kgas=$(g.kgas)")
end

# add methods programmatically
for myType in (
        :T_Conductivity_Whittington_parameterised,
        :T_Conductivity_Whittington,
        :ConstantConductivity,
        :TP_Conductivity,
        :T_Conductivity_Moitra,
        :T_Conductivity_Moitra_smoothed,
        :IAPWS_H2O_Conductivity,
        :GasMixture_Conductivity,
        :Maxwell_Conductivity,
    )
    @eval begin
        (s::$(myType))(args) = s(; args...)
        compute_conductivity(s::$(myType), args) = s(args)
    end
end

# Print info
function show(io::IO, g::TP_Conductivity)
    return if ustrip(Value(g.d)) == 0
        print(
            io,
            "T/P dependent conductivity: Name = $(String(collect(g.Name))), k = $(g.a.val) + $(g.b.val)/(T + $(g.c.val))  \n",
        )
    else
        print(
            io,
            "T/P dependent conductivity: Name = $(String(collect(g.Name))), k = ($(g.a.val) + $(g.b.val)/(T + $(g.c.val)))*(1 + $(g.d.val)*P)  \n",
        )
    end
end
#-------------------------------------------------------------------------

# Help info for the calculation routines
"""
    k = compute_conductivity(s::AbstractConductivity; P, T)
    k = compute_conductivity(s::AbstractConductivity, args)

Returns the thermal conductivity `k` for the parameterization `s`, evaluated at the temperature `T`
and pressure `P` passed either as keyword arguments or as a `NamedTuple` `args`.

Currently available:
- ConstantConductivity
- T\\_Conductivity_Whittington
- TP\\_Conductivity
- T\\_Conductivity_Moitra
- T\\_Conductivity_Moitra\\_smoothed
- IAPWS\\_H2O_Conductivity
- GasMixture\\_Conductivity
- Maxwell\\_Conductivity

# Example
```jldoctest
julia> k = T_Conductivity_Whittington();

julia> compute_conductivity(k; T = 1000.0)
1.9014576517265882
```


"""
compute_conductivity()

"""
    k = compute_conductivity(T::Any, s::AbstractConductivity)

Computes conductivity if only temperature (and not pressure) is specified
"""
# compute_conductivity(s::AbstractConductivity, T::AbstractArray{_T}) where _T =  compute_conductivity(s,similar(T),T)

"""
    k = compute_conductivity(s::ConstantConductivity)

Returns conductivity if we are sure that we will only employ constant values throughout the simulation
"""
#compute_conductivity(s::ConstantConductivity) =  compute_conductivity(s,0,0)

# Computational routines needed for computations with the MaterialParams structure
function compute_conductivity(s::AbstractMaterialParamsStruct, args)
    return s.Conductivity[1](args)
end

compute_conductivity(MatParam, arg, args::Vararg{Any, N}) where {N} = compute_param(compute_conductivity, MatParam, arg, args...)

"""
    compute_conductivity!(K::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},Temp::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of conductivity `K` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `Temp` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.

_________________________________________________________________________________________________________

    compute_conductivity!(k::AbstractArray{T,N}, PhaseRatios::AbstractArray{T, M}, P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of conductivity `k` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)
"""
compute_conductivity!(args::Vararg{Any, N}) where {N} = compute_param!(compute_conductivity, args...)

# extractor methods
for type in (ConstantConductivity, T_Conductivity_Whittington, T_Conductivity_Whittington_parameterised, TP_Conductivity, T_Conductivity_Moitra, T_Conductivity_Moitra_smoothed, IAPWS_H2O_Conductivity, GasMixture_Conductivity, Maxwell_Conductivity)
    @extractors(type, :Conductivity)
end

end
