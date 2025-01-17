module Conductivity

# This implements different methods to specify conductivity of rocks
#
# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, @extractors, add_extractor_functions
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractConductivity{T} <: AbstractMaterialParam end

export compute_conductivity,                       # calculation routines
    compute_conductivity!,
    param_info,
    ConstantConductivity,                       # constant
    T_Conductivity_Whittington,                 # T-dependent heat capacity
    T_Conductivity_Whittington_parameterised,   # parameterised T-dependent heat capacity
    TP_Conductivity,                            # TP dependent conductivity
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
@with_kw_noshow struct ConstantConductivity{T,U} <: AbstractConductivity{T}
    k::GeoUnit{T,U} = 3.0Watt / m / K
end
ConstantConductivity(args...) = ConstantConductivity(convert.(GeoUnit, args)...)

function param_info(s::ConstantConductivity) # info about the struct
    return MaterialParamsInfo(; Equation=L"k = cst")
end

# Calculation routine
function (s::ConstantConductivity)(; kwargs...)
    @unpack_val k = s

    return k
end

compute_conductivity(s::ConstantConductivity; kwargs...) = s()

function (s::ConstantConductivity)(I::Integer...)
    @unpack_val k = s

    return fill(k, I...)
end

compute_conductivity(s::ConstantConductivity, I::Integer...) = s(I...)

"""
    compute_conductivity(k_array::AbstractArray{<:AbstractFloat,N},P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, s::ConstantConductivity) where N

In-place routine to compute constant conductivity
"""
function compute_conductivity!(
    k_array::AbstractArray{_T,N}, s::ConstantConductivity; kwargs...
) where {_T, N} end

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
- d = 576.3m^2/s*K
- e = 0.062m^2/s
- f = 0.732m^2/s
- g = 0.000135m^2/s/K

This looks like:

![subet1](./assets/img/Conductivity_Whittington.png)

Example
===
```julia
julia> using GLMakie, GeoParams
julia> p=T_Conductivity_Whittington();
julia> T,k,plt = PlotConductivity(p)
```


"""
@with_kw_noshow struct T_Conductivity_Whittington{T,U1,U2,U3,U4,U5,U6,U7,U8,U9} <:
                       AbstractConductivity{T}
    # Note: the resulting curve of k was visually compared with Fig. 2 of the paper
    a0::GeoUnit{T,U1} = 199.5J / mol / K                # prefactor for low T       (T<= 846 K)
    a1::GeoUnit{T,U1} = 229.32J / mol / K               # prefactor for high T      (T>  846 K)
    b0::GeoUnit{T,U2} = 0.0857J / mol / K^2             # linear term for low T     (T<= 846 K)
    b1::GeoUnit{T,U2} = 0.0323J / mol / K^2             # linear term for high T    (T>  846 K)
    c0::GeoUnit{T,U3} = 5e6J / mol * K                  # quadratic term for low T  (T<= 846 K)
    c1::GeoUnit{T,U3} = 47.9e-6J / mol * K              # quadratic term for high T (T>  846 K)
    molmass::GeoUnit{T,U4} = 0.22178kg / mol               # average molar mass
    Tcutoff::GeoUnit{T,U5} = 846.0K                      # cutoff temperature
    rho::GeoUnit{T,U6} = 2700kg / m^3                  # Density they use for an average crust
    d::GeoUnit{T,U7} = 576.3 * 1e-6m^2 / s * K           # diffusivity parameterization
    e::GeoUnit{T,U8} = 0.062 * 1e-6m^2 / s             # diffusivity parameterization
    f::GeoUnit{T,U8} = 0.732 * 1e-6m^2 / s             # diffusivity parameterization
    g::GeoUnit{T,U9} = 0.000135 * 1e-6m^2 / s / K        # diffusivity parameterization
end
T_Conductivity_Whittington(args...) = T_Conductivity_Whittington(convert.(GeoUnit, args)...)

function param_info(s::T_Conductivity_Whittington) # info about the structwhere {_T,N}
    return MaterialParamsInfo(; Equation=L"k = f(T) ")
end

# Calculation routine
function (s::T_Conductivity_Whittington{_T})(; T=0e0, kwargs...) where {_T}
    if T isa Quantity
        @unpack_units a0, a1, b0, b1, c0, c1, molmass, Tcutoff, rho, d, e, f, g = s
    else
        @unpack_val a0, a1, b0, b1, c0, c1, molmass, Tcutoff, rho, d, e, f, g = s
    end

    if T ≤ Tcutoff
        return (a0 + b0 * T - c0 / T^2) / molmass * (d / T - e) * rho
    else
        return (a1 + b1 * T - c1 / T^2) / molmass * (f - g * T) * rho
    end
end

function compute_conductivity(s::T_Conductivity_Whittington{_T}; T=0e0) where {_T}
    return s(; T=T)
end

function (s::T_Conductivity_Whittington)(T::AbstractArray; kwargs...)
    if T isa Quantity
        @unpack_units a0, a1, b0, b1, c0, c1, molmass, Tcutoff, rho, d, e, f, g = s
    else
        @unpack_val a0, a1, b0, b1, c0, c1, molmass, Tcutoff, rho, d, e, f, g = s
    end

    k = similar(T)   #creating an array makes 1 allocation
    inv_molmass = 1 / molmass # multiplication is considerably faster than division
    @inbounds for i in eachindex(T)
        if T[i] <= Tcutoff
            a, b, c = a0, b0, c0
            κ = d / T[i] - e
        else
            a, b, c = a1, b1, c1
            κ = f - g * T[i]
        end
        #  κ = κ*1e-6

        Cp = (a + b * T[i] - c / T[i]^2) * inv_molmass # conductivity
        k[i] = κ * rho * Cp       # compute conductivity from diffusivity
    end

    return k
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
@with_kw_noshow struct T_Conductivity_Whittington_parameterised{T,U1,U2,U3,U4,U5} <:
                       AbstractConductivity{T}
    # Note: the resulting curve of k was visually compared with Fig. 2 of the paper
    a::GeoUnit{T,U1} = -2e-09Watt / m / K^4                #
    b::GeoUnit{T,U2} = 6e-06Watt / m / K^3              #
    c::GeoUnit{T,U3} = -0.0062Watt / m / K^2              #
    d::GeoUnit{T,U4} = 4Watt / m / K                #
    Ts::GeoUnit{T,U5} = 273.15 * K                    # Temperature shift [C->K]
end
function T_Conductivity_Whittington_parameterised(args...)
    return T_Conductivity_Whittington_parameterised(convert.(GeoUnit, args)...)
end

function param_info(s::T_Conductivity_Whittington_parameterised) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"k = a*(T[K]-Ts)^3 + b*(T[K]-Ts)^2  + c*(T[K]-Ts) + d, Ts=273.15"
    )
end

# Calculation routine
function (s::T_Conductivity_Whittington_parameterised{_T})(;
    T=0e0, kwargs...
) where {_T}
    if T isa Quantity
        @unpack_units a, b, c, d, Ts = s
    else
        @unpack_val a, b, c, d, Ts = s
    end

    T_C = T - Ts
    return a * T_C^3 + b * T_C^2 + c * T_C + d
end

function (s::T_Conductivity_Whittington_parameterised)(T::AbstractArray; kwargs...)
    if T isa Quantity
        @unpack_units a, b, c, d, Ts = s
    else
        @unpack_val a, b, c, d, Ts = s
    end

    k = similar(T)   # creating an array makes 1 allocation

    @inbounds for i in eachindex(T)
        # Note: in general, we operate with SI units or the non-dimensional equivalent of that
        # This parameterisation was developed
        T_C = T[i] - Ts
        k[i] = a * T_C^3 + b * T_C^2 + c * T_C + d
    end

    return k
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
    k = \\left( a_k +  {b_k \\over {T + c_k}} \\right) (1 + d_k P)
```

where ``k`` is the conductivity [``W/K/m``], and ``a_k,b_k,c_k,d_k`` are parameters that dependent on the temperature `T` and pressure `P`:
- ``a_k`` = 1.18Watt/K/m
- ``b_k`` = 474Watt/m
- ``c_k`` = 77K
- ``d_k`` = 0/MPa
"""
@with_kw_noshow struct TP_Conductivity{T,N,U1,U2,U3,U4} <: AbstractConductivity{T}
    Name::NTuple{N,Char} = ""                  # The name is encoded as a NTuple{Char} to make it isbits
    a::GeoUnit{T,U1} = 1.18Watt / K / m        # empirical fitting term
    b::GeoUnit{T,U2} = 474.0Watt / m         # empirical fitting term
    c::GeoUnit{T,U3} = 77.0K               # empirical fitting term
    d::GeoUnit{T,U4} = 0.0 / MPa             # empirical fitting term
end
function TP_Conductivity(args...)
    return TP_Conductivity(
        NTuple{length(args[1]),Char}(collect.(args[1])), convert.(GeoUnit, args[2:end])...
    )
end

function param_info(s::TP_Conductivity)
    name = String(collect(s.Name))
    eq = L"k = \left(a_k + {b_k/{T + c_k}} \right)*(1 + d_k*P) "
    if name == ""
        return MaterialParamsInfo(; Equation=eq)
    end
    return MaterialParamsInfo(; Equation=eq, Comment=TP_Conductivity_info[name][2].Comment)
end

"""
    Set_TP_Conductivity["Name of temperature(-pressure) dependent conductivity"]

This is a dictionary with pre-defined laws:
- "UpperCrust"
- "LowerCrust"
- "OceanicCrust"
- "Mantle"

# Example
```julia
julia> k=Set_TP_Conductivity["Mantle"]
T/P dependent conductivity: k = (0.73 W K⁻¹ m⁻¹ + 1293 W m⁻¹/(T + 77 K))*(1 + 4.0e-5 MPa⁻¹*P)
```

"""
Set_TP_Conductivity(name::String) = TP_Conductivity_info[name][1]

TP_Conductivity_info = Dict(
    [
        (
            "UpperCrust",
            (
                TP_Conductivity(;
                    Name="UpperCrust", a=0.64Watt / K / m, b=807Watt / m, c=77K, d=0 / MPa
                ),
                MaterialParamsInfo(;
                    Comment="Sediment/upper crust T-dependent conductivity, as listed in table 21.2 of Gerya et al. | Reference still to be verified!",
                ),
            ),
        )
        (
            "LowerCrust",
            (
                TP_Conductivity(;
                    Name="LowerCrust", a=1.18Watt / K / m, b=474Watt / m, c=77K, d=0 / MPa
                ),
                MaterialParamsInfo(;
                    Comment="Lower crust T-dependent conductivity, as listed in table 21.2 of Gerya et al. | Reference still to be verified!",
                ),
            ),
        )
        (
            "OceanicCrust",
            (
                TP_Conductivity(;
                    Name="OceanicCrust", a=1.18Watt / K / m, b=474Watt / m, c=77K, d=0 / MPa
                ),
                MaterialParamsInfo(;
                    Comment="Oceanic crust T-dependent conductivity, as listed in table 21.2 of Gerya et al. | Reference still to be verified!",
                ),
            ),
        )
        (
            "Mantle",
            (
                TP_Conductivity(;
                    Name="Mantle",
                    a=0.73Watt / K / m,
                    b=1293Watt / m,
                    c=77K,
                    d=0.00004 / MPa,
                ),
                MaterialParamsInfo(;
                    Comment="Mantle T-dependent conductivity, as listed in table 21.2 of Gerya et al. | Reference still to be verified!",
                ),
            ),
        )
    ],
)

# Calculation routine
function (s::TP_Conductivity{_T})(; P=0e0, T=0e0, kwargs...) where {_T}
    if T isa Quantity
        @unpack_units a, b, c, d = s
    else
        @unpack_val a, b, c, d = s
    end

    if ustrip(d) == 0
        return a + b / (T + c)
    else
        return (a + b / (T + c)) * (1 + d * P)
    end
end


function (s::TP_Conductivity{_T})(
    P::AbstractArray, T::AbstractArray; kwargs...
) where {_T}
    if T isa Quantity
        @unpack_units a, b, c, d = s
    else
        @unpack_val a, b, c, d = s
    end

    k = similar(T)

    @inbounds if ustrip(d) == 0
        for i in eachindex(T)
            k[i] = a + b / (T[i] + c)
        end
    else
        if size(T) != size(P)
            error("Size of P and T arrays should be the same")
        end

        for i in eachindex(T)
            k[i] = (a + b / (T[i] + c)) * (1 + d * P[i])
        end
    end

    return k
end


#(s::TP_Conductivity)(P::AbstractArray, T::AbstractArray) = s(P, T)
compute_conductivity(s::TP_Conductivity, P::AbstractArray, T::AbstractArray) = s(P, T)

# add methods programmatically
for myType in (
    :T_Conductivity_Whittington_parameterised,
    :T_Conductivity_Whittington,
    :ConstantConductivity,
    :TP_Conductivity,
)
    @eval begin
        (s::$(myType))(args) = s(; args...)
        compute_conductivity(s::$(myType), args) = s(args)
    end
end

# Print info
function show(io::IO, g::TP_Conductivity)
    if ustrip(Value(g.d)) == 0
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
    k = compute_conductivity(P, T, s:<AbstractConductivity)

Returns the thermal conductivity `k` at any temperature `T` and pressure `P` using any of the parameterizations implemented.

Currently available:
- ConstantConductivity
- T\\_Conductivity_Whittington
- TP\\_Conductivity

# Example
Using dimensional units
```julia
julia> T  = (250:100:1250)*K;
julia> cp = T_HeatCapacity_Whittington()
julia> Cp = ComputeHeatCapacity(0,T,cp)
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

"""
    compute_conductivity!(K::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},Temp::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of conductivity `K` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `Temp` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.

_________________________________________________________________________________________________________

compute_conductivity!(k::AbstractArray{T,N}, PhaseRatios::AbstractArray{T, M}, P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of conductivity `k` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)
"""
compute_conductivity(args::Vararg{Any, N}) where N = compute_param(compute_conductivity, args...)
compute_conductivity!(args::Vararg{Any, N}) where N = compute_param!(compute_conductivity, args...)

# extractor methods
for type in (ConstantConductivity, T_Conductivity_Whittington, T_Conductivity_Whittington_parameterised, TP_Conductivity)
    @extractors(type, :Conductivity)
end

end
