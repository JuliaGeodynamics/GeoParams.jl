module MeltingParam

# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful, ForwardDiff
using ..Units
using GeoParams:
    AbstractMaterialParam, PhaseDiagram_LookupTable, AbstractMaterialParamsStruct
import Base.show, GeoParams.param_info 
using ..MaterialParameters: MaterialParamsInfo
using Setfield # allows modifying fields in immutable struct

abstract type AbstractMeltingParam{T} <: AbstractMaterialParam end

export compute_meltfraction,
    compute_meltfraction!,    # calculation routines
    compute_meltfraction_ratio,
    compute_dϕdT,             # derivative of melt fraction versus
    compute_dϕdT!,
    param_info,
    MeltingParam_Caricchi,
    MeltingParam_Smooth3rdOrder,
    MeltingParam_4thOrder,
    MeltingParam_5thOrder,
    MeltingParam_Quadratic,
    MeltingParam_Assimilation,
    SmoothMelting

include("../Utils.jl")
include("../Computations.jl")

# Caricchi  -------------------------------------------------------
"""
    MeltingParam_Caricchi()

Implements the T-dependent melting parameterisation used by Caricchi, Simpson et al. (as for example described in Simpson)
```math
    \\theta = {(a - (T + c)) \\over b}
```
```math
    \\phi_{melt} = {1.0 \\over (1.0 + e^\\theta)}
```

Note that T is in Kelvin. As default parameters we employ:
```math
b=23K,  a=800K , c=273.15K
```
Which gives a reasonable fit to experimental data of granodioritic composition (Piwinskii and Wyllie, 1968):

![MeltParam_Carrichi](./assets/img/MeltingParam_Caricchi.png)

References
====
- Simpson G. (2017) Practical finite element modelling in Earth Sciences Using MATLAB.
"""
@with_kw_noshow struct MeltingParam_Caricchi{T,U} <: AbstractMeltingParam{T}
    a::GeoUnit{T,U} = 800.0K
    b::GeoUnit{T,U} = 23.0K
    c::GeoUnit{T,U} = 273.15K # shift from C to K
    apply_bounds::Bool = true
end
MeltingParam_Caricchi(args...) = MeltingParam_Caricchi(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_Caricchi) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"\phi = {1 \over {1 + \exp( {800-T[^oC] \over 23})}}"
    )
end

# Calculation routines
function (p::MeltingParam_Caricchi)(; T, kwargs...)
    @unpack_val a, b, c = p
    θ = (a - (T - c)) / b
    return 1.0 / (1.0 + exp(θ))
end

function compute_dϕdT(p::MeltingParam_Caricchi; T, kwargs...)
    @unpack_val a, b, c = p

    _b = inv(b)
    dϕdT = exp((a + c - T) * _b) / (b * ((1.0 + exp((a + c - T) * _b))^2))

    return dϕdT
end

# Print info
function show(io::IO, g::MeltingParam_Caricchi)
    return print(io, "Caricchi et al. melting parameterization")
end
#-------------------------------------------------------------------------

# Melnik  -------------------------------------------------------
"""
    MeltingParam_Smooth3rdOrder()

Implements the a smooth 3rd order T-dependent melting parameterisation (as used by Melnik and coworkers)
```math
    x = {  (T - 273.15) \\over 1000.0}
```
```math
    \\theta = { a + b * x + c * x^2 + d * x^3}
```
```math
    \\phi_{melt} = {1.0 \\over (1.0 + e^\\theta)}
```
Note that T is in Kelvin. 

As default parameters we employ:
```math
a=517.9,  b=-1619.0, c=1699.0, d = -597.4
```
which gives a reasonable fit to experimental data for basalt.

Data for rhyolite are:
```math
a=3043.0,  b=-10552.0, c=12204.9, d = -4709.0 
```

![MeltingParam_Smooth3rdOrder](./assets/img/MeltingParam_Melnik.png)
Red: Rhyolite, Blue: Basalt

References
====
"""
@with_kw_noshow struct MeltingParam_Smooth3rdOrder{T,U,U1} <: AbstractMeltingParam{T}
    a::GeoUnit{T,U} =  517.9NoUnits
    b::GeoUnit{T,U} = -1619.0NoUnits
    c::GeoUnit{T,U} = 1699.0NoUnits
    d::GeoUnit{T,U} = -597.4NoUnits
    T0::GeoUnit{T,U1} = 273.15K 
    Tchar::GeoUnit{T,U1} = 1000K # normalization
    apply_bounds::Bool = true
end
MeltingParam_Smooth3rdOrder(args...) = MeltingParam_Smooth3rdOrder(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_Smooth3rdOrder) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"\phi = f(T)"
    )
end

# Calculation routines
function (p::MeltingParam_Smooth3rdOrder)(; T, kwargs...)
    @unpack_val a, b, c, d, Tchar, T0 = p
    x = (T - T0)/ Tchar
    
    θ = min(evalpoly(x, (a, b, c, d)),200.0)

    ϕ = inv(1.0 + exp(θ))
    return ϕ
end

function compute_dϕdT(p::MeltingParam_Smooth3rdOrder; T, kwargs...)
    @unpack_val a, b, c, d, Tchar, T0 = p
    
    x =  (T - T0) / Tchar
    θ = min(evalpoly(x, (a, b, c, d)),200.0)

    dϕdT = -exp(θ)*(1 / (1.0 + exp(θ))^2 )*(b / Tchar + (3d*x^2) / Tchar + (2(T - T0)*c) / (Tchar^2))

    return dϕdT
end

function foo(T,T0,Tchar,a,b,c,d)
    x =  (T - T0) / Tchar
    θ = a + b * x + c * x^2 + d * x^3;
    ϕ = inv(1.0 + exp(θ))
    return x^2
end

# Print info
function show(io::IO, g::MeltingParam_Smooth3rdOrder)
    return print(io, "Smooth 3rd order melting parameterization")
end
#-------------------------------------------------------------------------

# MeltingParam_5thOrder  -------------------------------------------------------
"""
    MeltingParam_5thOrder(a,b,c,d,e,f,T_s,T_l)

Uses a 5th order polynomial to describe the melt fraction `phi` between solidus temperature `T_s` and liquidus temperature `T_l`
```math
    \\phi = a T^5 + b T^4 + c T^3 + d T^2 + e T + f  \\textrm{   for   } T_s ≤ T ≤ T_l
```
```math
    \\phi = 1  \\textrm{   if   } T>T_l
```
```math
    \\phi = 0  \\textrm{   if   } T<T_s
```
Temperature `T` is in Kelvin.

![MeltingParam_5thOrder](./assets/img/MeltingParam_5thorder.png)

The default values are for a composite liquid-line-of-descent:
- the upper part is for Andesite from: (Blatter, D. L. & Carmichael, I. S. (2001) Hydrous phase equilibria of a Mexican highsilica andesite: a candidate for a mantle origin? Geochim. Cosmochim. Acta 65, 4043–4065
- the lower part is extrapolated to the granitic minimum using the Marxer & Ulmer LLD for Andesite (Marxer, F. & Ulmer, P. (2019) Crystallisation and zircon saturation of calc-alkaline tonalite from the Adamello Batholith at upper crustal conditions: an experimental study. *Contributions Mineral. Petrol.* 174, 84)

"""
@with_kw_noshow struct MeltingParam_5thOrder{T,U,U1,U2,U3,U4,U5,U6} <:
                       AbstractMeltingParam{T}
    a::GeoUnit{T,U1} = 2.083291971482524e-12 / K^5
    b::GeoUnit{T,U2} = -1.239502833666574e-08 / K^4
    c::GeoUnit{T,U3} = 2.938887604687626e-05 / K^3
    d::GeoUnit{T,U4} = -0.034711533077108 / K^2
    e::GeoUnit{T,U5} = 20.425403874539178 / K
    f::GeoUnit{T,U6} = -4.790664658179178e+03 * NoUnits
    T_s::GeoUnit{T,U} = 963.15K
    T_l::GeoUnit{T,U} = 1388.2K
    apply_bounds::Bool = true
end
MeltingParam_5thOrder(args...) = MeltingParam_5thOrder(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_5thOrder) # info about the struct
    return MaterialParamsInfo(; Equation=L"\phi = aT^5 + bT^4 + cT^3 + dT^2 + eT + f")
end


# Calculation routines
function (p::MeltingParam_5thOrder)(; T, kwargs...)
    @unpack_val a, b, c, d, e, f, T_s, T_l = p

    coeffs = f, e, d, c, b, a
    ϕ = evalpoly(T, coeffs)
    if p.apply_bounds
        if T < T_s
            ϕ = 0.0
        elseif T > T_l
            ϕ = 1.0
        end
    end

    return ϕ
end

compute_dϕdT(p::MeltingParam_5thOrder, T, kwargs...) = compute_dϕdT(p; T, kwargs...) 

function compute_dϕdT(p::MeltingParam_5thOrder; T, kwargs...)
    @unpack_val a, b, c, d, e, T_s, T_l = p

    coeffs = e, 2*d, 3*c, 4*b, 5*a
    dϕdT = evalpoly(T, coeffs)
    
    if p.apply_bounds && (T < T_s || T > T_l)
        dϕdT = 0.0
    end

    return dϕdT
end

# Print info
function show(io::IO, g::MeltingParam_5thOrder)
    return print(
        io,
        "5th order polynomial melting curve: phi = $(NumValue(g.a)) T^5 + $(NumValue(g.b))T^4 + $(NumValue(g.c))T^3 + $(NumValue(g.d))T^2 + $(NumValue(g.e))T + $(NumValue(g.f))  $(Value(g.T_s)) ≤ T ≤ $(Value(g.T_l))",
    )
end
#-------------------------------------------------------------------------

# MeltingParam_4thOrder  -------------------------------------------------------
"""
    MeltingParam_4thOrder(b,c,d,e,f,T_s,T_l)

Uses a 4th order polynomial to describe the melt fraction `phi` between solidus temperature `T_s` and liquidus temperature `T_l`
```math
    \\phi = b T^4 + c T^3 + d T^2 + e T + f  \\textrm{   for   } T_s ≤ T ≤ T_l
```
```math
    \\phi = 1 \\textrm{   if   } T>T_l
```
```math
    \\phi = 0  \\textrm{   if   } T<T_s
```
Temperature `T` is in Kelvin.

![MeltingParam_4thOrder](./assets/img/MeltingParam_4thorder.png)

The default values are for Tonalite experiments from Marxer and Ulmer (2019):
- Marxer, F. & Ulmer, P. (2019) Crystallisation and zircon saturation of calc-alkaline tonalite from the Adamello Batholith at upper crustal conditions: an experimental study. *Contributions Mineral. Petrol.* 174, 84

"""
@with_kw_noshow struct MeltingParam_4thOrder{T,U,U2,U3,U4,U5,U6} <: AbstractMeltingParam{T}
    b::GeoUnit{T,U2} = -7.594512597174117e-10 / K^4
    c::GeoUnit{T,U3} = 3.469192091489447e-06 / K^3
    d::GeoUnit{T,U4} = -0.005923529809260 / K^2
    e::GeoUnit{T,U5} = 4.482855645604745 / K
    f::GeoUnit{T,U6} = -1.268730161921053e+03 * NoUnits
    T_s::GeoUnit{T,U} = 963.15K
    T_l::GeoUnit{T,U} = 1270.15K
    apply_bounds::Bool = true
end
MeltingParam_4thOrder(args...) = MeltingParam_4thOrder(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_4thOrder) # info about the struct
    return MaterialParamsInfo(; Equation=L"\phi = bT^4 + cT^3 + dT^2 + eT + f")
end

# Calculation routines
function (p::MeltingParam_4thOrder)(; T, kwargs...)
    @unpack_val b, c, d, e, f, T_s, T_l = p

    coeffs = f, e, d, c, b
    ϕ = evalpoly(T, coeffs)
    if p.apply_bounds
        if T < T_s
            ϕ = 0.0
        elseif T > T_l
            ϕ = 1.0
        end
    end

    return ϕ
end

function compute_dϕdT(p::MeltingParam_4thOrder; T, kwargs...)
    @unpack_val b, c, d, e, T_s, T_l = p

    coeffs = e, 2*d, 3*c, 4*b
    dϕdT = evalpoly(T, coeffs)
    if p.apply_bounds && (T < T_s || T > T_l)
        dϕdT = 0.0
    end
    return dϕdT
end

function compute_dϕdT!(
    dϕdT::AbstractArray, p::MeltingParam_4thOrder; T::AbstractArray, kwargs...
)
    for i in eachindex(T)
        @inbounds dϕdT[i] = compute_dϕdT(p; T=T[i])
    end

    return nothing
end

# Print info
function show(io::IO, g::MeltingParam_4thOrder)
    return print(
        io,
        "4th order polynomial melting curve: phi = $(NumValue(g.b))T^4 + $(NumValue(g.c))T^3 + $(NumValue(g.d))T^2 + $(NumValue(g.e))T + $(NumValue(g.f))  $(Value(g.T_s)) ≤ T ≤ $(Value(g.T_l))",
    )
end
#-------------------------------------------------------------------------

# MeltingParam_Quadratic  -------------------------------------------------------
"""
    MeltingParam_Quadratic(T_s,T_l)

Quadratic melt fraction parameterisation where melt fraction ``\\phi`` depends only on solidus (``T_s``) and liquidus (``T_l``) temperature:
```math
    \\phi = 1.0 - \\left( {T_l - T} \\over {T_l - T_s} \\right)^2
```
```math
    \\phi = 1.0 \\textrm{ if } T>T_l
```
```math
    \\phi = 0.0 \\textrm{ if } T<T_s
```
Temperature `T` is in Kelvin.

![MeltingParam_Quadratic](./assets/img/MeltingParam_Quadratic.png)

This was used, among others, in Tierney et al. (2016) Geology

"""
@with_kw_noshow struct MeltingParam_Quadratic{T,U} <: AbstractMeltingParam{T}
    T_s::GeoUnit{T,U} = 963.15K
    T_l::GeoUnit{T,U} = 1273.15K
    apply_bounds::Bool = true
end
MeltingParam_Quadratic(args...) = MeltingParam_Quadratic(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_Quadratic) # info about the struct
    return MaterialParamsInfo(; Equation=L"\phi = 1.0 - ((T_l - T)/(T_l - T_s))^2")
end

# Calculation routines
function (p::MeltingParam_Quadratic)(; T, kwargs...)
    @unpack_val T_s, T_l = p

    ϕ = 1.0 - ((T_l - T) / (T_l - T_s))^2
    if p.apply_bounds
        if T > T_l
            ϕ = 1.0
        elseif T < T_s
            ϕ = 0.0
        end
    end
    return ϕ
end

function compute_dϕdT(p::MeltingParam_Quadratic; T, kwargs...)
     @unpack_val T_s, T_l = p

    dϕdT = (2T_l - 2T) / ((T_l - T_s)^2)
    if p.apply_bounds && (T > T_l || T < T_s)
        dϕdT = 0.0
    end
    return dϕdT
end

# Print info
function show(io::IO, g::MeltingParam_Quadratic)
    return print(
        io,
        "Quadratic melting curve:  ϕ = 1.0 - ((Tₗ-T)/(Tₗ-Tₛ))² with Tₛ=$(Value(g.T_s)), Tₗ=$(Value(g.T_l)) ",
    )
end
#-------------------------------------------------------------------------

# MeltingParam_Assimilation  -------------------------------------------------------
"""
    MeltingParam_Assimilation(T_s,T_l,a)

Melt fraction parameterisation that takes the assimilation of crustal host rocks into account, as used by Tierney et al. (2016) based upon a parameterisation of Spera and Bohrson (2001)

Here, the fraction of molten and assimilated host rocks ``\\phi`` depends on the solidus (``T_s``) and liquidus (``T_l``) temperatures of the rocks, as well as on a parameter ``a=0.005``
```math
    X = \\left( {T - T_s} \\over {T_l - T_s} \\right)
```
```math
    \\phi = a \\cdot \\left( \\exp^{2ln(100)X} - 1.0 \\right) \\textrm{ if } X ≤ 0.5
```
```math
    \\phi = 1- a \\cdot \\exp^{2ln(100)(1-X)}  \\textrm{ if } X > 0.5
```
```math
    \\phi = 1.0 \\textrm{ if } T>T_l
```
```math
    \\phi = 0.0 \\textrm{ if } T<T_s
```
Temperature `T` is in Kelvin.

![MeltingParam_Assimilation](./assets/img/MeltingParam_Assimilation.png)

This was used, among others, in Tierney et al. (2016), who employed as default parameters:
```math
   T_s=973.15, T_l=1173.15, a=0.005
```

References
==========
- Spera, F.J., and Bohrson, W.A., 2001, Energy-Constrained Open-System Magmatic Processes I: General Model and Energy-Constrained Assimilation and Fractional Crystallization (EC- AFC) Formulation: Journal of Petrology, v. 42, p. 999–1018.
- Tierney, C.R., Schmitt, A.K., Lovera, O.M., de Silva, S.L., 2016. Voluminous plutonism during volcanic quiescence revealed by thermochemical modeling of zircon. Geology 44, 683–686. https://doi.org/10.1130/G37968.1

"""
@with_kw_noshow struct MeltingParam_Assimilation{T,U,U1} <: AbstractMeltingParam{T}
    T_s::GeoUnit{T,U} = 973.15K
    T_l::GeoUnit{T,U} = 1173.15K
    a::GeoUnit{T,U1} = 0.005NoUnits
    apply_bounds::Bool = true
end
MeltingParam_Assimilation(args...) = MeltingParam_Assimilation(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_Assimilation) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"\phi = f(T_S,T_l,a) taking crustal assimilation into account."
    )
end

# Calculation routines
function (p::MeltingParam_Assimilation)(; T, kwargs...)
    @unpack_val T_s, T_l, a = p

    X = (T - T_s) / (T_l - T_s)

    if X <= 0.5
        ϕ = a * (exp(2 * log(100) * X) - 1.0)
    else
        ϕ = 1.0 - a * exp(2 * log(100) * (1 - X))
    end
    if p.apply_bounds
        if T > T_l
            ϕ = 1.0
        elseif T < T_s
            ϕ = 0.0
        end
    end
    return ϕ
end

function compute_dϕdT(p::MeltingParam_Assimilation; T, kwargs...)
    @unpack_val T_s, T_l, a = p

    X = (T - T_s) / (T_l - T_s)
    dϕdT =
        (
            9.210340371976184 *
            a *
            exp((9.210340371976184 * T - 9.210340371976184 * T_s) / (T_l - T_s))
        ) / (T_l - T_s)
    if X > 0.5
        dϕdT =
            (
                9.210340371976184 *
                a *
                exp(
                    9.210340371976184 +
                    (9.210340371976184 * T_s - 9.210340371976184 * T) / (T_l - T_s),
                )
            ) / (T_l - T_s)
    end

    if p.apply_bounds && (T > T_l || T < T_s)
        dϕdT = 0.0
    end
    return dϕdT
end

# Print info
function show(io::IO, g::MeltingParam_Assimilation)
    return print(
        io, "Quadratic melting assimilation parameterisation after Spera & Bohrson (2001)"
    )
end

#-------------------------------------------------------------------------

# Smooth melting function ------------------------------------------------

"""
    SmoothMelting(; p=MeltingParam_4thOrder(), k_sol=0.2/K,  k_liq=0.2/K)

This smoothens the melting parameterisation ``p`` around the solidus ``T_{sol}`` and liquidus ``T_{liq}``
using a smoothened Heaviside step functions for the solidus:

```math
    H_{sol} =  {1.0 \\over { 1 + \\exp( -2 k_{sol} (T - T_{sol} - {2 \\over k_{sol}}) )  }}
```
and liquidus:
```math
    H_{liq} =  1.0 - {1.0 \\over { 1 + \\exp( -2 k_{liq} (T - T_{liq} + {2 \\over k_{liq}}) )  }}
```
The resulting melt fraction ``\\phi`` is computed from the original melt fraction ``\\phi_0`` (computed using one of the methods above) as:
```math
    \\phi =  \\phi_0 H_{sol} H_{liq} + 1.0 - H_{liq}
```
The width of the smoothening zones is controlled by ``k_{sol}, k_{liq}`` (larger values = sharper boundary).

This is important, as jumps in the derivative ``dϕ/dT`` can cause numerical instabilities in latent heat computations, which is prevented with this smoothening.

Example
====

Let's consider a 4th order parameterisation:
```julia
julia> using GLMakie, GeoParams
julia> p = MeltingParam_4thOrder();
julia> T= collect(650.0:1:1050.) .+ 273.15;
julia> T,phi,dϕdT =  PlotMeltFraction(p,T=T);
```

The same but with smoothening:
```julia
julia> p_s = SmoothMelting(p=MeltingParam_4thOrder(), k_liq=0.21/K);
4th order polynomial melting curve: phi = -7.594512597174117e-10T^4 + 3.469192091489447e-6T^3 + -0.00592352980926T^2 + 4.482855645604745T + -1268.730161921053  963.15 K ≤ T ≤ 1270.15 K with smooth Heaviside function smoothening using k_sol=0.1 K⁻¹·⁰, k_liq=0.11 K⁻¹·⁰
julia> T_s,phi_s,dϕdT_s =  PlotMeltFraction(p_s,T=T);
```

We can create plots of this with:
```julia
julia> plt1 = plot(T.-273.15, phi, ylabel="Melt Fraction ϕ", color=:red, label="original", xlabel="Temperature [C]")
julia> plt1 = plot(plt1, T.-273.15, phi_s,  color=:black, label="smoothened", legend=:bottomright)
julia> plt2 = plot(T.-273.15, dϕdT, ylabel="dϕ/dT", color=:red, label="original", xlabel="Temperature [C]")
julia> plt2 = plot(plt2, T.-273.15, dϕdT_s,  color=:black, label="smoothened", legend=:topright)
julia> plot!(plt1,plt2,   xlabel="Temperature [C]", layout=(2,1))
```
The derivative no longer has a jump now:

![MeltingParam_Smooth](./assets/img/MeltingParam_Smooth.png)

"""
struct SmoothMelting{P,T,U} <: AbstractMeltingParam{T}
    p::P
    k_sol::GeoUnit{T,U}
    k_liq::GeoUnit{T,U}
end

# Set default values:
function SmoothMelting(; p=MeltingParam_4thOrder(), k_sol=0.2 / K, k_liq=0.2 / K)
    k_sol = convert(GeoUnit, k_sol)
    k_liq = convert(GeoUnit, k_liq)
    p = @set p.apply_bounds = false
    return SmoothMelting(p, k_sol, k_liq)
end

SmoothMelting(p::AbstractMeltingParam) = SmoothMelting(; p=p)

# Calculation routines
function (param::SmoothMelting)(; T, kwargs...)
    @unpack_val k_sol, k_liq = param

    ϕ = param.p(; T, kwargs...)     # Melt fraction computed in usual manner

    T_s = param.p.T_s
    H_s =inv(1.0 + exp(-2 * k_sol * (T - T_s - (2 / k_sol))))

    T_l = param.p.T_l
    H_l = 1.0 - inv(1.0 + exp(-2 * k_liq * (T - T_l + (2 / k_liq))))

    # Apply heaviside smoothening above liquidus & below solidus
    ϕ = ϕ * H_s * H_l + 1.0 - H_l

    return ϕ
end


function compute_dϕdT(param::SmoothMelting; T, kwargs...)
    @unpack_val k_sol, k_liq = param

    # compute heaviside functions & derivatives of that vs. T
    T_s = param.p.T_s
    T_l = param.p.T_l

    f_s(T) = inv(1.0 + exp(-2 * k_sol * (T - T_s - (2 / k_sol))))
    f_l(T) = 1.0 - inv(1.0 + exp(-2 * k_liq * (T - T_l + (2 / k_liq))))

    H_s, dHs_dT = value_and_partial(f_s, T)
    H_l, dHl_dT = value_and_partial(f_l, T)

    # melt fraction & derivative
    dϕdT = compute_dϕdT(param.p; T=T)
    ϕ = param.p(; T, kwargs...)

    # The derivative of the function
    # ϕ = ϕ(T)*H_s(T)*H_l(T)  + 1.0 - H_l
    # versus T is
    dϕdT_tot = dϕdT * H_s * H_l + ϕ * dHs_dT * H_l + ϕ * H_s * dHl_dT - dHl_dT

    return dϕdT_tot
end

# Print info
function show(io::IO, g::SmoothMelting)
    param = show(io, g.p)
    return print(
        io,
        " with smooth Heaviside function smoothening using k_sol=$(Value(g.k_sol)), k_liq=$(Value(g.k_liq))",
    )
end
#-------------------------------------------------------------------------


"""
    compute_meltfraction(P,T, p::AbstractPhaseDiagramsStruct)

Computes melt fraction in case we use a phase diagram lookup table. The table should have the column `:meltFrac` specified.
"""
function compute_meltfraction(
    p::PhaseDiagram_LookupTable; P::_T, T::_T, kwargs...
) where {_T}
    return p.meltFrac.(T, P)
end

compute_meltfraction(p::PhaseDiagram_LookupTable, args) = compute_meltfraction(p; args...)
"""
    compute_meltfraction!(ϕ::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T:AbstractArray{<:AbstractFloat}, p::PhaseDiagram_LookupTable)

In-place computation of melt fraction in case we use a phase diagram lookup table. The table should have the column `:meltFrac` specified.
"""
function compute_meltfraction!(p::PhaseDiagram_LookupTable, args) end

"""
    compute_dϕdT(P,T, p::AbstractPhaseDiagramsStruct)

Computes derivative of melt fraction vs T in case we use a phase diagram lookup table. The table should have the column `:meltFrac` specified.
The derivative is computed by finite differencing.
"""
function compute_dϕdT(p::PhaseDiagram_LookupTable; P::_T, T::_T, kwargs...) where {_T}
    dT = (maximum(T) - minimum(T)) * 0.5 * 1e-6 + 1e-6   # 1e-6 of the average T value
    ϕ1 = p.meltFrac.(T .+ dT, P)
    ϕ0 = p.meltFrac.(T, P)
    dϕdT = (ϕ1 - ϕ0) / dT

    return dϕdT
end

compute_dϕdT(p::PhaseDiagram_LookupTable, args) = compute_dϕdT(p; args...)

"""
    compute_dϕdT!(dϕdT::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T:AbstractArray{<:AbstractFloat}, p::PhaseDiagram_LookupTable)

In-place computation of melt fraction in case we use a phase diagram lookup table. The table should have the column `:meltFrac` specified.
The derivative is computed by finite differencing.
"""
# function compute_dϕdT!(
#     dϕdT::AbstractArray{_T},
#     p::PhaseDiagram_LookupTable;
#     P::AbstractArray{_T},
#     T::AbstractArray{_T},
#     kwargs...,
# ) where {_T}
#     dT = (maximum(T) - minimum(T)) / 2.0 * 1e-6 + 1e-6   # 1e-6 of the average T value
#     ϕ1 = p.meltFrac.(T .+ dT, P)
#     ϕ0 = p.meltFrac.(T, P)
#     dϕdT = (ϕ1 - ϕ0) / dT

#     return nothing
# end

# compute_dϕdT!( ϕ::AbstractArray, p::PhaseDiagram_LookupTable, args) = compute_dϕdT!(p; args...)

# fill methods programmatically
for myType in (
    :MeltingParam_Caricchi,
    :MeltingParam_Smooth3rdOrder,
    :MeltingParam_5thOrder,
    :MeltingParam_4thOrder,
    :MeltingParam_Quadratic,
    :MeltingParam_Assimilation,
    :SmoothMelting,
)
    @eval begin
        (p::$(myType))(args) = p(; args...)
        compute_meltfraction(p::$(myType), args) = p(; args...)
        compute_dϕdT(p::$(myType), args) = compute_dϕdT(p; args...)
    end
end

# Computational routines needed for computations with the MaterialParams structure
function compute_meltfraction(s::AbstractMaterialParamsStruct, args)
    if isempty(s.Melting) #in case there is a phase with no melting parametrization
        return zero(typeof(args).types[1])
    else
        return compute_meltfraction(s.Melting[1], args)
    end
end

function compute_dϕdT(s::AbstractMaterialParamsStruct, args)
    if isempty(s.Melting) #in case there is a phase with no melting parametrization
        return zero(typeof(args).types[1])
    else
        return compute_dϕdT(s.Melting[1], args)
    end
end

"""
    ϕ = compute_meltfraction(Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

Computation of melt fraction ϕ for the whole domain and all phases, in case an array with phase properties `MatParam` is provided, along with `P` and `T` arrays.
"""
compute_meltfraction(args::Vararg{Any, N}) where N = compute_param(compute_meltfraction, args...)

"""
    compute_meltfraction(ϕ::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of melt fraction ϕ for the whole domain and all phases, in case an array with phase properties `MatParam` is provided, along with `P` and `T` arrays.
"""
compute_meltfraction!(args::Vararg{Any, N}) where N = compute_param!(compute_meltfraction, args...)

"""
    compute_meltfraction(ϕ::AbstractArray{<:AbstractFloat}, PhaseRatios::Union{NTuple{N,T}, SVector{N,T}}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

Computation of melt fraction ϕ for the whole domain and all phases, in case an array with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)
"""
compute_meltfraction_ratio(args::Vararg{Any, N}) where N = compute_param_times_frac(compute_meltfraction, args...)

"""
    ϕ = compute_dϕdT(Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

Computates the derivative of melt fraction ϕ versus temperature `T` for the whole domain and all phases, in case an array with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This is employed in computing latent heat terms in an implicit manner, for example
"""
compute_dϕdT(args::Vararg{Any, N}) where N = compute_param(compute_dϕdT, args...)

"""
    compute_dϕdT!(ϕ::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

Computates the derivative of melt fraction `ϕ` versus temperature `T`, ``{\\partial \\phi} \\over {\\partial T}`` for the whole domain and all phases, in case an array with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This is employed, for example, in computing latent heat terms in an implicit manner.
"""
compute_dϕdT!(args::Vararg{Any, N}) where N = compute_param!(compute_dϕdT, args...)

end
