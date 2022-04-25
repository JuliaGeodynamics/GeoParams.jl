module MeltingParam

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams:
    AbstractMaterialParam, PhaseDiagram_LookupTable, AbstractMaterialParamsStruct
import Base.show, GeoParams.param_info
using ..MaterialParameters: MaterialParamsInfo

abstract type AbstractMeltingParam{T} <: AbstractMaterialParam end

export  compute_meltfraction, 
        compute_meltfraction!,    # calculation routines
        compute_dϕdT,             # derivative of melt fraction versus 
        compute_dϕdT!,
        param_info,
        MeltingParam_Caricchi,    
        MeltingParam_4thOrder,
        MeltingParam_5thOrder,
        MeltingParam_Quadratic,
        MeltingParam_Assimilation
        
        
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
b=23
a=800
c=273.15
```
Which gives a reasonable fit to experimental data of granodioritic composition (Piwinskii and Wyllie, 1968)

References
====
- Simpson G. (2017) Practical finite element modelling in Earth Sciences Using MATLAB.
"""
@with_kw_noshow struct MeltingParam_Caricchi{T,U} <: AbstractMeltingParam{T}
    a::GeoUnit{T,U} = 800.0K
    b::GeoUnit{T,U} = 23.0K
    c::GeoUnit{T,U} = 273.15K # shift from C to K
end
MeltingParam_Caricchi(args...) = MeltingParam_Caricchi(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_Caricchi) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"\phi = {1 \over {1 + \exp( {800-T[^oC] \over 23})}}"
    )
end

# Calculation routines
function (p::MeltingParam_Caricchi)(; T, kwargs...)
    if T isa Quantity
        @unpack_units a, b, c = p
    else
        @unpack_val a, b, c = p
    end
    θ = (a - (T - c)) / b
    return 1.0 / (1.0 + exp(θ))
end

function compute_meltfraction!(ϕ::AbstractArray, p::MeltingParam_Caricchi; T, kwargs...)
    if T isa Quantity
        @unpack_units a, b, c = p
    else
        @unpack_val a, b, c = p
    end

    for i in eachindex(T)
        @inbounds ϕ[i] = p(; T=T[i])
    end
    return nothing
end

function compute_dϕdT(p::MeltingParam_Caricchi; T, kwargs...)
    if T isa Quantity
        @unpack_units a, b, c = p
    else
        @unpack_val a, b, c = p
    end

    dϕdT = exp((a + c - T) / b) / (b * ((1.0 + exp((a + c - T) / b))^2))

    return dϕdT
end

function compute_dϕdT!(dϕdT::AbstractArray, p::MeltingParam_Caricchi; T, kwargs...)
    if T isa Quantity
        @unpack_units a, b, c = p
    else
        @unpack_val a, b, c = p
    end

    @. dϕdT = exp((a + c - T) / b) / (b * ((1.0 + exp((a + c - T) / b))^2))

    return nothing
end

# Print info 
function show(io::IO, g::MeltingParam_Caricchi)
    return print(io, "Caricchi et al. melting parameterization")
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
end
MeltingParam_5thOrder(args...) = MeltingParam_5thOrder(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_5thOrder) # info about the struct
    return MaterialParamsInfo(; Equation=L"\phi = aT^5 + bT^4 + cT^3 + dT^2 + eT + f")
end

# Calculation routines
function (p::MeltingParam_5thOrder)(; T, kwargs...)
    if T isa Quantity
        @unpack_units a, b, c, d, e, f, T_s, T_l = p
    else
        @unpack_val a, b, c, d, e, f, T_s, T_l = p
    end

    ϕ = a * T^5 + b * T^4 + c * T^3 + d * T^2 + e * T + f
    if T < T_s
        ϕ = 0.0
    elseif T > T_l
        ϕ = 1.0
    end

    return ϕ
end

function compute_meltfraction!(
    ϕ::AbstractArray, p::MeltingParam_5thOrder; T::AbstractArray, kwargs...
)
    if T isa Quantity
        @unpack_units a, b, c, d, e, f, T_s, T_l = p
    else
        @unpack_val a, b, c, d, e, f, T_s, T_l = p
    end

    @. ϕ = a * T^5 + b * T^4 + c * T^3 + d * T^2 + e * T + f

    @views ϕ[T .< T_s] .= 0.0
    @views ϕ[T .> T_l] .= 1.0

    return nothing
end

function compute_dϕdT(p::MeltingParam_5thOrder; T::Real, kwargs...)
    if T isa Quantity
        @unpack_units a, b, c, d, e, T_s, T_l = p
    else
        @unpack_val a, b, c, d, e, T_s, T_l = p
    end
    
    dϕdT = 5 * a * T^4 + 4 * b * T^3 + 3 * c * T^2 + 2 * d * T + e
    if T < T_s || T > T_l
        dϕdT = 0.0
    end

    return dϕdT
end

function compute_dϕdT!(
    dϕdT::AbstractArray, p::MeltingParam_5thOrder; T::AbstractArray, kwargs...
)
    if T isa Quantity
        @unpack_units a, b, c, d, e, f, T_s, T_l = p
    else
        @unpack_val a, b, c, d, e, f, T_s, T_l = p
    end

    @. dϕdT = 5 * a * T^4 + 4 * b * T^3 + 3 * c * T^2 + 2 * d * T + e

    @views dϕdT[T .< T_s] .= 0.0
    @views dϕdT[T .> T_l] .= 0.0

    return nothing
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
end
MeltingParam_4thOrder(args...) = MeltingParam_4thOrder(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_4thOrder) # info about the struct
    return MaterialParamsInfo(; Equation=L"\phi = bT^4 + cT^3 + dT^2 + eT + f")
end

# Calculation routines
function (p::MeltingParam_4thOrder)(; T::Real, kwargs...)
    if T isa Quantity
        @unpack_units b, c, d, e, f, T_s, T_l = p
    else
        @unpack_val  b, c, d, e, f, T_s, T_l = p
    end

    ϕ = b * T^4 + c * T^3 + d * T^2 + e * T + f
    if T < T_s
        ϕ = 0.0
    elseif T > T_l
        ϕ = 1.0
    end

    return ϕ
end

function compute_meltfraction!(
    ϕ::AbstractArray, p::MeltingParam_4thOrder; T::AbstractArray, kwargs...
)
    if T isa Quantity
        @unpack_units b, c, d, e, f, T_s, T_l = p
    else
        @unpack_val b, c, d, e, f, T_s, T_l = p
    end

    @. ϕ = b * T^4 + c * T^3 + d * T^2 + e * T + f

    @views ϕ[T .< T_s] .= 0.0
    @views ϕ[T .> T_l] .= 1.0

    return nothing
end

function compute_dϕdT(p::MeltingParam_4thOrder; T::Real, kwargs...)
    if T isa Quantity
        @unpack_units  b, c, d, e, T_s, T_l = p
    else
        @unpack_val b, c, d, e, T_s, T_l = p
    end

    dϕdT = 4 * b * T^3 + 3 * c * T^2 + 2 * d * T + e
    if T < T_s || T > T_l
        dϕdT = 0.0
    end
    return dϕdT
end

function compute_dϕdT!(
    dϕdT::AbstractArray, p::MeltingParam_4thOrder; T::AbstractArray, kwargs...
)
    if T isa Quantity
        @unpack_units b, c, d, e, T_s, T_l = p
    else
        @unpack_val b, c, d, e, T_s, T_l = p
    end

    @. dϕdT = 4 * b * T^3 + 3 * c * T^2 + 2 * d * T + e

    @views dϕdT[T .< T_s] .= 0.0
    @views dϕdT[T .> T_l] .= 0.0

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
This was used, among others, in Tierney et al. (2016) Geology

"""
@with_kw_noshow struct MeltingParam_Quadratic{T,U} <: AbstractMeltingParam{T}
    T_s::GeoUnit{T,U} = 963.15K
    T_l::GeoUnit{T,U} = 1273.15K
end
MeltingParam_Quadratic(args...) = MeltingParam_Quadratic(convert.(GeoUnit, args)...)

function param_info(s::MeltingParam_Quadratic) # info about the struct
    return MaterialParamsInfo(; Equation=L"\phi = 1.0 - ((T_l - T)/(T_l - T_s))^2")
end

# Calculation routines
function (p::MeltingParam_Quadratic)(; T::Real, kwargs...)
    if T isa Quantity
        @unpack_units T_s, T_l = p
    else
        @unpack_val T_s, T_l = p
    end

    ϕ = 1.0 - ((T_l - T) / (T_l - T_s))^2
    if T > T_l
        ϕ = 1.0
    elseif T < T_s
        ϕ = 0.0
    end
    return ϕ
end

function compute_meltfraction!(
    ϕ::AbstractArray, p::MeltingParam_Quadratic; T::AbstractArray, kwargs...
)
    if T isa Quantity
        @unpack_units T_s, T_l = p
    else
        @unpack_val T_s, T_l = p
    end

    @. ϕ = 1.0 - ((T_l - T) / (T_l - T_s))^2
    @views ϕ[T .< T_s] .= 0.0
    @views ϕ[T .> T_l] .= 1.0

    return nothing
end

function compute_dϕdT(p::MeltingParam_Quadratic; T::Real, kwargs...)
    if T isa Quantity
        @unpack_units T_s, T_l = p
    else
        @unpack_val T_s, T_l = p
    end

    dϕdT = (2T_l - 2T) / ((T_l - T_s)^2)
    if T > T_l || T < T_s
        dϕdT = 0.0
    end
    return dϕdT
end

function compute_dϕdT!(
    dϕdT::AbstractArray, p::MeltingParam_Quadratic; T::AbstractArray, kwargs...
)
    if T isa Quantity
        @unpack_units T_s, T_l = p
    else
        @unpack_val T_s, T_l = p
    end

    @. dϕdT = (2T_l - 2T) / ((T_l - T_s)^2)
    @views dϕdT[T .< T_s] .= 0.0
    @views dϕdT[T .> T_l] .= 0.0

    return nothing
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
    T_s::GeoUnit{T,U}   =   973.15K
    T_l::GeoUnit{T,U}   =   1173.15K
    a::GeoUnit{T,U1}    =   0.005NoUnits
end
MeltingParam_Assimilation(args...) = MeltingParam_Assimilation(convert.(GeoUnit,args)...)

function param_info(s::MeltingParam_Assimilation) # info about the struct
    return MaterialParamsInfo(Equation =  L"\phi = f(T_S,T_l,a) taking crustal assimilation into account.")
end

# Calculation routines
function (p::MeltingParam_Assimilation)(; T::Real, kwargs...)
    if T isa Quantity
        @unpack_units T_s,T_l,a   = p
    else
        @unpack_val T_s,T_l,a   = p
    end

    X = (T - T_s)/(T_l - T_s)

    if X <= 0.5
        ϕ   =    a * (exp(2*log(100)*X) - 1.0 )
    else
        ϕ   =    1.0 - a * exp(2*log(100)*(1-X) )
    end
   
    if T>T_l
        ϕ = 1.0
    elseif T<T_s
        ϕ = 0.0
    end
    return ϕ
end

function compute_meltfraction!(
    ϕ::AbstractArray, p::MeltingParam_Assimilation; T::AbstractArray, kwargs...
)
    if T isa Quantity
        @unpack_units T_s,T_l,a   = p
    else
        @unpack_val T_s,T_l,a   = p
    end       

    X =   zero(T)
    @. X = (T - T_s)/(T_l - T_s)
    @. ϕ = a * (exp(2*log(100)*X) - 1.0 )

    @views ϕ[X .> 0.5] .=  1.0 .- a .* exp.( 2.0*log(100).*(1.0 .- X[X.>0.5] )) 

    @views ϕ[T.<T_s] .= 0.
    @views ϕ[T.>T_l] .= 1.

    return nothing
end

function compute_dϕdT(p::MeltingParam_Assimilation; T::Real, kwargs...)
    if T isa Quantity
        @unpack_units T_s,T_l,a   = p
    else
        @unpack_val T_s,T_l,a   = p
    end       
    
    X      =   (T - T_s)/(T_l - T_s)
    dϕdT   =   (9.210340371976184*a*exp((9.210340371976184*T - 9.210340371976184*T_s) / (T_l - T_s))) / (T_l - T_s)
    if X>0.5
        dϕdT   = (9.210340371976184*a*exp(9.210340371976184 + (9.210340371976184*T_s - 9.210340371976184*T) / (T_l - T_s))) / (T_l - T_s)
    end

    if T>T_l || T<T_s
        dϕdT = 0.0
    end
    return dϕdT
end

function compute_dϕdT!(
    dϕdT::AbstractArray, p::MeltingParam_Assimilation; T::AbstractArray, kwargs...
)
    if T isa Quantity
        @unpack_units T_s,T_l,a   = p
    else
        @unpack_val T_s,T_l,a   = p
    end       
            
    X         =    zero(T)
    @. X      =   (T .- T_s)./(T_l .- T_s)
    @. dϕdT   =   (9.210340371976184*a*exp((9.210340371976184*T - 9.210340371976184*T_s) / (T_l - T_s))) / (T_l - T_s)

    @views dϕdT[X.>0.5]  = (9.210340371976184.*a.*exp.(9.210340371976184 .+ (9.210340371976184.*T_s .- 9.210340371976184.*T[X.>0.5]) ./ (T_l - T_s))) ./ (T_l - T_s)
    @views dϕdT[T.<T_s] .= 0.
    @views dϕdT[T.>T_l] .= 0.

    return nothing
end

# Print info 
function show(io::IO, g::MeltingParam_Assimilation)
    return print(
        io,
        "Quadratic melting assimilation parameterisation after Spera & Bohrson (2001)",
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
function compute_meltfraction!(
    ϕ::AbstractArray{_T},
    p::PhaseDiagram_LookupTable;
    P::AbstractArray{_T},
    T::AbstractArray{_T},
    kwargs...,
) where {_T}
    ϕ[:] = p.meltFrac.(T, P)

    return nothing
end

compute_meltfraction!( ϕ::AbstractArray, p::PhaseDiagram_LookupTable, args) = compute_meltfraction!(p; args...)

"""
    compute_dϕdT(P,T, p::AbstractPhaseDiagramsStruct)

Computes derivative of melt fraction vs T in case we use a phase diagram lookup table. The table should have the column `:meltFrac` specified.
The derivative is computed by finite differencing.
"""
function compute_dϕdT(p::PhaseDiagram_LookupTable; P::_T, T::_T, kwargs...) where {_T}
    dT = (maximum(T) - minimum(T)) / 2.0 * 1e-6 + 1e-6   # 1e-6 of the average T value
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
function compute_dϕdT!(
    dϕdT::AbstractArray{_T},
    p::PhaseDiagram_LookupTable;
    P::AbstractArray{_T},
    T::AbstractArray{_T},
    kwargs...,
) where {_T}
    dT = (maximum(T) - minimum(T)) / 2.0 * 1e-6 + 1e-6   # 1e-6 of the average T value
    ϕ1 = p.meltFrac.(T .+ dT, P)
    ϕ0 = p.meltFrac.(T, P)
    dϕdT = (ϕ1 - ϕ0) / dT

    return nothing
end

compute_dϕdT!( ϕ::AbstractArray, p::PhaseDiagram_LookupTable, args) = compute_dϕdT!(p; args...)

# fill methods programatically
for myType in (
    :MeltingParam_Caricchi,
    :MeltingParam_5thOrder,
    :MeltingParam_4thOrder,
    :MeltingParam_Quadratic,
    :MeltingParam_Assimilation,
)
    @eval begin
        (p::$(myType))(args) = p(; args...)
        compute_meltfraction(p::$(myType), args) = p(; args...)
        function compute_meltfraction!(ϕ::AbstractArray, p::$(myType), args)
            return compute_meltfraction!(ϕ, p; args...)
        end
        compute_dϕdT(p::$(myType), args) = compute_dϕdT(p; args...)
        function compute_dϕdT!(dϕdT::AbstractArray, p::$(myType), args)
            return compute_dϕdT!(dϕdT, p; args...)
        end
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
compute_meltfraction(args...) = compute_param(compute_meltfraction, args...)

"""
    compute_meltfraction(ϕ::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of melt fraction ϕ for the whole domain and all phases, in case an array with phase properties `MatParam` is provided, along with `P` and `T` arrays.
"""
compute_meltfraction!(args...) = compute_param!(compute_meltfraction, args...)

"""
    ϕ = compute_dϕdT(Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

Computates the derivative of melt fraction ϕ versus temperature `T` for the whole domain and all phases, in case an array with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This is employed in computing latent heat terms in an implicit manner, for example
"""
compute_dϕdT(args...) = compute_param(compute_dϕdT, args...)

"""
    compute_dϕdT!(ϕ::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

Computates the derivative of melt fraction `ϕ` versus temperature `T`, ``\\partial \\phi} \\over {\\partial T}`` for the whole domain and all phases, in case an array with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This is employed in computing latent heat terms in an implicit manner, for example
"""
compute_dϕdT!(args...) = compute_param!(compute_dϕdT, args...)

end
