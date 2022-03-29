module MeltingParam

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, PhaseDiagram_LookupTable, AbstractMaterialParamsStruct
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
        MeltingParam_Quadratic
        
include("../Utils.jl")
include("../Computations.jl")

# Caricchi  -------------------------------------------------------
"""
    MeltingParam_Caricchi()
    
Implements the T-dependent melting parameterisation used by Caricchi et al 
```math  
    \\theta = {(800.0 - (T + 273.15)) \\over 23.0} 
```
```math  
    \\phi_{melt} = {1.0 \\over (1.0 + e^\\theta)}
```

Note that T is in Kelvin.

"""
@with_kw_noshow struct MeltingParam_Caricchi{T,U} <: AbstractMeltingParam{T}
    a::GeoUnit{T,U}              =   800.0K              
    b::GeoUnit{T,U}              =   23.0K
    c::GeoUnit{T,U}              =   273.15K # shift from C to K
end
MeltingParam_Caricchi(args...) = MeltingParam_Caricchi(convert.(GeoUnit,args)...)

function param_info(s::MeltingParam_Caricchi) # info about the struct
    return MaterialParamsInfo(Equation =  L"\phi = {1 \over {1 + \exp( {800-T[^oC] \over 23})}}")
end

# Calculation routines
function compute_meltfraction(p::MeltingParam_Caricchi{_T}, P::Quantity, T::Quantity) where _T
    @unpack_units a,b,c   = p

    θ       =   (a - (T - c))/b
    ϕ       =   1.0./(1.0 .+ exp.(θ))

    return ϕ
end


function compute_meltfraction(p::MeltingParam_Caricchi{_T}, P::_T, T::_T ) where _T
    @unpack_val a,b,c   = p

    θ       =   (a - (T - c))/b
    return 1.0/(1.0 + exp(θ))
end


function compute_meltfraction!(ϕ::AbstractArray{_T}, p::MeltingParam_Caricchi{_T}, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    @unpack_val a,b,c   = p
    
    @. ϕ = 1.0/(1.0 + exp((a-(T-c))/b)) 

    return nothing
end

function compute_dϕdT(p::MeltingParam_Caricchi{_T}, P::Quantity, T::Quantity) where _T
    @unpack_units a,b,c   = p

    # analytically computed derivative (using Symbolics.jl)
    dϕdT    =   exp((a + c - T) / b) / (b*((1.0 + exp((a + c - T) / b))^2))

    return dϕdT
end

function compute_dϕdT(p::MeltingParam_Caricchi{_T}, P::_T, T::_T ) where _T
    @unpack_val a,b,c   = p

    dϕdT    =   exp((a + c - T) / b) / (b*((1.0 + exp((a + c - T) / b))^2))

    return dϕdT
end

function compute_dϕdT!(dϕdT::AbstractArray{_T}, p::MeltingParam_Caricchi{_T}, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    @unpack_val a,b,c   = p
    
    @. dϕdT = exp((a + c - T) / b) / (b*((1.0 + exp((a + c - T) / b))^2))

    return nothing
end

# Print info 
function show(io::IO, g::MeltingParam_Caricchi)  
    print(io, "Caricchi et al. melting parameterization")  
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
@with_kw_noshow struct MeltingParam_5thOrder{T,U,U1,U2,U3,U4,U5,U6} <: AbstractMeltingParam{T}
    a::GeoUnit{T,U1}    =   2.083291971482524e-12/K^5              
    b::GeoUnit{T,U2}    =  -1.239502833666574e-08/K^4
    c::GeoUnit{T,U3}    =   2.938887604687626e-05/K^3 
    d::GeoUnit{T,U4}    =  -0.034711533077108/K^2 
    e::GeoUnit{T,U5}    =   20.425403874539178/K 
    f::GeoUnit{T,U6}    =   -4.790664658179178e+03*NoUnits 
    T_s::GeoUnit{T,U}   =   963.15K
    T_l::GeoUnit{T,U}   =   1388.2K 
end
MeltingParam_5thOrder(args...) = MeltingParam_5thOrder(convert.(GeoUnit,args)...)

function param_info(s::MeltingParam_5thOrder) # info about the struct
    return MaterialParamsInfo(Equation =  L"\phi = aT^5 + bT^4 + cT^3 + dT^2 + eT + f")
end

# Calculation routines
function compute_meltfraction(p::MeltingParam_5thOrder{_T}, P::Quantity, T::Quantity) where _T
    @unpack_units a,b,c,d,e,f,T_s,T_l   = p

    ϕ   =   a*T^5 + b*T^4 + c*T^3 + d*T^2 + e*T + f
    if T<T_s
        ϕ = 0.
    elseif  T>T_l
        ϕ = 1.
    end

    return ϕ
end


function compute_meltfraction(p::MeltingParam_5thOrder{_T}, P::_T, T::_T ) where _T
    @unpack_val a,b,c,d,e,f,T_s,T_l   = p

    ϕ   =   a*T^5 + b*T^4 + c*T^3 + d*T^2 + e*T + f
    if T<T_s
        ϕ = 0.
    elseif  T>T_l
        ϕ = 1.
    end
    return  ϕ

end

function compute_meltfraction!(ϕ::AbstractArray{_T}, p::MeltingParam_5thOrder{_T}, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    @unpack_val a,b,c,d,e,f,T_s,T_l   = p

    @. ϕ   =   a*T^5 + b*T^4 + c*T^3 + d*T^2 + e*T + f
    
    ϕ[T.<T_s] .= 0.
    ϕ[T.>T_l] .= 1.

    return nothing
end

function compute_dϕdT(p::MeltingParam_5thOrder{_T}, P::Quantity, T::Quantity) where _T
    @unpack_units a,b,c,d,e,T_s,T_l   = p

    dϕdT   =   5*a*T^4 + 4*b*T^3 + 3*c*T^2 + 2*d*T + e
    if T<T_s || T>T_l
        dϕdT = 0.
    end

    return dϕdT
end

function compute_dϕdT(p::MeltingParam_5thOrder{_T}, P::_T, T::_T ) where _T
    @unpack_val a,b,c,d,e,T_s,T_l   = p

    dϕdT   =   5*a*T^4 + 4*b*T^3 + 3*c*T^2 + 2*d*T + e
    if T<T_s || T>T_l
        dϕdT = 0.
    end

    return  dϕdT
end

function compute_dϕdT!(dϕdT::AbstractArray{_T}, p::MeltingParam_5thOrder{_T}, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    @unpack_val a,b,c,d,e,f,T_s,T_l   = p

    @. dϕdT   =   5*a*T^4 + 4*b*T^3 + 3*c*T^2 + 2*d*T + e
    
    dϕdT[T.<T_s] .= 0.
    dϕdT[T.>T_l] .= 0.

    return nothing
end

# Print info 
function show(io::IO, g::MeltingParam_5thOrder)  
    print(io, "5th order polynomial melting curve: phi = $(NumValue(g.a)) T^5 + $(NumValue(g.b))T^4 + $(NumValue(g.c))T^3 + $(NumValue(g.d))T^2 + $(NumValue(g.e))T + $(NumValue(g.f))  $(Value(g.T_s)) ≤ T ≤ $(Value(g.T_l))")  
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
    b::GeoUnit{T,U2}    =  -7.594512597174117e-10/K^4
    c::GeoUnit{T,U3}    =   3.469192091489447e-06/K^3 
    d::GeoUnit{T,U4}    =  -0.005923529809260/K^2 
    e::GeoUnit{T,U5}    =   4.482855645604745/K 
    f::GeoUnit{T,U6}    =   -1.268730161921053e+03*NoUnits 
    T_s::GeoUnit{T,U}   =   963.15K
    T_l::GeoUnit{T,U}   =   1270.15K
end
MeltingParam_4thOrder(args...) = MeltingParam_4thOrder(convert.(GeoUnit,args)...)

function param_info(s::MeltingParam_4thOrder) # info about the struct
    return MaterialParamsInfo(Equation =  L"\phi = bT^4 + cT^3 + dT^2 + eT + f")
end

# Calculation routines
function compute_meltfraction(p::MeltingParam_4thOrder{_T}, P::Quantity, T::Quantity) where _T
    @unpack_units b,c,d,e,f,T_s,T_l   = p

    ϕ   =   b*T^4 + c*T^3 + d*T^2 + e*T + f
    if T<T_s
        ϕ = 0.
    elseif  T>T_l
        ϕ = 1.
    end

    return ϕ
end

function compute_meltfraction(p::MeltingParam_4thOrder{_T}, P::_T, T::_T ) where _T
    @unpack_val b,c,d,e,f,T_s,T_l   = p

    ϕ   =    b*T^4 + c*T^3 + d*T^2 + e*T + f
    if T<T_s
        ϕ = 0.
    elseif  T>T_l
        ϕ = 1.
    end
    return  ϕ

end

function compute_meltfraction!(ϕ::AbstractArray{_T}, p::MeltingParam_4thOrder{_T}, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    @unpack_val b,c,d,e,f,T_s,T_l   = p

    @. ϕ   =    b*T^4 + c*T^3 + d*T^2 + e*T + f
    
    ϕ[T.<T_s] .= 0.
    ϕ[T.>T_l] .= 1.

    return nothing
end

function compute_dϕdT(p::MeltingParam_4thOrder{_T}, P::Quantity, T::Quantity) where _T
    @unpack_units b,c,d,e,T_s,T_l   = p

    dϕdT   =   4*b*T^3 + 3*c*T^2 + 2*d*T + e    
    if T<T_s || T>T_l
        dϕdT = 0.
    end

    return dϕdT
end


function compute_dϕdT(p::MeltingParam_4thOrder{_T}, P::_T, T::_T ) where _T
    @unpack_val b,c,d,e,T_s,T_l   = p

    dϕdT   =   4*b*T^3 + 3*c*T^2 + 2*d*T + e    
    if T<T_s || T>T_l
        dϕdT = 0.
    end
    return  dϕdT

end

function compute_dϕdT!(dϕdT::AbstractArray{_T}, p::MeltingParam_4thOrder{_T}, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    @unpack_val b,c,d,e,T_s,T_l   = p

    @.  dϕdT   =   4*b*T^3 + 3*c*T^2 + 2*d*T + e    
    
    dϕdT[T.<T_s] .= 0.
    dϕdT[T.>T_l] .= 0.

    return nothing
end


# Print info 
function show(io::IO, g::MeltingParam_4thOrder)  
    print(io, "4th order polynomial melting curve: phi = $(NumValue(g.b))T^4 + $(NumValue(g.c))T^3 + $(NumValue(g.d))T^2 + $(NumValue(g.e))T + $(NumValue(g.f))  $(Value(g.T_s)) ≤ T ≤ $(Value(g.T_l))")  
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
    T_s::GeoUnit{T,U}   =   963.15K
    T_l::GeoUnit{T,U}   =   1273.15K
end
MeltingParam_Quadratic(args...) = MeltingParam_Quadratic(convert.(GeoUnit,args)...)

function param_info(s::MeltingParam_Quadratic) # info about the struct
    return MaterialParamsInfo(Equation =  L"\phi = 1.0 - ((T_l - T)/(T_l - T_s))^2")
end

# Calculation routines
function compute_meltfraction(p::MeltingParam_Quadratic{_T}, P::Quantity, T::Quantity) where _T
    @unpack_units T_s,T_l   = p

    ϕ   =    1.0 -  ((T_l - T)/(T_l - T_s))^2
    if T>T_l
        ϕ = 1.0
    elseif T<T_s
        ϕ = 0.0
    end
    return ϕ
end

function compute_meltfraction(p::MeltingParam_Quadratic{_T}, P::_T, T::_T ) where _T
    @unpack_val T_s,T_l   = p

    ϕ   =    1.0 -  ((T_l - T)/(T_l - T_s))^2
    if T>T_l
        ϕ = 1.0
    elseif T<T_s
        ϕ = 0.0
    end
    return  ϕ

end

function compute_meltfraction!(ϕ::AbstractArray{_T}, p::MeltingParam_Quadratic{_T}, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    @unpack_val T_s,T_l   = p

    @. ϕ   =   1.0 -  ((T_l - T)/(T_l - T_s))^2
    ϕ[T.<T_s] .= 0.
    ϕ[T.>T_l] .= 1.

    return nothing
end

function compute_dϕdT(p::MeltingParam_Quadratic{_T}, P::Quantity, T::Quantity) where _T
    @unpack_units T_s,T_l   = p

    dϕdT   =    (2T_l - 2T) / ((T_l - T_s)^2)
    if T>T_l || T<T_s
        dϕdT = 0.0
    end
    return dϕdT
end

function compute_dϕdT(p::MeltingParam_Quadratic{_T}, P::_T, T::_T ) where _T
    @unpack_val T_s,T_l   = p

    dϕdT   =    (2T_l - 2T) / ((T_l - T_s)^2)
    if T>T_l || T<T_s
        dϕdT = 0.0
    end
    return  dϕdT
end

function compute_dϕdT!(dϕdT::AbstractArray{_T}, p::MeltingParam_Quadratic{_T}, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    @unpack_val T_s,T_l   = p

    @. dϕdT = (2T_l - 2T) / ((T_l - T_s)^2)
    dϕdT[T.<T_s] .= 0.
    dϕdT[T.>T_l] .= 0.

    return nothing
end

# Print info 
function show(io::IO, g::MeltingParam_Quadratic)  
    print(io, "Quadratic melting curve:  ϕ = 1.0 - ((Tₗ-T)/(Tₗ-Tₛ))² with Tₛ=$(Value(g.T_s)), Tₗ=$(Value(g.T_l)) ")  
end
#-------------------------------------------------------------------------


"""
    compute_meltfraction(P,T, p::AbstractPhaseDiagramsStruct)

Computes melt fraction in case we use a phase diagram lookup table. The table should have the column `:meltFrac` specified.
"""
function compute_meltfraction(p::PhaseDiagram_LookupTable, P::_T,T::_T) where _T
   return p.meltFrac.(T,P)
end

"""
    compute_meltfraction!(ϕ::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T:AbstractArray{<:AbstractFloat}, p::PhaseDiagram_LookupTable)

In-place computation of melt fraction in case we use a phase diagram lookup table. The table should have the column `:meltFrac` specified.
"""
function compute_meltfraction!(ϕ::AbstractArray{_T}, p::PhaseDiagram_LookupTable, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    ϕ[:]    =   p.meltFrac.(T,P)

    return nothing
end

"""
    compute_dϕdT(P,T, p::AbstractPhaseDiagramsStruct)

Computes derivative of melt fraction vs T in case we use a phase diagram lookup table. The table should have the column `:meltFrac` specified.
The derivative is computed by finite differencing.
"""
function compute_dϕdT(p::PhaseDiagram_LookupTable, P::_T,T::_T) where _T  
    
   dT = (maximum(T) - minimum(T))/2.0*1e-6  + 1e-6   # 1e-6 of the average T value
   ϕ1 = p.meltFrac.(T .+ dT ,P)
   ϕ0 = p.meltFrac.(T  ,P)
   dϕdT = (ϕ1-ϕ0)/dT

   return dϕdT
end

"""
    compute_dϕdT!(dϕdT::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T:AbstractArray{<:AbstractFloat}, p::PhaseDiagram_LookupTable)

In-place computation of melt fraction in case we use a phase diagram lookup table. The table should have the column `:meltFrac` specified.
The derivative is computed by finite differencing.
"""
function compute_meltfraction!(dϕdT::AbstractArray{_T}, p::PhaseDiagram_LookupTable, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    
    dT = (maximum(T) - minimum(T))/2.0*1e-6 + 1e-6   # 1e-6 of the average T value
    ϕ1 = p.meltFrac.(T .+ dT ,P)
    ϕ0 = p.meltFrac.(T  ,P)
    dϕdT = (ϕ1-ϕ0)/dT

    return nothing
end

# Computational routines needed for computations with the MaterialParams structure 
function compute_meltfraction(s::AbstractMaterialParamsStruct, P::_T=zero(_T),T::_T=zero(_T)) where {_T}
    if isempty(s.Melting) #in case there is a phase with no melting parametrization
        return zero(_T)
    end
    return compute_meltfraction(s.Melting[1], P,T)
end

function compute_dϕdT(s::AbstractMaterialParamsStruct, P::_T=zero(_T),T::_T=zero(_T)) where {_T}
    if isempty(s.Melting) #in case there is a phase with no melting parametrization
        return zero(_T)
    end
    return compute_dϕdT(s.Melting[1], P,T)
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