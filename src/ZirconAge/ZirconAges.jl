module ZirconSaturation

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
import Base.show, GeoParams.param_info
using ..MaterialParameters: MaterialParamsInfo

abstract type AbstractZirconSaturation{T} <: AbstractMaterialParam end

export  compute_zirconsaturation, 
        compute_zirconsaturation!,    # calculation routines
        param_info,
        Tierney    
        
include("../Utils.jl")
include("../Computations.jl")

# Tierney et al.  -------------------------------------------------------
"""
    Tierney_etal()
    
The cummulative fraction of zircon crystallizing ``F_{zr}`` is given by 

```math  
    F_{zrs} = a - b \\exp \\left(-c \\over T \\right)  
```
```math  
    F_{zrs} = 0 \\textrm{  if  } T<T_s
```
```math  
    F_{zrs} = 1 \\textrm{  if  }T>T_{zsat}
```
where ``a``,``b``,``c`` are constants, ``T`` the temperature in Kevin, ``T_s`` is the solidus temperature, and ``T_{zsat}`` the zircon saturation temperature.

The default parameters are from Tierney et al. (2016): ``a=1.62``, ``b=1.8*10^{-4}``, ``c=1e4K``, ``T_s=963.15 K``, ``T_{zrs}=1098.15K``
Note that the supplement of Tierney et al., (2016) has a typo in the value of `b`, and that it is given in percentage, whereas here ``0 \\le F_{zrs} \\le 1``. 

"""
@with_kw_noshow struct Tierney{T,U,U1} <: AbstractZirconSaturation{T}
    a::GeoUnit{T,U}         =   1.62NoUnits              
    b::GeoUnit{T,U}         =   1.8e4NoUnits 
    c::GeoUnit{T,U1}        =   1e4K
    T_zrs::GeoUnit{T,U1}    =   1098.15K    # 825 C (see Weber et al.  Nat Comm 2020)
    T_s::GeoUnit{T,U1}      =   963.15K     # 690 C (see Weber et al.  Nat Comm 2020)
end
Tierney(args...) = Tierney(convert.(GeoUnit,args)...)

function param_info(s::Tierney) # info about the struct
    return MaterialParamsInfo(Equation =  L"F_{zrs} = a - b*exp(-c/T)")
end

# Calculation routine
function compute_zirconsaturation(p::Tierney{_T}, P::Quantity, T::Quantity) where _T
    @unpack_units a,b,c,T_zrs,T_s   = p

    Fzrs = a-b*exp(-c/T)
    if T<T_s
        Fzrs = 1.
    elseif  T>T_zrs
        Fzrs = 0.
    end
    return Fzrs
end

function compute_zirconsaturation(p::Tierney{_T},  P::_T , T::_T ) where _T
    @unpack_val a,b,c,T_zrs,T_s  = p

    Fzrs = a-b*exp(-c/T)
    if T<T_s
        Fzrs = 1.
    elseif  T>T_zrs
        Fzrs = 0.
    end
    return Fzrs
end

function compute_zirconsaturation!(Fzrs::AbstractArray{_T}, p::Tierney{_T}, P::AbstractArray{_T}, T::AbstractArray{_T}) where _T
    @unpack_val a,b,c,T_zrs,T_s   = p
    
    @. Fzrs = a-b*exp(-c/T)

    Fzrs[T.<T_s  ] .= 1.
    Fzrs[T.>T_zrs] .= 0.

    return nothing
end

# Print info 
function show(io::IO, g::Tierney)  
    print(io, "Tierney et al. zircon saturation parameterization: Fzrs = a - b*exp(-c/T)")  
end
#-------------------------------------------------------------------------

# Computational routines needed for computations with the MaterialParams structure 
function compute_zirconsaturation(s::AbstractMaterialParamsStruct, P::_T=zero(_T), T::_T=zero(_T)) where {_T}
    if isempty(s.ZirconSaturation) #in case there is a phase with no zircon saturation
        return zero(_T)
    end

    return compute_zirconsaturation(s.ZirconSaturation[1], P, T)
end

"""
    compute_zirconsaturation!(Fzrs::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of zircon saturation `Fzrs` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `T` arrays.
"""
compute_zirconsaturation(args...)  = compute_param(compute_zirconsaturation, args...)
compute_zirconsaturation!(args...) = compute_param!(compute_zirconsaturation, args...)


end