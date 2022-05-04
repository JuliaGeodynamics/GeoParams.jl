module Plasticity

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractPlasticity{T} <: AbstractMaterialParam end

export  compute_yieldfunction,      # calculation routines
        compute_yieldfunction!,
        param_info,
        DruckerPrager               # constant
        
include("../Computations.jl")
include("../Utils.jl")
        
# DruckerPrager  -------------------------------------------------------
"""
    DruckerPrager(ϕ=30, Ψ=0, C=10e6Pa, FluidPressure=false)

Sets parameters for Drucker-Prager plasticity, where the yield stress ``\\sigma_{y}`` is computed by
```math  
    \\sigma_{y} = (P-P_f)\\tan(ϕ) + C
```
with ``\\phi`` being the friction angle (in degrees), ``C`` cohesion, ``P`` dynamic pressure and ``P_f`` the fluid pressure (both positive under compression).  

*Yielding* occurs when the second invariant of the deviatoric stress tensor, ``\\tau_{II}=(0.5\\tau_{ij}\\tau_{ij})^{0.5}`` touches the yield stress. 
This can be computed with the yield function ``F`` and the plastic flow potential ``Q``, which are respectively given by 
```math  
    F = \\tau_{II} - \\cos(ϕ)C - \\sin(ϕ)(P-P_f)
```
```math  
    Q = \\tau_{II} - \\sin(Ψ)(P-P_f) 
```
Here, Ψ is the dilation angle, which must be zero for incompressible setups.

Plasticity is activated when ``F(\\tau_{II}^{trial})`` (the yield function computed with a trial stress) is >0. In that case, plastic strainrate ``\\dot{\\varepsilon}^{pl}_{ij}`` is computed by:
```math  
    \\dot{\\varepsilon}^{pl}_{ij} =\\dot{\\lambda} {\\partial Q \\over \\partial \\sigma_{ij}}
```
where ``\\dot{\\lambda}`` is a (scalar) that is nonzero and chosen such that the resuling stress gives ``F(\\tau_{II}^{final})=0``, and ``\\sigma_{ij}=-P + \\tau_{ij}`` denotes the total stress tensor.   
        
"""
@with_kw_noshow struct DruckerPrager{T,U,U1} <: AbstractPlasticity{T} 
    ϕ::GeoUnit{T,U}         =   30NoUnits      # Friction angle
    Ψ::GeoUnit{T,U}         =   0NoUnits        # Dilation angle
    C::GeoUnit{T,U1}        =   10e6Pa          # Cohesion
end
DruckerPrager(args...) = DruckerPrager(convert.(GeoUnit,args)...)

function param_info(s::DruckerPrager) # info about the struct
    return MaterialParamsInfo(Equation = L"F = \\tau_{II} - \\cos(ϕ)C - \\sin(ϕ)(P-P_f); Q=\\tau_{II} - \\sin(Ψ)(P-P_f)")
end

# Calculation routines
function (s::DruckerPrager{_T,U,U1})(; P::_T=zero(_T), τII::_T=zero(_T),  Pf::_T=zero(_T), kwargs...) where {_T,U,U1}
    @unpack_val ϕ, C   = s
    sinϕ, cosϕ =  sincosd(ϕ)
    
    F = τII - cosϕ*C - sinϕ*(P-Pf)   # with fluid pressure (set to zero by default)
    
    return F
end

"""
    compute_yieldfunction(s::DruckerPrager; P, τII_old, Pf, kwargs...) 

Computes the plastic yield function ``F`` for a given second invariant of the deviatoric stress tensor `τII`,  `P` the pressure, and `Pf` fluid pressure.
"""
compute_yieldfunction(s::DruckerPrager{_T}; P::_T=zero(_T), τII::_T=zero(_T), Pf::_T=zero(_T)) where _T = s(; P=P, τII=τII, Pf=Pf)

"""
    compute_yieldfunction!(F::AbstractArray{_T,N}, s::DruckerPrager{_T}; P::AbstractArray{_T,N}, τII::AbstractArray{_T,N}, kwargs...) 

Computes the plastic yield function ``F`` for Drucker-Prager plasticity in an in-place manner.
Required input arrays are pressure ``P`` and the second invariant of the deviatoric stress tensor ``τII`` at every point.
"""
function compute_yieldfunction!(F::AbstractArray{_T,N}, s::DruckerPrager{_T}; P::AbstractArray{_T,N}, τII::AbstractArray{_T,N}, kwargs...) where {N,_T}
        
    @inbounds for i in eachindex(P)
        F[i] = compute_yieldfunction(s, P=P[i], τII=τII[i])
    end

    return nothing
end

"""
    compute_yieldfunction!(F::AbstractArray{_T,N}, s::DruckerPrager{_T}; P::AbstractArray{_T,N}, τII::AbstractArray{_T,N},  Pf::AbstractArray{_T,N}, kwargs...) 

Computes the plastic yield function ``F`` for Drucker-Prager plasticity in an in-place manner when fluid pressure is provided.
Required input arrays are pressure ``P``, fluid pressure ``Pf`` and the second invariant of the deviatoric stress tensor ``τII`` at every point.
"""
function compute_yieldfunction!(F::AbstractArray{_T,N}, s::DruckerPrager{_T}; P::AbstractArray{_T,N}, τII::AbstractArray{_T,N},  Pf::AbstractArray{_T,N}, kwargs...) where {N,_T}
        
    @inbounds for i in eachindex(P)
        F[i] = compute_yieldfunction(s, P=P[i], τII=τII[i], Pf=Pf[i])
    end

    return nothing
end


# Print info 
function show(io::IO, g::DruckerPrager) 
    print(io, "Drucker-Prager plasticity with: C = $(UnitValue(g.C)), ϕ = $(UnitValue(g.ϕ))ᵒ, Ψ = $(UnitValue(g.Ψ))ᵒ" )  
end   
#-------------------------------------------------------------------------


# Computational routines needed for computations with the MaterialParams structure 
function compute_yieldfunction(s::AbstractMaterialParamsStruct, args) 
    if isempty(s.Plasticity)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return s.Plasticity[1](args)
    end
end

# add methods programmatically
for myType in (:DruckerPrager,)
    @eval begin
        (s::$(myType))(args)= s(; args...)
        compute_yieldfunction(s::$(myType), args) = s(args)
        compute_yieldfunction!(H::AbstractArray{_T,N}, s::$(myType){_T}, args) where {_T,N} = compute_yieldfunction!(H, s; args...)
    end
end

compute_yieldfunction(args...)  = compute_param(compute_yieldfunction, args...)
compute_yieldfunction!(args...) = compute_param!(compute_yieldfunction, args...)


end