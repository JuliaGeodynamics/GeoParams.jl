module Elasticity

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractElasticity{T} <: AbstractMaterialParam end

export  compute_elastic_shear_strainrate,      # calculation routines
        compute_elastic_shear_strainrate!,
        param_info,
        ConstantElasticity                      # constant

include("../Computations.jl")
include("../Utils.jl")
        
# ConstantElasticity  -------------------------------------------------------

"""
    ConstantElasticity(ϕ=30, Ψ=0, C=10e6Pa, FluidPressure=false)

Structure that holds parameters for constant, isotropic, linear elasticity.
"""
@with_kw_noshow struct ConstantElasticity{T,U,U1} <: AbstractElasticity{T} 
    G::GeoUnit{T,U}     =   5e10Pa          # Elastic shear modulus
    ν::GeoUnit{T,U1}    =   0.5NoUnits      # Poisson ratio
    K::GeoUnit{T,U}     =   1e11Pa          # Elastic bulk modulus
    E::GeoUnit{T,U}     =   1e11Pa          # Elastic Young's modulus
end
ConstantElasticity(args...) = ConstantElasticity(convert.(GeoUnit,args)...)

function param_info(s::ConstantElasticity) # info about the struct
    return MaterialParamsInfo(Equation = L"Constant elasticity")
end

# Calculation routines
function (s::ConstantElasticity{_T})(; τII::_T=zero(_T), τII_old::_T=zero(_T),  dt::_T=1.0, kwargs...) where _T
    @unpack_val G   = s
    ε_el = (τII-τII_old)/(2.0 * G * dt)    

    return ε_el
end
compute_elastic_shear_strainrate(s::ConstantElasticity{_T}; τII::_T=zero(_T), τII_old::_T=zero(_T), dt::_T=1.0) where _T = s(; τII = τII, τII_old=τII_old, dt=dt)

# Calculation routine
## obsolete?
function (s::ConstantElasticity{_T})(τII::AbstractArray{_T,N}, τII_old::AbstractArray{_T,N}, dt::_T; kwargs...) where {_T,N}
    ε_el = similar(τII) 
    @unpack_val G   = s
        
    @.  ε_el = (τII-τII_old)/(2.0 * G * dt)   

    return ε_el
end


(s::ConstantElasticity{_T})(τII::AbstractArray{_T,N}, τII_old::AbstractArray{_T,N}, dt::_T, args...) where {_T,N} = s(τII, τII_old, dt; args...)
compute_elastic_shear_strainrate(s::ConstantElasticity{_T}, τII::AbstractArray{_T,N}, τII_old::AbstractArray{_T,N}, dt::_T, args...) where {_T,N} = s(τII, τII_old, dt; args...)
## obsolete?

"""
    compute_elastic_shear_strainrate!(ε_el::AbstractArray{_T,N}, s::ConstantElasticity{_T}; τII::AbstractArray{_T,N}, τII_old::AbstractArray{_T,N}, dt::_T, kwargs...) 

Computes the elastic shear strainrate for given deviatoric stress invariants at the previous (`τII_old`) and new (`τII`) timestep, as well as the timestep `dt`  
"""
function compute_elastic_shear_strainrate!(ε_el::AbstractArray{_T,N}, p::ConstantElasticity{_T}; τII::AbstractArray{_T,N}, τII_old::AbstractArray{_T,N}, dt::_T, kwargs...) where {N,_T}
    for i in eachindex(τII)
        @inbounds ε_el[i] = compute_elastic_shear_strainrate(p, τII=τII[i], τII_old=τII_old[i], dt=dt)
    end
    return nothing
end

# Print info 
function show(io::IO, g::ConstantElasticity) 

    print(io, "Linear elasticity with shear modulus: G = $(UnitValue(g.G)), Poison's ratio: ν = $(UnitValue(g.ν)), bulk modulus: K = $(UnitValue(g.K)) and Young's module: E=$(UnitValue(g.E))")  
  
end   
#-------------------------------------------------------------------------


# Computational routines needed for computations with the MaterialParams structure 
function compute_elastic_shear_strainrate(s::AbstractMaterialParamsStruct, args) 
    #@show s, args
    if isempty(s.Elasticity)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return s.Elasticity[1](args)
    end
end

# add methods programmatically
for myType in (:ConstantElasticity,)
    @eval begin
        (s::$(myType))(args)= s(; args...)
        compute_elastic_shear_strainrate(s::$(myType), args) = s(args)
        compute_elastic_shear_strainrate!(ε_el::AbstractArray{_T,N}, s::$(myType){_T}, args) where {_T,N} = compute_elastic_shear_strainrate!(ε_el, s; args...)
    end
end

compute_elastic_shear_strainrate(args...)  = compute_param(compute_elastic_shear_strainrate, args...)
compute_elastic_shear_strainrate!(args...) = compute_param!(compute_elastic_shear_strainrate, args...)


end