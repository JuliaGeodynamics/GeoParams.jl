# This holds structures and computational routines for compositional rheologies
using StaticArrays
using Setfield

export CompositeRheology, Parallel, create_rheology_string, print_rheology_matrix
export time_τII_0D, compute_εII_harmonic, compute_τII_AD, isplastic, compute_p_τII, local_iterations_εvol, compute_p_harmonic
export computeViscosity_εII, computeViscosity_εII_AD, compute_yieldfunction, compute_elements_εII
export isvolumetricplastic
import Base.getindex

import GeoParams.Units: nondimensionalize, dimensionalize
import GeoParams: nreduce

include("Parallel.jl")              # all related to the Parallel struct
include("CompositeRheology.jl")     # all related to CompositeRheology struct
include("NonlinearIterations.jl")   # nonlinear local iterations
include("Viscosity.jl")             # viscosity computations

# Define rules to nondimensionalise this 
function nondimensionalize(MatParam::Union{Parallel,CompositeRheology}, g::GeoUnits{TYPE}) where {TYPE}
    field_new = ();
    field = MatParam.elements;
    for i=1:length(MatParam.elements)
        field_nd  = nondimensionalize(field[i], g) 
        field_new =  tuple(field_new..., field_nd)
    end
    MatParam = set(MatParam, Setfield.PropertyLens{:elements}(), field_new)

    return MatParam
end

function dimensionalize(MatParam::Union{Parallel,CompositeRheology}, g::GeoUnits{TYPE}) where {TYPE}
    field_new = ();
    field = MatParam.elements;
    for i=1:length(MatParam.elements)
        field_nd  = dimensionalize(field[i], g) 
        field_new =  tuple(field_new..., field_nd)
    end
    MatParam = set(MatParam, Setfield.PropertyLens{:elements}(), field_new)

    return MatParam
end

@inline function compute_τII!(
    τII::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    εII::AbstractArray{T,nDim},
    args;
    tol=1e-6, verbose=false
    ) where {T,nDim,N}
    for I in eachindex(τII)
        τII[I] = first(compute_τII(v, εII[I], (; zip(keys(args), getindex.(values(args), I))...)))
    end
end


# VISCOSITY COMPUTATIONS FROM COMPOSITE_RHEOLOGY & FRIENDS
""" 
    η = computeViscosity_εII(v::Union{Parallel{T,N}, CompositeRheology{T,N}, AbstractConstitutiveLaw}, εII::_T, args; tol=1e-6, verbose=false)

This computes the effective viscosity for a given input rheology `v` and strainrate `εII`
"""
function computeViscosity_εII(
    v::Union{Parallel, CompositeRheology}, 
    εII::_T, 
    args;
    tol=1e-6, verbose=false
) where {_T}
    τII, = compute_τII(v, εII, args; tol=tol, verbose=verbose)
    η    = _T(0.5) * τII * inv(εII)
    return η
end

function computeViscosity_εII(v::T, εII::_T, args; tol=1e-6, verbose=false) where {T<:AbstractConstitutiveLaw,_T}
    τII, = compute_τII(v, εII, args)
    η    = 0.5 * τII * inv(εII)
    return η
end

""" 
    η = computeViscosity_εII_AD(v::Union{Parallel{T,N}, CompositeRheology{T,N}, AbstractConstitutiveLaw}, εII::_T, args; tol=1e-6, verbose=false)

This computes the effective viscosity for a given input rheology `v` and strainrate `εII`, while using AD if necessary
"""
function computeViscosity_εII_AD(
    v::Union{Parallel, CompositeRheology, AbstractConstitutiveLaw}, 
    εII::_T, 
    args;
    tol=1e-6, verbose=false
) where {_T}
    τII = compute_τII_AD(v, εII, args; tol=tol, verbose=verbose)
    η   = _T(0.5) * τII * inv(εII)
    return η
end

function computeViscosity_εII_AD(v::T, εII::_T, args; tol=1e-6, verbose=false) where {T<:AbstractConstitutiveLaw,_T}
    return computeViscosity_εII(v, εII, args) 
end

@inline function elastic_ε(v::Union{CompositeRheology, Parallel}, τij_old, dt)
    return nreduce(vi -> _elastic_ε(vi, τij_old, dt), v.elements)
end

@inline function elastic_ε(v::NTuple{N,AbstractMaterialParamsStruct}, τij_old, dt, phase::Int64) where N
    return nphase(vi -> _phase_elastic_ε(vi.CompositeRheology[1].elements, τij_old, dt), phase, v)
end

# we do need this kernel because one can't have nested generators/closures in a
# @generated function (i.e. a nreduce call inside nphase breaks the function 
# purity and doesnt work)
@generated function _phase_elastic_ε(v::NTuple{N,Any}, τij_old, dt) where {N}
    Base.@_inline_meta
    quote
        val = 0.0
        Base.Cartesian.@nexprs $N i -> @inbounds val += _elastic_ε(v[i], τij_old, dt)
        return val
    end
end
