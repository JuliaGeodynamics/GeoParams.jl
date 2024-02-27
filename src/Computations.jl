using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
using ..Units
using Parameters, Unitful
using StaticArrays

# Computational routines needed for computations with the MaterialParams structure 

# with tuple & vector - apply for all phases in MatParam
# function compute_param!(
#     fn::F, rho::AbstractVector, MatParam::NTuple{N,AbstractMaterialParamsStruct}, args
# ) where {F,N}
#     return rho .= map(x -> fn(x, P, T), MatParam)
# end

# each individual calculation 
# function compute_param(
#     fn::F, MatParam::NTuple{N,AbstractMaterialParamsStruct}, args
# ) where {F,N}
#     return map(x -> fn(x, args), MatParam)
# end

#---------------------------------------------------------------------------------------------------------------------------#
# Computational routines for Phases

# performs computation given a single Phase
@generated function compute_param(
    fn::F, MatParam::NTuple{N,AbstractMaterialParamsStruct}, Phase::Int64, args::Vararg{Any,NA}
) where {F,N,NA}
    quote
        Base.@_inline_meta
        Base.Cartesian.@nexprs $N i ->
            @inbounds (MatParam[i].Phase == Phase) && return fn(MatParam[i], args...)
        return 0.0
    end
end

@generated function compute_param(
    fn::F, MatParam::NTuple{N,AbstractMaterialParamsStruct}, phase_ratios::Union{SVector{N,T}, NTuple{N,T}}, args::Vararg{Any,NA}
) where {F,N,T,NA}
    quote
        Base.@_inline_meta
        x = zero($T)
        Base.Cartesian.@nexprs $N i ->
            @inbounds  x+= fn(MatParam[i], args...) * phase_ratios[i]
        return x
    end
end

@inline function compute_param(
    fn::F, MatParam::AbstractVector{AbstractMaterialParamsStruct}, Phase::Union{SArray, Int64}, args::Vararg{Any,NA}
) where {F,NA}
    return compute_param(fn, Tuple(MatParam), Phase, args...)
end

@inline function compute_param(fn::F, MatParam::AbstractMaterialParam, args::Vararg{Any,NA}) where {F,NA}
    return fn(MatParam, args...)
end

@inline function compute_param(fn::F, MatParam::AbstractMaterialParamsStruct, args::Vararg{Any,NA}) where {F,NA}
    return fn(MatParam, args...)
end

@inline function compute_param!(
    fn::F, rho::AbstractArray, MatParam::AbstractMaterialParam, args
) where {F}
    @inbounds for I in eachindex(rho)
        k = keys(args)
        v = getindex.(values(args), I)      # works for scalars & arrays thanks to overload (above)
        argsi = (; zip(k, v)...)
        rho[I] = compute_param(fn, MatParam, argsi)
    end
end

@inline function compute_param!(
    fn::F,
    rho::AbstractArray{T, ndim},
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    Phases::AbstractArray{M,ndim},
    args,
) where {F, T, ndim, N, M<:Integer}
    @inbounds for I in eachindex(Phases)
        k = keys(args)
        v = getindex.(values(args), I)      # works for scalars & arrays thanks to overload (above)
        argsi = (; zip(k, v)...)
        rho[I] = compute_param(fn, MatParam, Int64(Phases[I]), argsi)
    end
end

function compute_param!(
    fn::F,
    rho::AbstractArray{T, ndim},
    MatParam::Vector{AbstractMaterialParamsStruct},
    Phases::AbstractArray{M,ndim},
    args,
) where {F, T, ndim, M<:Integer}
    return compute_param!(fn, rho, Tuple(MatParam), Int64(Phases), args)
end

#Phase PhaseRatios

function compute_param!(
    fn::F,
    rho::AbstractArray{_T,N},
    MatParam::NTuple{K,AbstractMaterialParamsStruct},
    PhaseRatio::AbstractArray{_T,M},
    args,
) where {F,_T,N,M,K}
    if M != (N + 1)
        error("The PhaseRatios array should have one dimension more than the other arrays")
    end
    #I1 = last(rho)
    @inbounds for I in CartesianIndices(rho)
        k = keys(args)
        v = getindex.(values(args), Tuple(I)...)   # works for scalars & arrays thanks to overload (above)
        argsi = (; zip(k, v)...)
        # unroll the dot product
        val = Ref(zero(_T))
        ntuple(Val(K)) do i
            Base.@_inline_meta
            val[] += PhaseRatio[I.I..., i] * fn(MatParam[i], argsi)
        end
        rho[I] = val[]
    end
end

#Multiplies parameter with the fraction of a phase
@generated function compute_param_times_frac(
    fn::F, PhaseRatios::Union{NTuple{N,T}, SVector{N,T}}, MatParam::NTuple{N,AbstractMaterialParamsStruct}, argsi
) where {F,N,T}
    # # Unrolled dot product
    quote
        val = zero($T)
        Base.Cartesian.@nexprs $N i -> val += @inbounds PhaseRatios[i] * fn(MatParam[i], argsi)
        return val
    end
end
