using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
using ..Units
using Parameters, Unitful

# Computational routines needed for computations with the MaterialParams structure 

# with tuple & vector - apply for all phases in MatParam
function compute_param!(
    fn::F,
    rho::AbstractVector,
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    args,
) where {F,N}
    return rho .= map(x -> fn(x, P, T), MatParam)
end

# each individual calcuation 
function compute_param(
    fn::F, MatParam::NTuple{N,AbstractMaterialParamsStruct}, args
) where {F,N}
    return map(x -> fn(x, args), MatParam)
end

#---------------------------------------------------------------------------------------------------------------------------#
#Computational routines for Phases

# performs computation given a single Phase
@inline function compute_param(
    fn::F, MatParam::Tuple{N,AbstractMaterialParamsStruct}, Phase::Int64, args
) where {F,N}
    Phase_tup = ntuple(i -> MatParam[i].Phase, Val(N))
    idx = find_ind(Phase_tup, Phase)
    T = isempty(args) ? 0.0 : zero(typeof(args).types[1])
    out = ntuple(Val(N)) do i
        Base.@_inline_meta
        if i == idx
            return fn(MatParam[i], args)
        else
            return T
        end
    end
    return out[idx]
end

function compute_param(
    fn::F, MatParam::AbstractVector{AbstractMaterialParamsStruct}, Phase::Int64, args
) where F
    return compute_param(fn, Tuple(MatParam), Phase, args)
end

@inline function compute_param!(
    fn::F,
    rho::AbstractArray,
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    Phases::AbstractArray{<:Integer,ndim},
    args,
) where {F,ndim,N}
    @inbounds for I in eachindex(Phases)
        k = keys(args)
        v = getindex.(values(args), I)
        argsi = (; zip(k, v)...)
        rho[I] = compute_param(fn, MatParam, Phases[I], argsi)
    end
end

function compute_param!(
    fn::F,
    rho::AbstractArray,
    MatParam::AbstractVector{AbstractMaterialParamsStruct},
    Phases::AbstractArray{<:Integer,ndim},
    args,
) where {F,ndim}
    return compute_param!(fn, rho, Tuple(MatParam), Phases, args)
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

    @inbounds for I in CartesianIndices(rho)
        k = keys(args)
        v = getindex.(values(args), Tuple(I)...)
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
function compute_param_times_frac(
    fn::F,
    PhaseRatios::NTuple{N,T},
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    argsi,
) where {F,N,T}
    # Unrolled dot product
    val = Ref(zero(T))
    ntuple(Val(N)) do i
        Base.@_inline_meta
        val[] += PhaseRatios[i] * fn(MatParam[i], argsi)
    end
    return val[]
end
