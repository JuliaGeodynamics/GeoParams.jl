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
@inline @generated function compute_param(
    fn::F,MatParam::NTuple{N,AbstractMaterialParamsStruct}, Phase::Integer, args
) where {F, N}
    quote
        T = isempty(args) ? 0.0 : zero(typeof(args).types[1])
        out = T
        Base.Cartesian.@nexprs $N i-> out += MatParam[i].Phase == Phase ? fn(MatParam[i], args) : T
    end
end

function compute_param(
    fn::F, MatParam::AbstractVector{AbstractMaterialParamsStruct}, Phase::Integer, args
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
        v = getindex.(values(args), I)      # works for scalars & arrays thanks to overload (above)
        argsi = (; zip(k, v)...)
        rho[I] = compute_param(fn, MatParam, Phases[I], argsi)
    end
end

function compute_param!(
    fn::F,
    rho::AbstractArray,
    MatParam::AbstractVector{AbstractMaterialParamsStruct},
    Phases::AbstractArray{Integer,ndim},
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
