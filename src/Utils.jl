# Various helper functions (mostly for internal use)

# Finds index in an allocation-free manner
@generated function find_ind(x::NTuple{N,Integer}, k::Integer) where {N}
    quote
        Base.Cartesian.@nexprs $N i -> x[i] == k && return i
        return 0
    end
end

# Find max element in a tuple
@generated function find_max_tuple(x::NTuple{N,T}) where {N,T}
    quote
        max = x[1]
        Base.Cartesian.@nexprs $N i -> max = (i > 1 && x[i] > max) ? x[i] : max
        return max
    end
end

# Find inner tuple of maximum length
function max_length_tuple(t::NTuple{N,Tuple}) where {N}
    return find_max(ntuple(x -> length(t[x]), Val(N)))
end

# broadcast getindex() to NamedTuples
@inline function ntuple_idx(args::NamedTuple, I::Vararg{Integer,N}) where {N}
    k = keys(args)
    v = getindex.(values(args), I...)
    return (; zip(k, v)...)
end

# fast exponential
@inline fastpow(x::Number, n::Integer) = x^n

@inline function fastpow(x::Number, n::AbstractFloat)
    isinteger(n) && return x^Int(n)
    x > 0 && return exp(log(x) * n)
    return x^n
end

@inline function fastpow(x::Quantity, n::AbstractFloat)
    isinteger(n) && return x^Int(n)
    return x^n
end

@inline function pow_check(x::T, n) where T
    if isone(x) || isone(n)
        x
    elseif iszero(n)
        one(T)
    else
        fastpow(x, n)
    end
end

# Tuple iterators
@generated function nreduce(f::F, v::NTuple{N,Any}) where {N,F}
    Base.@_inline_meta
    quote
        val = 0.0
        Base.Cartesian.@nexprs $N i -> val += @inbounds f(v[i])
        return val
    end
end

@generated function nreduce(
    f::F, v::NTuple{N,Any}, id_args::NTuple{N,Integer}, args::NTuple{NT,Any}
) where {N,NT,F<:Function}
    Base.@_inline_meta
    quote
        val = 0.0
        Base.Cartesian.@nexprs $N i -> val += @inbounds f(v[i], args[id_args[i]])
        return val
    end
end

@generated function nphase(f::F, phase::Integer, v::NTuple{N,AbstractMaterialParamsStruct}) where {N,F<:Function}
    Base.@_inline_meta
    quote
        Base.Cartesian.@nexprs $N i -> @inbounds v[i].Phase === phase && return f(v[i])
        return 0.0
    end
end

@generated function nphase_ratio(f::F, phase_ratio::NTuple{N,T}, v::NTuple{N,AbstractMaterialParamsStruct}) where {N,F,T}
    Base.@_inline_meta
    quote
        val = 0.0
        Base.Cartesian.@nexprs $N i -> val += @inbounds f(v[i]) * phase_ratio[i]
        return val
    end
end

# Creates tuple without branching
make_tuple(x) = (x,)
make_tuple(x::Tuple) = x

# Macros 
macro print(a1, a2)
    return :($(esc(a1)) === true ? println($(esc(a2))) : nothing)
end
