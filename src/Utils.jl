# Various helper functions (mosty for internal use)

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
function ntuple_idx(args::NamedTuple, I::Vararg{Integer,N}) where {N}
    k = keys(args)
    v = getindex.(values(args), Tuple(I)...)
    return (; zip(k, v)...)
end

# fast exponential
function fastpow(x::Number, n::Integer)
    n > 3 && x > 0 && return exp(log(x) * n)
    return x^n
end

function fastpow(x::Number, n::AbstractFloat)
    x > 0 && return exp(log(x) * n)
    return x^n
end
