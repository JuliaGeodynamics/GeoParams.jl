# Various helper functions (mosty for internal use)

# Finds index in an allocation-free manner
function find_ind(x::NTuple{N,_I}, k::_I) where {N, _I<:Integer}
    @inbounds for i in 1:N
        if x[i] == k
            return i
        end
    end
    return 0
end

# Find max element in a tuple
function find_max_tuple(t::NTuple{N,T}) where {N,T}
    max = t[1]
    @inbounds for i in 2:N
        if t[i] > max
            max = t[i]
        end
    end
    max
end

# Find inner tuple of maximum length
function max_length_tuple(t::NTuple{N, Tuple}) where N
    find_max(ntuple(x->length(t[x]), Val(N)))
end

# broadcast getindex() to NamedTuples
function ntuple_idx(args::NamedTuple, I::Integer...)
    k = keys(args)
    v = getindex.(values(args), Tuple(I)...)
    return (; zip(k, v)...)
end

# fast exponential
function fastpow(x::Number, n::Integer)
    n > 3 && x > 0 && return exp(log(x)*n)
    return x^n
end

function fastpow(x::Number, n::AbstractFloat) 
    x > 0 && return exp(log(x)*n)
    return x^n
end