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

macro pow(ex)
    substitute_walk(ex)
    esc(:($ex))
end

@inline function substitute_walk(ex::Expr)
    for (i, arg) in enumerate(ex.args)
        new_arg = substitute_walk(arg)
        if !isnothing(new_arg)
            ex.args[i] = new_arg
        end 
    end
end
@inline substitute_walk(sym::Symbol) = sym == :(^) ? :(pow_check) : sym
@inline substitute_walk(x) = x

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

"""
    str2tuple(x::String) 
Converts a string to a tuple with fixed length
"""
function str2tuple(str::String) 
    N = 100
    if length(str) > N
        error("Name String is too long; max. allowed length=$N")
    end
    str_padded = rpad(str, N)
    cu_iterator = codeunits(str_padded)

    return ntuple(i -> cu_iterator[i], Val(N))
end

str2tuple(not_a_str) = not_a_str

function uint2str(x::Vector{UInt8})
    idx = 101-findfirst(x[i] != 0x20 for i in 100:-1:1)
    return String(x[1:idx])
end

@inline uint2str(x::NTuple{100, UInt8}) = uint2str([x...])
# uint2str(x::AbstractCreepLaw) = uint2str(x.Name)

# Creates tuple without branching
make_tuple(x) = (x,)
make_tuple(x::Tuple) = x

# Macros 
macro print(a1, a2)
    return :($(esc(a1)) === true ? println($(esc(a2))) : nothing)
end

# Deriving the given function f with initial guess x using ForwardDiff 
# while returning the value for the function and its derivative in df as Dual number
@inline function value_and_partial(f::F, x::R) where {F, R<:Real}     
    T = typeof(ForwardDiff.Tag(f, R))     
    df =  f(ForwardDiff.Dual{T}(x, one(x)), )     
    return df.value, ForwardDiff.extract_derivative(T, df) 
end

macro extractors(type, field)
    esc(quote
        add_extractor_functions($type, $field)
    end)
end

function add_extractor_functions(::Type{_T}, param_field) where _T
    fields = fieldnames(_T)
    for f in fields
        fun = Symbol("get_$(string(f))")
        checker = !isdefined(GeoParams, fun)
        @eval GeoParams begin 
            $fun(a::$_T) = a.$(f).val
            if $checker
                $fun(a::AbstractMaterialParamsStruct) = isempty(a.$(param_field)) ? 0.0 : $(fun)(a.$(param_field)[1])
                $fun(a::NTuple{N, AbstractMaterialParamsStruct}, phase) where N = nphase($(fun), phase, a)
            end
        end
    end
end