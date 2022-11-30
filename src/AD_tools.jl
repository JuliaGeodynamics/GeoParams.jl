# Hack into ForwardDiff.jl to return both the f'(x) and f(x)
function derivative(f::F, x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    res = f(ForwardDiff.Dual{T}(x, one(x)))
    return res.value, res.partials.values[1]
end

# Return vector-evaluated function and its Jacobian  
function jacobian(f, x::Abstractrray)
    T = typeof(ForwardDiff.Tag(f, eltype(x)))
    result = ForwardDiff.static_dual_eval(T, f, x)
    J = ForwardDiff.extract_jacobian(T, result, x)
    f = extract_value(result)
    return f, J
end

extract_value(result::SVector{N, ForwardDiff.Dual{Tag, T, N}}) where {N,T,Tag} = SVector{N,T}(result[i].value for i in 1:N)
extract_value(result::AbstractArray) = [result[i].value for i in eachindex]

