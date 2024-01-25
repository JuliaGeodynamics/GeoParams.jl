
abstract type AbstractSoftening end

struct NoSoftening <: AbstractSoftening end

@inline (softening::NoSoftening)(max_value, ::Vararg{Any, N}) where N = max_value
  
struct LinearSoftening{T} <: AbstractSoftening
    hi::T
    lo::T
    max_value::T
    min_value::T
    slope::T

    function LinearSoftening(min_value::T, max_value::T, lo::T, hi::T) where T
        slope = (max_value - min_value) / (hi - lo)
        new{T}(hi, lo, max_value, min_value, slope)
    end
end

LinearSoftening(min_max_values::NTuple{2, T}, lo_hi::NTuple{2, T}) where T = LinearSoftening(min_max_values..., lo_hi...)

@inline (softening::LinearSoftening)(args::Vararg{Any, N}) where N = softening(promote(args...)...)

@inline function (softening::LinearSoftening)(max_value::T, softening_var::T) where T

    softening_var ≥ softening.hi && return softening.min_value
    softening_var ≤ softening.lo && return max_value
    
    return fma((one(T) - softening_var), softening.slope, softening.min_value)
end
