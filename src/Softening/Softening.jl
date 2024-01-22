
abstract type AbstractSoftening end

struct NoSoftening <: AbstractSoftening end

@inline (softening::NoSoftening)(max_value, ::Vararg{Any, N}) where N = max_value
  
struct LinearSoftening{T} <: AbstractSoftening
    hi::T
    lo::T
    max_value::T
    min_value::T
    damage::T

    function LinearSoftening(min_value::T, max_value::T, lo::T, hi::T) where T
        damage = max_value + (max_value - min_value) / (hi - lo)
        new{T}(hi, lo, max_value, min_value, damage)
    end
end

LinearSoftening(min_max_values::NTuple{2, T}, lo_hi::NTuple{2, T}) where T = @inline LinearSoftening(min_max_values..., lo_hi...)

function (softening::LinearSoftening)(max_value, softening_var::T) where T
    (; hi, lo, min_value, damage) = softening

    softening_var ≥ hi && return min_value
    softening_var ≤ lo && return T(max_value)
    
    return softening_var * damage
end
