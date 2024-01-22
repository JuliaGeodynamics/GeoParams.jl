
abstract type AbstractSoftening end

struct NoSoftening{T} <: AbstractSoftening end

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

function (softening::LinearSoftening)(softening_var)
    (; hi, lo, max_value, min_value, damage) = softening

    softening_var ≥ hi && return min_value
    softening_var ≤ lo && return max_value
    
    return softening_var * damage
end
