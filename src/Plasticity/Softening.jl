abstract type AbstractSoftening end

struct NoSoftening{T} <: AbstractSoftening end

struct LinearSoftening{T} <: AbstractSoftening
    hi::T
    lo::T
    max_value::T
    min_value::T
    damage::T

    function LinearSoftening(hi::T, lo::T, max_value::T, min_value::T) where T
        damage = max_value + (max_value - min_value) / (hi - lo)
        new{T}(hi, lo, max_value, min_value, damage)
    end
end

function (softening::LinearSoftening{T})(x, softening_var) where T
    (; hi, lo, min_value, damage) = softening

    softening_var ≥ hi && return min_value
    softening_var ≤ lo && return T(x)
    
    return softening_var * damage
end

# soft = LinearSoftening(1.0, 0.0, 30e0, 15e0)

# soft(30e0, 0.5)
