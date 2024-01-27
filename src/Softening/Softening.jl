
abstract type AbstractSoftening end

struct NoSoftening <: AbstractSoftening end

@inline (softening::NoSoftening)(max_value, ::Vararg{Any, N}) where N = max_value
  
## Linear softening
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

## Non linear softening 
# (Thibault et al 2021; https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2021GC009675)
import SpecialFunctions.erfc

@with_kw struct NoLinearSoftening{T} <: AbstractSoftening
    ξ₀::T= 0.0 # maximum value
    Δ::T = 0.0 # amplitude of the softening (i.e. minimum value)
    μ::T = 1.0 # mean of the softening
    σ::T = 0.5 # standard deviation of the softening
end

NoLinearSoftening(args::Vararg{Any, N}) where N = NoLinearSoftening(promote(args...)...)

@inline (softening::NoLinearSoftening)(args::Vararg{Any, N}) where N = softening(promote(args...)...)

@inline function (softening::NoLinearSoftening)(softening_var::T) where T
    return softening.ξ₀ - 0.5 * softening.Δ * erfc(- (softening_var - softening.μ) / softening.σ)
end
