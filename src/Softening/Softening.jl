abstract type AbstractSoftening <: AbstractMaterialParam end
abstract type AbstractNoSoftening <: AbstractMaterialParam end

# struct NoSoftening end
struct NoSoftening <: AbstractSoftening end

@inline (softening::NoSoftening)(::Any, max_value, ::Vararg{Any, N}) where {N} = max_value

## Linear softening
struct LinearSoftening{T1, T2, T3} <: AbstractSoftening
    min_value::T2
    max_value::T2
    lo::T1
    hi::T1
    slope::T3

    function LinearSoftening(min_value::T1, max_value::T1, lo::T2, hi::T2, slope::T3) where {T1, T2, T3}
        hi_GU = convert(GeoUnit, hi)
        lo_GU = convert(GeoUnit, lo)
        max_value_GU = convert(GeoUnit, max_value)
        min_value_GU = convert(GeoUnit, min_value)
        slope_GU = convert(GeoUnit, slope)

        return new{typeof(hi_GU), typeof(max_value_GU), typeof(slope_GU)}(
            min_value_GU, max_value_GU, lo_GU, hi_GU, slope_GU
        )
    end
end

function LinearSoftening(min_value::T1, max_value::T1, lo::T2, hi::T2) where {T1, T2}
    slope = if T1 <: AbstractFloat
        (max_value - min_value) / (hi - lo)
    else
        (max_value - min_value).val / (hi - lo)
    end
    return LinearSoftening(min_value, max_value, lo, hi, slope)
end

LinearSoftening(min_max_values::NTuple{2, T1}, lo_hi::NTuple{2, T2}) where {T1, T2} = LinearSoftening(min_max_values..., lo_hi...)

@inline function (softening::LinearSoftening)(softening_var, max_value)

    @unpack_val lo, hi, min_value, slope = softening

    softening_var ≥ hi && return min_value
    softening_var ≤ lo && return max_value

    return @muladd (1 - softening_var) * slope + min_value
end

## Non linear softening
# (Duretz et al 2021; https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2021GC009675)
using SpecialFunctions

@with_kw_noshow struct NonLinearSoftening{T, U1, U2} <: AbstractSoftening
    ξ₀::GeoUnit{T, U1} = 0.0NoUnits # maximum value
    Δ::GeoUnit{T, U1} = 0.0NoUnits # amplitude of the softening (i.e. minimum value)
    μ::GeoUnit{T, U2} = 1.0NoUnits # mean of the softening
    σ::GeoUnit{T, U2} = 0.5NoUnits # standard deviation of the softening
end

NonLinearSoftening(args::Vararg{Any, N}) where {N} = NonLinearSoftening(convert.(GeoUnit, promote(args...))...)

@inline function (softening::NonLinearSoftening)(softening_var, ::Vararg{Any, N}) where {N}
    @unpack_val ξ₀, Δ, μ, σ = softening
    return ξ₀ - 0.5 * Δ * erfc(- (softening_var - μ) / σ)
end

# Non linear softening from Taras
@with_kw_noshow struct DecaySoftening{T, U1, U2} <: AbstractSoftening
    εref::GeoUnit{T, U1} = 1.0e-13
    n::GeoUnit{T, U2} = 0.1
end

DecaySoftening(args::Vararg{Any, N}) where {N} = DecaySoftening(convert.(GeoUnit, promote(args...))...)

@inline function (softening::DecaySoftening)(softening_var::T, max_value::T) where {T}
    @unpack_val εref, n = softening

    return @pow max_value * inv((softening_var / εref + 1)^n)
end
