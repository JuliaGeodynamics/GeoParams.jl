module Interpolations

using Adapt

export LinearInterpolator, interpolate

"""
    LinearInterpolator{T, A}

A custom 2D linear interpolation object that works on both CPU and GPU.
Stores the knots (grid points) and coefficients (data values) for interpolation.
"""
struct LinearInterpolator{T, A <: AbstractArray{T, 2}}
    knots::Tuple{AbstractArray{T,1}, AbstractArray{T, 1}}  # (T_vec, P_vec)
    coefs::A                            # 2D data array
end

# Make LinearInterpolator adaptable for GPU arrays
Adapt.@adapt_structure LinearInterpolator


"""
    lerp(a, b, t)

Linear interpolation between values `a` and `b` with weight `t`.
"""
@inline lerp(a, b, t) = a + t * (b - a)


"""
    get_corners(F, i, j)

Get the four corner values needed for bilinear interpolation.
"""
@inline Base.@propagate_inbounds function get_corners(F::AbstractArray{T, 2}, i::Int, j::Int) where T
    i1, j1 = i + 1, j + 1
    @inbounds begin
        a = F[i, j]
        b = F[i1, j]
        c = F[i, j1]
        d = F[i1, j1]
    end
    return a, b, c, d
end

"""
    interpolate(knots::Tuple{AbstractArray{T, 1}, AbstractArray{T, 1}}, data::AbstractArray{T, 2})

Create a 2D linear interpolation object similar to Interpolations.jl's linear_interpolation function.
"""
function interpolate(knots::Tuple{AbstractArray{T, 1}, AbstractArray{T, 1}}, data::AbstractArray{T, 2}) where T
    return LinearInterpolator(knots, data)
end

"""
    (itp::LinearInterpolator)(x, y)

Evaluate the interpolation at point (x, y) using bilinear interpolation.
"""
function (itp::LinearInterpolator{T})(x::Real, y::Real) where T
    # Promote input types to the interpolator's element type
    x_promoted = T(x)
    y_promoted = T(y)

    x_knots, y_knots = itp.knots
    data = itp.coefs

    nx = length(x_knots)
    ny = length(y_knots)

    # Clamp inputs to valid range (Flat extrapolation)
    x_clamped = clamp(x_promoted, x_knots[1], x_knots[end])
    y_clamped = clamp(y_promoted, y_knots[1], y_knots[end])

    # Find indices for clamped values
    i = find_interval(x_knots, x_clamped)
    j = find_interval(y_knots, y_clamped)

    # Ensure indices are within bounds
    i = clamp(i, 1, nx - 1)
    j = clamp(j, 1, ny - 1)

    # Calculate interpolation weights
    t = if x_clamped ≥ x_knots[end]
        one(T)
    else
        (x_clamped - x_knots[i]) / (x_knots[i+1] - x_knots[i])
    end

    s = if y_clamped ≥ y_knots[end]
        one(T)
    else
        (y_clamped - y_knots[j]) / (y_knots[j+1] - y_knots[j])
    end

    # get corner values
    a, b, c, d = get_corners(data, i, j)

    # Interpolate in x direction first
    v0 = lerp(a, b, t)  # Interpolate bottom edge
    v1 = lerp(c, d, t)  # Interpolate top edge

    return lerp(v0, v1, s)  # Final interpolation in y direction
end

# Vectorized evaluation
function (itp::LinearInterpolator{T})(x::AbstractArray, y::AbstractArray) where T
    return itp.(x, y)
end

"""
    find_interval(knots, x)

Find the interval index i such that knots[i] <= x < knots[i+1].
Returns an index that may be out of bounds for extrapolation handling.
"""
function find_interval(knots::AbstractVector{T}, x::T) where T
    n = length(knots)

    # Handle boundary cases
    if x ≥ knots[end]
        return n - 1  # Return n-1 for the last valid interval
    end

    # Binary search
    left, right = 1, n
    while right - left > 1
        mid = (left + right) ÷ 2
        if knots[mid] ≤ x
            left = mid
        else
            right = mid
        end
    end

    return left
end

end
