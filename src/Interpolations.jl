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

    # Find indices for x (temperature)
    i = find_interval(x_knots, x_promoted)
    # Find indices for y (pressure)
    j = find_interval(y_knots, y_promoted)

    # Get local coordinates
    if i ≥ nx
        i = nx - 1
        t = one(T)
    else
        t = (x - x_knots[i]) / (x_knots[i+1] - x_knots[i])
    end

    if j ≥ ny
        j = length(y_knots) - 1
        s = one(T)
    else
        s = (y - y_knots[j]) / (y_knots[j+1] - y_knots[j])
    end

    # Bilinear interpolation
    v00 = data[i, j]
    v10 = data[i+1, j]
    v01 = data[i, j+1]
    v11 = data[i+1, j+1]

    # Interpolate in x direction first
    v0 = (one(T) - t) * v00 + t * v10
    v1 = (one(T) - t) * v01 + t * v11

    # Then interpolate in y direction
    return (one(T) - s) * v0 + s * v1
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

    # Handle out of bounds
    if x <= knots[1]
        return 1
    elseif x >= knots[end]
        return n
    end

    # Binary search for the interval
    left, right = 1, n
    while right - left > 1
        mid = (left + right) ÷ 2
        if knots[mid] <= x
            left = mid
        else
            right = mid
        end
    end

    return left
end

end
