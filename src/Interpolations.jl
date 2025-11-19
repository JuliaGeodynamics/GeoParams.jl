module Interpolations

using Adapt

export LinearInterpolator, interpolate

"""
    LinearInterpolator{T, A}

A custom 2D linear interpolation object that works on both CPU and GPU.
Stores the knots (grid points) and coefficients (data values) for interpolation.
"""
struct LinearInterpolator{T, A <: AbstractArray{T, 2}}
    # knots::Tuple{AbstractArray{T,1}, AbstractArray{T, 1}}  # (T_vec, P_vec)
    T0::T                           # Starting value in T direction
    dT::T                           # Spacing in T direction
    numT::Int                       # Number of knots in T direction
    Tmax::T                         # Maximum value in T direction
    P0::T                           # Starting value in P direction
    dP::T                           # Spacing in P direction
    numP::Int                       # Number of knots in P direction
    Pmax::T                         # Maximum value in P direction
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
function interpolate(T0::T, dT::T, numT::Int, Tmax::T, P0::T, dP::T, numP::Int, Pmax::T, data::AbstractArray{T, 2}) where T
    return LinearInterpolator(T0, dT, numT, Tmax, P0, dP, numP, Pmax, data)
end

"""
    (itp::LinearInterpolator)(x, y)

Evaluate the interpolation at point (x, y) using bilinear interpolation.
"""
function (itp::LinearInterpolator{T})(x::Real, y::Real) where T
    # Promote input types to the interpolator's element type
    x_promoted = T(x)
    y_promoted = T(y)

    # Extract parameters from struct (replaces x_knots, y_knots = itp.knots)
    T0, dT, numT, Tmax = itp.T0, itp.dT, itp.numT, itp.Tmax
    P0, dP, numP, Pmax = itp.P0, itp.dP, itp.numP, itp.Pmax
    data = itp.coefs

    # Clamp inputs to valid range (Flat extrapolation)
    x_clamped = clamp(x_promoted, T0, Tmax)      # T0 instead of x_knots[1]
    y_clamped = clamp(y_promoted, P0, Pmax)      # P0 instead of y_knots[1]

    # Find indices
    i = clamp(floor(Int, (x_clamped - T0) / dT) + 1, 1, numT - 1)
    j = clamp(floor(Int, (y_clamped - P0) / dP) + 1, 1, numP - 1)

    # Calculate knot values on-the-fly
    x_i = T0 + (i - 1) * dT
    x_i1 = x_i + dT
    y_j = P0 + (j - 1) * dP
    y_j1 = y_j + dP

    # Calculate interpolation weights
    t = if x_clamped ≥ Tmax
        one(T)
    else
        (x_clamped - x_i) / (x_i1 - x_i)
    end

    s = if y_clamped ≥ Pmax
        one(T)
    else
        (y_clamped - y_j) / (y_j1 - y_j)
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
