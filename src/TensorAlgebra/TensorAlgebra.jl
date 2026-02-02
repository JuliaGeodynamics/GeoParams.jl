import Base: (:)

# convert tensor form to Voigt notation
@inline tensor2voigt(A::T) where {T} = tensor2voigt(T, A)
@inline tensor2voigt(::Type{SMatrix{2, 2, T, 4}}, A) where {T} = @SVector [A[1, 1], A[2, 2], A[1, 2]]
@inline tensor2voigt(::Type{SMatrix{3, 3, T, 9}}, A) where {T} = @SVector [A[1, 1], A[2, 2], A[3, 3], A[2, 3], A[1, 3], A[1, 2]]

@inline  function tensor2voigt(::Type{<:AbstractMatrix}, A)
    if size(A) == (2, 2)
        return [A[1, 1], A[2, 2], A[1, 2]]

    elseif size(A) == (3, 3)
        return [A[1, 1], A[2, 2], A[3, 3], A[2, 3], A[1, 3], A[1, 2]]
    end
end

# convert from Voigt to tensor notation
@inline voigt2tensor_2x2(A) = @SMatrix [A[1] A[3]; A[3] A[2]]
@inline voigt2tensor_3x3(A) = @SMatrix [A[1] A[6] A[5]; A[6] A[2] A[4]; A[5] A[4] A[3]]

@inline voigt2tensor(A::AbstractMatrix) = A
@inline voigt2tensor(A::Union{SVector{3, T}, NTuple{3, T}}) where {T} = voigt2tensor_2x2(A)
@inline voigt2tensor(A::Union{SVector{6, T}, NTuple{6, T}}) where {T} = voigt2tensor_3x3(A)
@inline function voigt2tensor(A::AbstractVector)
    if length(A) == 3
        return voigt2tensor_2x2(A)

    elseif length(A) == 6
        return voigt2tensor_3x3(A)
    end
end

# tensor averages in staggered grids

@inline average_pow2(x::Union{NTuple{N}, SVector{N}}) where {N} = mapreduce(x -> x^2, +, x) / N

@inline staggered_tensor_average(x) = x
function staggered_tensor_average(x::NTuple{N, Union{Any, NTuple{4}}}) where {N}
    return ntuple(Val(N)) do i
        Base.@_inline_meta
        _staggered_tensor_average(x[i])
    end
end

@inline _staggered_tensor_average(x::NTuple{N}) where {N} = sum(x) / N
@inline _staggered_tensor_average(x::Number) = x

# Methods to compute the invariant of a tensor

@inline function doubledot(A::T, B::T) where {T <: Union{SVector{3}, NTuple{3}}}
    # Include the third diagonal component Tzz = -Txx - Tyy for 2D plane strain
    Azz = -A[1] - A[2]  # Azz = -Axx - Ayy
    Bzz = -B[1] - B[2]  # Bzz = -Bxx - Byy
    return (A[1] * B[1] + A[2] * B[2] + Azz * Bzz) + 2 * (A[3] * B[3])
end

@inline function doubledot(A::T, B::T) where {T <: Union{SVector{6}, NTuple{6}}}
    return (A[1] * B[1] + A[2] * B[2] + A[3] * B[3]) +
        2 * (A[4] * B[4] + A[5] * B[5] + A[6] * B[6])
end

@inline doubledot(A::SMatrix, B::SMatrix) = sum(A .* B)

@inline second_invariant(A::NTuple) = √(0.5 * doubledot(A, A))
@inline second_invariant(A::SMatrix) = √(0.5 * doubledot(A, A))
@inline second_invariant(A::SVector) = √(0.5 * doubledot(A, A))
@inline second_invariant(A::Matrix) = √(0.5 * sum(Ai * Ai for Ai in A))
# So that is differentiable...
@inline second_invariant(xx, yy, xy) = √(0.5 * (xx^2 + yy^2 + (-xx - yy)^2) + xy^2)
@inline second_invariant(xx, yy, zz, yz, xz, xy) = √(0.5 * (xx^2 + yy^2 + zz^2) + xy^2 + yz^2 + xz^2)

"""
    second_invariant_staggered(Aii::NTuple{2,T}, Axy::NTuple{4,T}) where {T}

Computes the second invariant of the 2D tensor `A` when its off-diagonal components
need to be mapped from cell center to cell vertex.  `Aii` is a tuple containinig the diagonal
terms of `A` at the i-th vertex, and `Axy` is a tuple that contains `A_xy` at the cell centers
around the i-th vertex.
"""
@inline function second_invariant_staggered(Aii::NTuple{2}, Axy::NTuple{4})
    return √(0.5 * (Aii[1]^2 + Aii[2]^2 + (-Aii[1] - Aii[2])^2) + average_pow2(Axy))
end

@inline function second_invariant_staggered(Axx, Ayy, Axy::NTuple{4})
    return second_invariant_staggered((Axx, Ayy), Axy)
end
@inline function second_invariant(Axx, Ayy, Axy::NTuple{4})
    return second_invariant_staggered(Axx, Ayy, Axy)
end

"""
    second_invariant_staggered(Aii::NTuple{3,T}, Ayz::NTuple{4,T}, Axz::NTuple{4,T}, Axy::NTuple{4,T}) where {T}

Computes the second invariant of the 2D tensor `A` when its off-diagonal components
need to be mapped from cell center to cell vertex. `Aii` is a tuple containinig the diagonal
terms of `A` at the i-th vertex, and `Ayz`, `Axz`, and `Axy` are tuples that contain the off-diagonal components of the tensor
at the cell centers around the i-th vertex.
"""
@inline function second_invariant_staggered(
        Aii::NTuple{3}, Ayz::NTuple{4}, Axz::NTuple{4}, Axy::NTuple{4}
    )
    return √(
        0.5 * (Aii[1]^2 + Aii[2]^2 + Aii[3]^2) +
            average_pow2(Ayz) +
            average_pow2(Axz) +
            average_pow2(Axy),
    )
end

@inline function second_invariant_staggered(
        Axx::T, Ayy::T, Azz::T, Ayz::NTuple{4}, Axz::NTuple{4}, Axy::NTuple{4}
    ) where {T <: Number}
    return second_invariant_staggered((Axx, Ayy, Azz), Ayz, Axz, Axy)
end
@inline function second_invariant(
        Axx::T, Ayy::T, Azz::T, Ayz::NTuple{4}, Axz::NTuple{4}, Axy::NTuple{4}
    ) where {T <: Number}
    return second_invariant_staggered(Axx, Ayy, Azz, Ayz, Axz, Axy)
end

"""
    second_invariant_staggered(Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Axy::Number) where {T}

Computes the second invariant of the 2D tensor `A` when its diagonal components
need to be mapped from cell center to cell vertex. `Axx`, and `Ayy` are tuples containinig the diagonal
terms of `A` at the cell centers around the i-th vertex., and `Axy` is the xy component at the i-th vertex.
"""
@inline function second_invariant_staggered(
        Axx::NTuple{4}, Ayy::NTuple{4}, Axy::Number
    )
    # Compute Azz = -Axx - Ayy element-wise for plane strain
    Azz = ntuple(i -> -Axx[i] - Ayy[i], Val(4))
    return √(0.5 * (average_pow2(Axx) + average_pow2(Ayy) + average_pow2(Azz)) + Axy^2)
end

"""
    second_invariant_staggered(Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Azz::NTuple{4,T}, Aij::NTuple{3,T}) where {T}

Computes the second invariant of the 2D tensor `A` when its diagonal components
need to be mapped from cell center to cell vertex. `Axx`, `Ayy`, and `Azz` are tuples containinig the diagonal
terms of `A` at the cell centers around the i-th vertex., and `Aij` is a tuple that contains the off-diagonal
components at the i-th vertex.
"""
@inline function second_invariant_staggered(
        Axx::NTuple{4}, Ayy::NTuple{4}, Azz::NTuple{4}, Aij::NTuple{3}
    )
    return √(
        0.5 * (average_pow2(Axx) + average_pow2(Ayy) + average_pow2(Azz)) +
            Aij[1]^2 +
            Aij[2]^2 +
            Aij[3]^2,
    )
end

@inline function second_invariant_staggered(
        Axx::NTuple{4}, Ayy::NTuple{4}, Azz::NTuple{4}, Ayz::T, Axz::T, Axy::T
    ) where {T <: Number}
    return second_invariant_staggered(Axx, Ayy, Azz, (Ayz, Axz, Axy))
end

@inline function second_invariant(
        Axx::NTuple{4}, Ayy::NTuple{4}, Azz::NTuple{4}, Ayz::T, Axz::T, Axy::T
    ) where {T <: Number}
    return second_invariant_staggered(Axx, Ayy, Azz, Ayz, Axz, Axy)
end


# Methods to rotate the elastic stress

@inline rotate_elastic_stress(ω, τ, dt) = _rotate_elastic_stress(ω, staggered_tensor_average(τ), dt)

@inline _rotate_elastic_stress(ω::Union{AbstractVector, NTuple}, τ, dt) = rotate_elastic_stress3D(ω, τ, dt)
@inline _rotate_elastic_stress(ω, τ, dt) = rotate_elastic_stress2D(ω, τ, dt)

"""
    rotate_elastic_stress2D(ω, τ::T, dt) where T

Bi-dimensional rotation of the elastic stress where τ is in the Voig notation
and ω = 1/2(dux/dy - duy/dx)
"""
@inline Base.@propagate_inbounds function rotate_elastic_stress2D(ω, τ, dt)
    θ = ω * dt
    # NOTE: inlining sincos speeds up considerably this kernel but breaks for < 1.8
    s, c = @inline @fastmath sincos(θ)
    # rotate tensor
    return tensor_rotation(τ, s, c)
end

@inline function tensor_rotation(τ, s, c)
    τxx = τ[1]
    τyy = τ[2]
    τxy = τ[3]

    s2, c2 = s^2, c^2
    sc2    = 2 * c * s

    xx = @muladd c2 * τxx + s2 * τyy - sc2 * τxy
    yy = @muladd s2 * τxx + c2 * τyy + sc2 * τxy
    xy = @muladd (c2 - s2) * τxy + c * s * (τxx - τyy)
   
    return xx, yy, xy
end

"""
    rotate_elastic_stress3D(ω, τ::T, dt) where T

Trii-dimensional rotation of the elastic stress where τ is in the Voig notation and
ω = [duz/dy - duy/dz, dux/dz - duz/dx, duy/dx - dux/dy]
"""
# from Anton's talk
@inline Base.@propagate_inbounds function rotate_elastic_stress3D(ωi, τ, dt)
    # vorticity
    ω = √(sum(x^2 for x in ωi))
    # unit rotation axis
    n = SVector{3, Float64}(inv(ω) * ωi[i] for i in 1:3)
    # integrate rotation angle
    θ = dt * 0.5 * ω
    # Euler Rodrigues rotation matrix
    R = rodrigues_euler(θ, n)
    # rotate tensor
    τij = voigt2tensor(τ)
    τij_rot = R * (τij * R')
    return tensor2voigt(τij_rot)
end

# Euler Rodrigues rotation matrix
@inline function rodrigues_euler(θ, n)
    sinθ, cosθ = sincos(θ)
    c0 = sinθ * cosθ
    Base.@nexprs 3 i -> c_i = sinθ * n[i]
    R1 = @SMatrix [
        c0   -c_3   c_2
        c_3    c0  -c_1
        -c_2   c_1    c0
    ]
    R2 = (1.0 - cosθ) .* (n * n')
    return R1 + R2
end
