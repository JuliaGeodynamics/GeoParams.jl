import Base: (:)

# convert tensor form to Voigt notation
@inline tensor2voigt(A::T) where T = tensor2voigt(T, A)
@inline tensor2voigt(::Type{NTuple{3, T}}, A) where T = A[1,1] , A[2,2], A[1,2]
@inline tensor2voigt(::Type{NTuple{6, T}}, A) where T = A[1,1], A[2,2], A[3,3], A[2,3], A[1,3], A[1,2]
@inline tensor2voigt(::Type{SVector{3, T}}, A) where T = @SVector [A[1,1], A[2,2], A[1,2]]
@inline tensor2voigt(::Type{SVector{6, T}}, A) where T = @SVector [A[1,1], A[2,2], A[3,3], A[2,3], A[1,3], A[1,2]]

@inline  function tensor2voigt(::Type{AbstractMatrix}, A)
    if size(A) == (2, 2)
        return [A[1,1], A[2,2], A[1,2]]

    elseif size(A) == (3, 3)
        return [A[1,1], A[2,2], A[3,3], A[2,3], A[1,3], A[1,2]]
    end
end

# convert from Voigt to tensor notation
@inline voigt2tensor_2x2(A) = @SMatrix [A[1] A[3]; A[3] A[2]]
@inline voigt2tensor_3x3(A) = @SMatrix [A[1] A[6] A[5]; A[6] A[2] A[4]; A[5] A[4] A[3]]

@inline voigt2tensor(A::AbstractMatrix) = A
@inline voigt2tensor(A::Union{SVector{3, T}, NTuple{3, T}} ) where T = voigt2tensor_2x2(A)
@inline voigt2tensor(A::Union{SVector{6, T}, NTuple{6, T}} ) where T = voigt2tensor_3x3(A)
@inline function voigt2tensor(A::AbstractVector)
    if length(A) == 3
        return voigt2tensor_2x2(A)

    elseif length(A) == 6
        return voigt2tensor_3x3(A)
    end
end

# tensor averages in staggered grids

@inline average_pow2(x::Union{NTuple{N,T}, SVector{N,T}}) where {N,T} = mapreduce(x->x^2, +, x) / N

@inline staggered_tensor_average(x) = x
function staggered_tensor_average(x::NTuple{N,Union{T,NTuple{4,T}}}) where {N,T}
    ntuple(Val(N)) do i
        Base.@_inline_meta
        _staggered_tensor_average(x[i])
    end
end

@inline _staggered_tensor_average(x::NTuple{N,T}) where {N,T} = sum(x) / N
@inline _staggered_tensor_average(x::T) where {T<:Number} = x

# Methods to compute the invariant of a tensor

@inline function (:)(
    A::Union{SVector{3,T},NTuple{3,T}}, B::Union{SVector{3,T},NTuple{3,T}}
) where {T}
    return (A[1] * B[1] + A[2] * B[2]) + T(2.0) * (A[3] * B[3])
end

@inline function (:)(
    A::Union{SVector{6,T},NTuple{6,T}}, B::Union{SVector{6,T},NTuple{6,T}}
) where {T}
    return (A[1] * B[1] + A[2] * B[2] + A[3] * B[3]) +
           T(2.0) * (A[4] * B[4] + A[5] * B[5] + A[6] * B[6])
end

@inline (:)(A::SMatrix{M,M,T,N}, B::SMatrix{M,M,T,N}) where {M,N,T} = sum(A .* B)

@inline second_invariant(A::NTuple{N,T}) where {N,T} = √(0.5 * (A:A))
@inline second_invariant(A::SMatrix) = √(0.5 * (A:A))
@inline second_invariant(A::SVector) = √(0.5 * (A:A))
@inline second_invariant(A::Matrix{T}) where {T} = √(0.5 * sum(Ai * Ai for Ai in A))
# So that is differentiable...
@inline second_invariant(xx, yy, xy) = √(0.5*(xx^2 + yy^2) + xy^2)
@inline second_invariant(xx, yy, zz, yz, xz, xy) = √(0.5*(xx^2 + yy^2 + zz^2) + xy^2 + yz^2 + xz^2)

"""
    second_invariant_staggered(Aii::NTuple{2,T}, Axy::NTuple{4,T}) where {T} 

Computes the second invariant of the 2D tensor `A` when its off-diagonal components 
need to be mapped from cell center to cell vertex.  `Aii` is a tuple containinig the diagonal
terms of `A` at the i-th vertex, and `Axy` is a tuple that contains `A_xy` at the cell centers
around the i-th vertex.
"""
@inline function second_invariant_staggered(Aii::NTuple{2,T}, Axy::NTuple{4,T}) where {T}
    return √(0.5 * (Aii[1]^2 + Aii[2]^2) + average_pow2(Axy))
end

@inline function second_invariant_staggered(Axx::T, Ayy::T, Axy::NTuple{4,T}) where {T}
    return second_invariant_staggered((Axx, Ayy), Axy)
end
@inline function second_invariant(Axx::T, Ayy::T, Axy::NTuple{4,T}) where {T}
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
    Aii::NTuple{3,T}, Ayz::NTuple{4,T}, Axz::NTuple{4,T}, Axy::NTuple{4,T}
) where {T}
    return √(
        0.5 * (Aii[1]^2 + Aii[2]^2 + Aii[3]^2) +
        average_pow2(Ayz) +
        average_pow2(Axz) +
        average_pow2(Axy),
    )
end

@inline function second_invariant_staggered(
    Axx::T, Ayy::T, Azz::T, Ayz::NTuple{4,T}, Axz::NTuple{4,T}, Axy::NTuple{4,T}
) where {T}
    return second_invariant_staggered((Axx, Ayy, Azz), Ayz, Axz, Axy)
end
@inline function second_invariant(
    Axx::T, Ayy::T, Azz::T, Ayz::NTuple{4,T}, Axz::NTuple{4,T}, Axy::NTuple{4,T}
) where {T}
    return second_invariant_staggered(Axx, Ayy, Azz, Ayz, Axz, Axy)
end

"""
    second_invariant_staggered(Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Axy::Number) where {T} 

Computes the second invariant of the 2D tensor `A` when its diagonal components 
need to be mapped from cell center to cell vertex. `Axx`, and `Ayy` are tuples containinig the diagonal
terms of `A` at the cell centers around the i-th vertex., and `Axy` is the xy component at the i-th vertex.
"""
@inline function second_invariant_staggered(
    Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Axy::Number
) where {T}
    return √(0.5 * (average_pow2(Axx) + average_pow2(Ayy)) + Axy^2)
end

"""
    second_invariant_staggered(Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Azz::NTuple{4,T}, Aij::NTuple{3,T}) where {T} 

Computes the second invariant of the 2D tensor `A` when its diagonal components 
need to be mapped from cell center to cell vertex. `Axx`, `Ayy`, and `Azz` are tuples containinig the diagonal
terms of `A` at the cell centers around the i-th vertex., and `Aij` is a tuple that contains the off-diagonal
components at the i-th vertex.
"""
@inline function second_invariant_staggered(
    Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Azz::NTuple{4,T}, Aij::NTuple{3,T}
) where {T}
    return √(
        0.5 * (average_pow2(Axx) + average_pow2(Ayy) + average_pow2(Azz)) +
        Aij[1]^2 +
        Aij[2]^2 +
        Aij[3]^2,
    )
end

@inline function second_invariant_staggered(
    Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Azz::NTuple{4,T}, Ayz::T, Axz::T, Axy::T
) where {T}
    return second_invariant_staggered(Axx, Ayy, Azz, (Ayz, Axz, Axy))
end

@inline function second_invariant(
    Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Azz::NTuple{4,T}, Ayz::T, Axz::T, Axy::T
) where {T}
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
    sinθ, cosθ = sincos(θ) 
    # rotate tensor
    tensor_rotation(τ, cosθ, sinθ)
end

@inline function tensor_rotation(A, cosθ, sinθ)
    A11 = A[1]
    A22 = A[2]
    A12 = A[3]
    xx  = muladd( sinθ, muladd( sinθ, A22, cosθ * A12), cosθ * muladd( sinθ, A12, cosθ * A11))
    yy  = muladd(-sinθ, muladd(-sinθ, A11, cosθ * A12), cosθ * muladd(-sinθ, A12, cosθ * A22))
    xy  = muladd(-sinθ, muladd( sinθ, A12, cosθ * A11), cosθ * muladd( sinθ, A22, cosθ * A12))
    return xx, yy, xy
end


"""
    rotate_elastic_stress3D(ω, τ::T, dt) where T

Trii-dimensional rotation of the elastic stress where τ is in the Voig notation and 
ω = [duz/dy - duy/dz, dux/dz - duz/dx, duy/dx - dux/dy]
"""
# from Anton's talk
@inline Base.@propagate_inbounds function rotate_elastic_stress3D(ωi, τ::T, dt) where T
    # vorticity
    ω = √(sum(x^2 for x in ωi))
    # unit rotation axis
    n = inv(ω) .* ωi
    # integrate rotation angle
    θ = dt * 0.5 * ω
    # Euler Rodrigues rotation matrix
    R = rodrigues_euler(θ, n)
    # rotate tensor
    τij = voigt2tensor(τ)
    τij_rot = R * τij * R'
    tensor2voigt(T, τij_rot)
end

# Euler Rodrigues rotation matrix
@inline function rodrigues_euler(θ, n)
    sinθ, cosθ = sincos(θ)
    R1 = cosθ*I
    R2 = sinθ.*(@SMatrix [
        0     -n[3]   n[2]
        n[3]   0     -n[1]
       -n[3]   n[1]   0
    ])
    c = (1-cosθ)
    R3 = SMatrix{3, 3, Float64}(c*n[i]*n[j] for i in 1:3, j in 1:3)
    return R1 + R2 + R3
end
