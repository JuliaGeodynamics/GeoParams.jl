import Base: (:)

@inline average_pow2(x::NTuple{N,T}) where {N,T} = sum(xi^2 for xi in x) / N

@inline function (:)(
    A::Union{SVector{3,T},NTuple{3,T}}, B::Union{SVector{3,T},NTuple{3,T}}
) where {T}
    return (A[1] * B[1] + A[2] * B[2]) + T(2) * (A[3] * B[3])
end

@inline function (:)(
    A::Union{SVector{6,T},NTuple{6,T}}, B::Union{SVector{6,T},NTuple{6,T}}
) where {T}
    return (A[1] * B[1] + A[2] * B[2] + A[3] * B[3]) +
           T(2) * (A[4] * B[4] + A[5] * B[5] + A[6] * B[6])
end

@inline (:)(A::SMatrix{M,M,T,N}, B::SMatrix{M,M,T,N}) where {M,N,T} = sum(A .* B)

@inline second_invariant(A::NTuple{N,T}) where {N,T} = √(0.5 * (A:A))
@inline second_invariant(A::Vararg{N,T}) where {N,T} = second_invariant(A)
@inline second_invariant(A::SMatrix) = √(0.5 * (A:A))
@inline second_invariant(A::SVector) = √(0.5 * (A:A))
@inline second_invariant(A::Matrix{T}) where {T} = √(0.5 * sum(Ai * Ai for Ai in A))

"""
    second_invariant_staggered(Aii::NTuple{2,T}, Axy::NTuple{4,T}) where {T} 

Computes the second invariant of the 2D tensor `A` when its off-diagonal components 
need to be maped from cell center to cell vertex.  `Aii` is a tuple containinig the diagonal
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
need to be maped from cell center to cell vertex. `Aii` is a tuple containinig the diagonal
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
need to be maped from cell center to cell vertex. `Axx`, and `Ayy` are tuples containinig the diagonal
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
need to be maped from cell center to cell vertex. `Axx`, `Ayy`, and `Azz` are tuples containinig the diagonal
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

## Kernels for effective strain (Eij_eff = Eij + Tij/(2 G dt) ) 

@inline _elastic_ε(v::ConstantElasticity, τij_old, dt) = τij_old / (2 * v.G * dt)
@inline _elastic_ε(v::Vararg{Any,N}) where {N} = 0.0

@inline function elastic_ε(v::CompositeRheology, τij_old, dt)
    return nreduce(vi -> _elastic_ε(vi, τij_old, dt), v.elements)
end

@inline effective_ε(εij::T, v, τij_old::T, dt) where {T} = εij + elastic_ε(v, τij_old, dt)

# By letting ::NTuple{N, Union{T, Any}} it can recursively call itself when strain and stress are nested NamedTuples
# i.e. in staggered grids
@inline function effective_ε(
    εij::NTuple{N,Union{T,Any}}, v, τij_old::NTuple{N,Union{T,Any}}, dt
) where {N,T}
    return ntuple(i -> effective_ε(εij[i], v, τij_old[i], dt), Val(N))
end

# 2D wrapper 
function effective_ε(εxx, εyy, εxy, v, τxx_old, τyy_old, τxy_old, dt)
    return effective_ε((εxx, εyy, εxy), v, (τxx_old, τyy_old, τxy_old), dt)
end

function effective_εII(εxx, εyy, εxy, v, τxx_old, τyy_old, τxy_old, dt)
    εxx, εyy, εxy = effective_ε(εxx, εyy, εxy, v, τxx_old, τyy_old, τxy_old, dt)
    εII = second_invariant(εxx, εyy, εxy)
    return εII
end

# 3D wrapper 
function effective_ε(
    εxx,
    εyy,
    εzz,
    εyz,
    εxz,
    εxy,
    v,
    τxx_old,
    τyy_old,
    τzz_old,
    τyz_old,
    τxz_old,
    τxy_old,
    dt,
)
    return effective_ε(
        (εxx, εyy, εzz, εyz, εxz, εxy),
        v,
        (τxx_old, τyy_old, τzz_old, τyz_old, τxz_old, τxy_old),
        dt,
    )
end

function effective_εII(
    εxx,
    εyy,
    εzz,
    εyz,
    εxz,
    εxy,
    v,
    τxx_old,
    τyy_old,
    τzz_old,
    τyz_old,
    τxz_old,
    τxy_old,
    dt,
)
    εxx, εyy, εzz, εyz, εxz, εxy = effective_ε(
        εxx,
        εyy,
        εzz,
        εyz,
        εxz,
        εxy,
        v,
        τxx_old,
        τyy_old,
        τzz_old,
        τyz_old,
        τxz_old,
        τxy_old,
        dt,
    )
    εII = second_invariant(εxx, εyy, εzz, εyz, εxz, εxy)
    return εII
end
