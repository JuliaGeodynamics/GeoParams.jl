import Base: (:)

@generated function average_pow2(x::NTuple{N, T}) where {N,T}
    quote
        Base.@_inline_meta
        val = zero($T)
        Base.Cartesian.@nexprs $N i -> val = (@inbounds xi=x[i]; muladd(xi, xi, val))
        return (1/$N)*val
    end 
end

@inline (:)(A::NTuple{4,T}, B::NTuple{4,T}) where {T} = (A[1]*B[1] + A[2]*B[2]) + T(2)*(A[3]*B[3] + A[4]*B[4])
@inline (:)(A::NTuple{6,T}, B::NTuple{6,T}) where {T} = (A[1]*B[1] + A[2]*B[2] + A[3]*B[3]) + T(2)*(A[4]*B[4] + A[5]*B[5] + A[6]*B[6])
@inline (:)(A::SMatrix{M, M, T, N}, B::SMatrix{M, M, T, N}) where {M,N,T} = sum(A.*B)

second_invariant(A::NTuple{N,T}) where {N, T} = √(0.5*A:A)
second_invariant(A::SMatrix) = √(0.5*A:A)
second_invariant(A::Matrix{T}) where T = √(0.5*sum(Ai*Ai for Ai in A))

"""
    second_invariant_staggered(Aii::NTuple{2,T}, Axy::NTuple{4,T}) where {T} 

    Computes the second invariant of the 2D tensor `A` when its off-diagonal components 
    need to be maped from cell center to cell vertex.  `Aii` is a tuple containinig the diagonal
    terms of `A` at the i-th vertex, and `Axy` is a tuple that contains `A_xy` at the cell centers
    around the i-th vertex.
"""
function second_invariant_staggered(Aii::NTuple{2,T}, Axy::NTuple{4,T}) where {T} 
    return √(0.5*(Aii[1]^2 + Aii[2]^2) + average_pow2(Axy))
end

"""
    second_invariant_staggered(Aii::NTuple{3,T}, Ayz::NTuple{8,T}, Axz::NTuple{8,T}, Axy::NTuple{8,T}) where {T} 

    Computes the second invariant of the 2D tensor `A` when its off-diagonal components 
    need to be maped from cell center to cell vertex. `Aii` is a tuple containinig the diagonal
    terms of `A` at the i-th vertex, and `Ayz`, `Axz`, and `Axy` are tuples that contain the off-diagonal components of the tensor
    at the cell centers around the i-th vertex.
"""
function second_invariant_staggered(Aii::NTuple{3,T}, Ayz::NTuple{8,T}, Axz::NTuple{8,T}, Axy::NTuple{8,T}) where {T} 
    return √(0.5*(Aii[1]^2 + Aii[2]^2 + Aii[3]^2) + average_pow2(Ayz) + average_pow2(Axz) + average_pow2(Axy))
end

"""
    second_invariant_staggered(Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Axy::Number) where {T} 

    Computes the second invariant of the 2D tensor `A` when its diagonal components 
    need to be maped from cell center to cell vertex. `Axx`, and `Ayy` are tuples containinig the diagonal
    terms of `A` at the cell centers around the i-th vertex., and `Axy` is the xy component at the i-th vertex.
"""
function second_invariant_staggered(Axx::NTuple{4,T}, Ayy::NTuple{4,T}, Axy::Number) where {T} 
    return √(0.5*(average_pow2(Axx) + average_pow2(Ayy)) + Axy^2)
end

"""
second_invariant_staggered(Axx::NTuple{8,T}, Ayy::NTuple{8,T}, Azz::NTuple{8,T}, Aij::NTuple{3,T}) where {T} 

    Computes the second invariant of the 2D tensor `A` when its diagonal components 
    need to be maped from cell center to cell vertex. `Axx`, `Ayy`, and `Azz` are tuples containinig the diagonal
    terms of `A` at the cell centers around the i-th vertex., and `Aij` is a tuple that contains the off-diagonal
    components at the i-th vertex.
"""
function second_invariant_staggered(Axx::NTuple{8,T}, Ayy::NTuple{8,T}, Azz::NTuple{8,T}, Aij::NTuple{3,T}) where {T} 
    return √(0.5*(average_pow2(Axx) + average_pow2(Ayy) + average_pow2(Azz)) + Aij[1]^2 + Aij[2]^2 + Aij[3]^2)
end