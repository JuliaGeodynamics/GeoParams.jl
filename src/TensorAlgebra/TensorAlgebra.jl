import Base: (:)

@inline (:)(A::NTuple{4,T}, B::NTuple{4,T}) where {T} = (A[1]*B[1] + A[2]*B[2]) + T(2)*(A[3]*B[3] + A[4]*B[4])
@inline (:)(A::NTuple{6,T}, B::NTuple{6,T}) where {T} = (A[1]*B[1] + A[2]*B[2] + A[3]*B[3]) + T(2)*(A[4]*B[4] + A[5]*B[5] + A[6]*B[6])
@inline (:)(A::SMatrix{M, M, T, N}, B::SMatrix{M, M, T, N}) where {M,N,T} = sum(A.*B)

second_invariant(A::NTuple{N,T}) where {N, T} = √(0.5*A:A)
second_invariant(A::SMatrix) = √(0.5*A:A)
second_invariant(A::Matrix{T, 2}) where T = √(0.5*sum(Ai*Ai for Ai in A))