# Stress tensor computations
using StaticArrays
export compute_τij

"""
    τij = compute_τij(εij::NTuple{3,T}, τij_old::NTuple{3,T}, v, args)

computes `τij` for given strain rate values `εij`
"""
# Single material phase, collocated grid
function compute_τij(v, εij::NTuple{N,T}, args, τij_old::NTuple{N,T}) where {T,N}

    # Second invariant of effective strainrate (taking elasticity into account)
    #ε_eff = εij .+ 0.5.*τij_old./(1.0*args.dt)
    ε_eff = effective_ε(εij, v, τij_old, args.dt)
    εII = second_invariant(ε_eff)

    args = merge(args, (τII_old=0,))
    τII = first(compute_τII(v, εII, args))
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff

    return τij, τII
end

# Single material phase, staggered grid
function compute_τij(v, εij::NTuple{N,Union{T,NTuple{4,T}}}, args, τij_old::NTuple{N,Union{T,NTuple{4,T}}}) where {N,T}

    # Second invariant of effective strainrate (taking elasticity into account)
    ε_eff = effective_ε(εij, v, τij_old, args.dt)
    εII = second_invariant(ε_eff...)
    ε_eff_averaged = staggered_tensor_average(ε_eff)

    args = merge(args, (τII_old=0,))
    τII = first(compute_τII(v, εII, args))
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff_averaged

    return τij, τII
end

# Multiple material phases, collocated grid
function compute_τij(
    v::NTuple{N1,AbstractMaterialParamsStruct},
    εij::NTuple{N2,T},
    args,
    τij_old::NTuple{N2,T},
    phase::I,
) where {T,N1,N2,I<:Integer}

    # Second invariant of effective strainrate (taking elasticity into account)
    ε_eff = effective_ε(εij, v, τij_old, args.dt, phase)
    εII = second_invariant(ε_eff)

    args = merge(args, (τII_old=0,))
    τII = nphase(vi -> first(compute_τII(vi, εII, args)), phase, v)
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff

    return τij, τII, η_eff
end

# Multiple material phases, staggered grid
function compute_τij(
    v::NTuple{N1,AbstractMaterialParamsStruct},
    εij::NTuple{N2,Union{T,NTuple{4,T}}},
    args,
    τij_old::NTuple{N2,Union{T,NTuple{4,T}}},
    phases::NTuple{N2,Union{I,NTuple{4,I}}},
) where {T,N1,N2,I<:Integer}

    # Second invariant of effective strainrate (taking elasticity into account)
    ε_eff = effective_ε(εij, v, τij_old, args.dt, phases)
    εII = second_invariant(ε_eff...)
    ε_eff_averaged = staggered_tensor_average(ε_eff)

    τII = nphase(vi -> first(compute_τII(vi, εII, args)), phases[1], v)
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff_averaged

    return τij, τII, η_eff
end

## Hellper functions
@inline function staggered_tensor_average(x::NTuple{N,Union{T,NTuple{4,T}}}) where {N,T}
    ntuple(Val(N)) do i
        Base.@_inline_meta
        staggered_tensor_average(x[i])
    end
end

staggered_tensor_average(x::NTuple{N,T}) where {N,T} = sum(x) / N
staggered_tensor_average(x::T) where {T<:Number} = x
