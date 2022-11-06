# Stress tensor computations
using StaticArrays
export compute_τij

"""
    τij = compute_τij(εij::NTuple{3,T}, τij_old::NTuple{3,T}, v, args)

computes `τij` for given strain rate values `εij`
"""
function compute_τij(v, εij::NTuple{3,T}, args, τij_old::NTuple{3,T}, G::NTuple{3,T}) where T

    dt    = args.dt
    ε_eff = εij .+ 0.5.*τij_old./(G.*dt)     # effective strain rate
    εII   = second_invariant(ε_eff)
    
    args  = merge(args, (τ0=0,))
    τII   = compute_τII(v, εII, args)
    η_eff = 0.5*τII/εII
    τij   = 2*η_eff.*ε_eff

    return τij, τII
end