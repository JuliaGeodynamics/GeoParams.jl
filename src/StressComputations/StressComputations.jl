# Stress tensor computations
using StaticArrays
export compute_τij

"""
    τij = compute_τij(εij::NTuple{3,T}, τij_old::NTuple{3,T}, v, args)

computes `τij` for given strain rate values `εij`
"""
function compute_τij(v, εij::NTuple{3,T}, args, τij_old::NTuple{3,T}=(0.0,0.0,0.0)) where T

    # Second invariant of effective strainrate (taking elasticity into account)
    #ε_eff = εij .+ 0.5.*τij_old./(1.0*args.dt)
    ε_eff = effective_ε(εij[1], εij[2], εij[3], v, τij_old[1], τij_old[2], τij_old[3], args.dt)
    εII   = second_invariant(ε_eff)

    args  = merge(args, (τII_old=0,))    
    τII   = compute_τII(v, εII, args)
    η_eff = 0.5*τII/εII
    τij   = 2*η_eff.*ε_eff

    return τij, τII
end