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
    εII = second_invariant(ε_eff...)

    # args = merge(args, (τII_old=0,))
    τII = first(compute_τII(v, εII, args))
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff

    return τij, τII
end

# Single material phase, staggered grid
function compute_τij(
    v, εij::NTuple{N,Union{T,NTuple{4,T}}}, args, τij_old::NTuple{N,Union{T,NTuple{4,T}}}
) where {N,T}

    # Second invariant of effective strainrate (taking elasticity into account)
    ε_eff = effective_ε(εij, v, τij_old, args.dt)
    εII = second_invariant(ε_eff...)
    ε_eff_averaged = staggered_tensor_average(ε_eff)

    # args = merge(args, (τII_old=0,))
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
    εII = second_invariant(ε_eff...)

    # args = merge(args, (τII_old=0,))
    τII = nphase(vi -> first(compute_τII(vi.CompositeRheology[1], εII, args)), phase, v)
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

    τII = nphase(vi -> first(compute_τII(vi.CompositeRheology[1], εII, args)), phases[1], v)
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff_averaged

    return τij, τII, η_eff
end

# in-place stress calculation routines

# collocated grid
function compute_τij!(
    Txx, Tyy, Txy, Tii, Txx_o, Tyy_o, Txy_o, Exx, Eyy, Exy, η_vep, P, phase, MatParam, dt
)
    Threads.@threads for j in axes(Exx, 2)
        for i in axes(Exx, 1)
            @inbounds Txx[i, j], Tyy[i, j], Txy[i, j], Tii[i, j], η_vep[i, j] = _compute_τij(
                Txx_o[i,j],
                Tyy_o[i,j],
                Txy_o[i,j],
                Exx[i,j],
                Eyy[i,j],
                Exy[i,j],
                P[i,j],
                phase[i,j],
                MatParam,
                dt,
            )
        end
    end
end

function _compute_τij(Txx_o, Tyy_o, Txy_o, Exx, Eyy, Exy, P, phase, MatParam, dt)
    args = (; dt=dt, P=P, τII_old=0.0)
    εij = (Exx, Eyy, Exy)
    τij_o = (Txx_o, Tyy_o, Txy_o)
    Tij, Tii, η_vep = compute_τij(MatParam, εij, args, τij_o, phase)

    return Tij[1], Tij[2], Tij[3], Tii, η_vep
end

# staggered grid
function compute_τij!(
    Txx, Tyy, Txy, Tii, Txx_o, Tyy_o, Txyv_o, Exx, Eyy, Exyv, η_vep, P, phase_center, phase_vertex, MatParam, dt
)
    Threads.@threads for j in axes(Exx, 2)
        for i in axes(Exx, 1)
            @inbounds Txx[i, j], Tyy[i, j], Txy[i, j], Tii[i, j], η_vep[i, j] = _compute_τij(
                Txx_o[i,j],
                Tyy_o[i,j],
                Txyv_o,
                Exx[i,j],
                Eyy[i,j],
                Exyv,
                P[i,j],
                phase_center,
                phase_vertex,
                MatParam,
                dt,
                i,
                j
            )
        end
    end
end

function _compute_τij(
    Txx_o,
    Tyy_o,
    Txyv_o,
    Exx,
    Eyy,
    Exyv,
    Pt,
    phase_center,
    phase_vertex,
    MatParam,
    dt,
    i,
    j
)

    args = (; dt=dt, P=P, τII_old=0.0)
    # gather strain rate
    εij_v = (Exyv[i, j], Exyv[i + 1, j], Exyv[i, j + 1], Exyv[i + 1, j + 1]) # gather vertices around ij center
    εij = (Exx, Eyy, εij_v)
    # gather deviatoric stress
    τij_v = (Txyv_o[i, j], Txyv_o[i + 1, j], Txyv_o[i, j + 1], Txyv_o[i + 1, j + 1]) # gather vertices around ij center
    τij_o = (Txx_o, Tyy_o, τij_v)
    # gather material phases
    phases_v = (phase_vertex[i, j], phase_vertex[i + 1, j], phase_vertex[i, j + 1], phase_vertex[i + 1, j + 1]) # gather vertices around ij center
    phases = (phase_center[i, j], phase_center[i, j], phases_v)
    # update stress and effective viscosity
    Tij, Tii, η_vep = compute_τij(MatParam, εij, args, τij_o, phase)

    return Tij[1], Tij[2], Tij[3], Tii, η_vep

end

# ----------------------------------------------------------------------------------------

## Hellper functions
@inline function staggered_tensor_average(x::NTuple{N,Union{T,NTuple{4,T}}}) where {N,T}
    ntuple(Val(N)) do i
        Base.@_inline_meta
        staggered_tensor_average(x[i])
    end
end

staggered_tensor_average(x::NTuple{N,T}) where {N,T} = sum(x) / N
staggered_tensor_average(x::T) where {T<:Number} = x