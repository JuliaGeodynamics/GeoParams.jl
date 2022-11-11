# Stress tensor computations
using StaticArrays
export compute_τij,
    compute_p_τij,
    compute_τij!,
    compute_p_τij!,
    compute_τij_stagcenter!,
    compute_p_τij_stagcenter!

"""
    τij,τII = compute_τij(v, εij::NTuple{n,T}, args, τij_old::NTuple{N,T})

Computes deviatoric stress `τij` for given deviatoric strain rate `εij`, old stress `τij_old`, rheology `v` and arguments `args`.
This is for a case that all points are collocated and we have a single phase. 
"""
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

"""
    p,τij,τII = compute_p_τij(v, εij::NTuple{n,T}, P_old::T, args,  τij_old::NTuple{N,T})

Computes pressure `p` and deviatoric stress `τij` for given deviatoric strain rate `εij`, old stress `τij_old`, rheology `v` and arguments `args`.
This is for a case that all points are collocated and we have a single phase. 
"""
function compute_p_τij(
    v, εij::NTuple{N,T}, P_old::T, args, τij_old::NTuple{N,T}
) where {T,N}

    # Second invariant of effective strainrate (taking elasticity into account)
    return ε_eff = effective_ε(εij, v, τij_old, args.dt)
    # εII = second_invariant(ε_eff)
    # εvol = volumetric_strainrate(εij)    # Volumetric strainrate

    # args = merge(args, (τII_old=0,P_old=P_old))
    # P,τII = compute_p_τII(v, εII, εvol, args)
    # η_eff = 0.5 * τII / εII
    # τij = 2 * η_eff .* ε_eff

    # return P, τij, τII, η_eff
end

"""
    τij, τII = compute_τij(v, εij::NTuple{N,Union{T,NTuple{4,T}}}, args, τij_old::NTuple{3,T})

Computes deviatoric stress `τij` for given deviatoric strain rate `εij`, old stress `τij_old`, rheology `v` and arguments `args`.
This is for a staggered grid case with a single phase. 
"""
function compute_τij(
    v, εij::NTuple{N,Union{T,NTuple{4,T}}}, args, τij_old::NTuple{N,Union{T,NTuple{4,T}}}
) where {N,T}

    # Second invariant of effective strainrate (taking elasticity into account)
    ε_eff = effective_ε(εij, v, τij_old, args.dt)
    εII = second_invariant_staggered(ε_eff...)
    ε_eff_averaged = staggered_tensor_average(ε_eff)

    # args = merge(args, (τII_old=0,))
    τII = first(compute_τII(v, εII, args))
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff_averaged

    return τij, τII, η_eff
end

"""
    p, τij, τII = compute_p_τij(v, εij::NTuple{N,Union{T,NTuple{4,T}}}, P_old::T, args, τij_old::NTuple{3,T})

Computes deviatoric stress `τij` for given deviatoric strain rate `εij`, old stress `τij_old`, rheology `v` and arguments `args`.
This is for a staggered grid case with a single phase. 
"""
function compute_p_τij(
    v,
    εij::NTuple{N,Union{T,NTuple{4,T}}},
    P_old::T,
    args,
    τij_old::NTuple{N,Union{T,NTuple{4,T}}},
) where {N,T}

    # Second invariant of effective strainrate (taking elasticity into account)
    ε_eff = effective_ε(εij, v, τij_old, args.dt)
    εII = second_invariant_staggered(ε_eff...)
    ε_eff_averaged = staggered_tensor_average(ε_eff)
    εvol = volumetric_strainrate(staggered_tensor_average(εij))

    args = merge(args, (P_old=P_old, τII_old=0.0))
    P, τII, = compute_p_τII(v, εII, εvol, args)
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff_averaged

    return P, τij, τII, η_eff
end

# Multiple material phases, collocated grid
"""
    τij, τII = compute_τij(v::NTuple{N1,AbstractMaterialParamsStruct}, εij::NTuple{N2,T}, args, τij_old::NTuple{3,T}, phase::I)

Computes deviatoric stress `τij` for given deviatoric strain rate `εij`, old stress `τij_old`, rheology `v` and arguments `args`.
This is for a collocated grid case with a single phase `phase`. 
"""
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

"""
    P,τij, τII = compute_p_τij(v::NTuple{N1,AbstractMaterialParamsStruct}, εij::NTuple{N2,T}, P_old::T, args, τij_old::NTuple{3,T}, phase::I)

Computes deviatoric stress `τij` for given deviatoric strain rate `εij`, old stress `τij_old`, rheology `v` and arguments `args`.
This is for a collocated grid case with a single phase `phase`. 
"""
function compute_p_τij(
    v::NTuple{N1,AbstractMaterialParamsStruct},
    εij::NTuple{N2,T},
    P_old::T,
    args,
    τij_old::NTuple{N2,T},
    phase::I,
) where {T,N1,N2,I<:Integer}

    # Second invariant of effective strainrate (taking elasticity into account)
    ε_eff = effective_ε(εij, v, τij_old, args.dt, phase)
    εII = second_invariant(ε_eff)
    εvol = volumetric_strainrate(εij)

    args = merge(args, (τII_old=0, P_old=P_old))
    P, τII = nphase(vi -> compute_p_τII(vi.CompositeRheology[1], εII, εvol, args), phase, v)
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff

    return P, τij, τII, η_eff
end

# Multiple material phases, staggered grid
"""
    τij, τII, η_eff = compute_τij(v::NTuple{N1,AbstractMaterialParamsStruct}, εij::NTuple, args, τij_old::NTuple, phases::NTuple)

This computes deviatoric stress components `τij`, their second invariant `τII`, and effective viscosity `η_eff` for given deviatoric strainrates `εij`, old stresses `τij_old`, `phases` (integer) for every point and arguments `args`.
This handles various staggered grid arrangements; if staggered components are given as `NTuple{4,T}`, they will be averaged. Note that the phase of all staggered points should be the same.
"""
function compute_τij(
    v::NTuple{N1,AbstractMaterialParamsStruct},
    εij::NTuple{N2,Union{T,NTuple{4,T}}},
    args,
    τij_old::NTuple{N2,Union{T,NTuple{4,T}}},
    phases::NTuple{N2,Union{I,NTuple{4,I}}},
) where {T,N1,N2,I<:Integer}

    # Second invariant of effective strainrate (taking elasticity into account)
    ε_eff = effective_ε(εij, v, τij_old, args.dt, phases)
    εII = second_invariant_staggered(ε_eff...)
    ε_eff_averaged = staggered_tensor_average(ε_eff)

    τII = nphase(
        vi -> first(compute_τII(vi.CompositeRheology[1], εII, args)), phases[1][1], v
    )
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff_averaged

    return τij, τII, η_eff
end

"""
    τij, τII, η_eff = compute_p_τij(v::NTuple{N1,AbstractMaterialParamsStruct}, εij::NTuple, args, τij_old::NTuple, phases::NTuple)

This computes pressure `p` and deviatoric stress components `τij`, their second invariant `τII`, and effective viscosity `η_eff` for given deviatoric strainrates `εij`, old stresses `τij_old`, `phases` (integer) for every point and arguments `args`.
It handles staggered grids of various taste
"""
function compute_p_τij(
    v::NTuple{N1,AbstractMaterialParamsStruct},
    εij::NTuple{N2,Union{T,NTuple{4,T}}},
    P_old::T,
    args,
    τij_old::NTuple{N2,Union{T,NTuple{4,T}}},
    phases::NTuple{N2,Union{I,NTuple{4,I}}},
) where {T,N1,N2,I<:Integer}

    # Second invariant of effective strainrate (taking elasticity into account)
    ε_eff = effective_ε(εij, v, τij_old, args.dt, phases)
    εII = second_invariant_staggered(ε_eff...)
    ε_eff_averaged = staggered_tensor_average(ε_eff)
    εvol = volumetric_strainrate(staggered_tensor_average(εij))

    args = merge(args, (τII_old=0, P_old=P_old))
    P, τII = nphase(
        vi -> compute_p_τII(vi.CompositeRheology[1], εII, εvol, args), phases[1][1], v
    )  # note: assumes phases of all staggered points to be the same! 
    η_eff = 0.5 * τII / εII
    τij = 2 * η_eff .* ε_eff_averaged

    return P, τij, τII, η_eff
end

# in-place stress calculation routines

# collocated grid
"""
    compute_τij!(Txx, Tyy, Txy, Tii, η_vep, Exx, Eyy, Exy, P, Txx_o, Tyy_o, Txy_o, phase, MatParam, dt)

Computes 2D deviatoric stress components `(Txx,Tyy,Txy)` given deviatoric strainrate components `(Exx,Eyy,Exy)` and old deviatoric stresses `(Txx_o, Tyy_o, Txy_o)` (only used for viscoelastic cases).
Also returned are `Tii` (second invariant of the deviatoric stress tensor), and `η_vep` the viscoelastoplastic effective viscosity. 
Also required as input is `MatParam`, the material parameters for every phase and `phase`, an integer array of `size(Exx)` that indicates the phase of every point.

This function assumes that strainrate points are collocated and that `Exx`,`Eyy`,`Exy` are at the same points.
# """

function compute_τij!(
    Txx::AbstractArray{_T},
    Tyy::AbstractArray{_T},
    Txy::AbstractArray{_T},
    Tii::AbstractArray{_T},
    η_vep::AbstractArray{_T},
    Exx::AbstractArray{_T},
    Eyy::AbstractArray{_T},
    Exy::AbstractArray{_T},
    P::AbstractArray{_T},
    Txx_o::AbstractArray{_T},
    Tyy_o::AbstractArray{_T},
    Txy_o::AbstractArray{_T},
    phase::AbstractArray{I},
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    dt::_T,
) where {_T<:Number,I<:Integer,N}
    Threads.@threads for i in eachindex(Exx)
        @inbounds Txx[i], Tyy[i], Txy[i], Tii[i], η_vep[i] = _compute_τij_2D(
            Exx[i],
            Eyy[i],
            Exy[i],
            P[i],
            Txx_o[i],
            Tyy_o[i],
            Txy_o[i],
            phase[i],
            MatParam,
            dt,
        )
    end
end

# # Internal computation array
function _compute_τij_2D(Exx, Eyy, Exy, P, Txx_o, Tyy_o, Txy_o, phase, MatParam, dt)
    args = (; dt=dt, P=P, τII_old=0.0)
    εij = (Exx, Eyy, Exy)
    τij_o = (Txx_o, Tyy_o, Txy_o)
    Tij, Tii, η_vep = compute_τij(MatParam, εij, args, τij_o, phase)

    return Tij[1], Tij[2], Tij[3], Tii, η_vep
end

"""
    compute_p_τij!(Txx, Tyy, Txy, Tii, η_vep, P, Exx, Eyy, Exy, P_o, Txx_o, Tyy_o, Txy_o, phase, MatParam, dt)

Computes 2D deviatoric stress components `(Txx,Tyy,Txy)` and pressure `P` given deviatoric strainrate components `(Exx,Eyy,Exy)` and old deviatoric stresses `(Txx_o, Tyy_o, Txy_o)` and old pressure `P_o` (only used for viscoelastic cases).
Also returned are `Tii` (second invariant of the deviatoric stress tensor), and `η_vep` the viscoelastoplastic effective viscosity. 
Also required as input is `MatParam`, the material parameters for every phase and `phase`, an integer array of `size(Exx)` that indicates the phase of every point.

This function assumes that strainrate points are collocated and that `Exx`,`Eyy`,`Exy` are at the same points.
"""
function compute_p_τij!(
    Txx::AbstractArray{_T},
    Tyy::AbstractArray{_T},
    Txy::AbstractArray{_T},
    Tii::AbstractArray{_T},
    η_vep::AbstractArray{_T},
    P::AbstractArray{_T},
    Exx::AbstractArray{_T},
    Eyy::AbstractArray{_T},
    Exy::AbstractArray{_T},
    P_o::AbstractArray{_T},
    Txx_o::AbstractArray{_T},
    Tyy_o::AbstractArray{_T},
    Txy_o::AbstractArray{_T},
    phase::AbstractArray{I},
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    dt::_T,
) where {_T<:Number,I<:Integer,N}
    Threads.@threads for i in eachindex(Exx)
        @inbounds P[i], Txx[i], Tyy[i], Txy[i], Tii[i], η_vep[i] = _compute_p_τij_2D(
            Exx[i],
            Eyy[i],
            Exy[i],
            P[i],
            Txx_o[i],
            Tyy_o[i],
            Txy_o[i],
            P_o[i],
            phase[i],
            MatParam,
            dt,
        )
    end
end

# Internal computation array
function _compute_p_τij_2D(Exx, Eyy, Exy, P, Txx_o, Tyy_o, Txy_o, P_o, phase, MatParam, dt)
    args = (; dt=dt, P=P, τII_old=0.0)
    εij = (Exx, Eyy, Exy)
    τij_o = (Txx_o, Tyy_o, Txy_o)
    P, Tij, Tii, η_vep = compute_p_τij(MatParam, εij, P_o, args, τij_o, phase)

    return P, Tij[1], Tij[2], Tij[3], Tii, η_vep
end

"""
    compute_τij_stagcenter!(Txx, Tyy, Txy, Tii, η_vep,  Exx, Eyy, Exyv, P, Txx_o, Tyy_o, Txyv_o, phase_center, phase_vertex, MatParam, dt)

Updates deviatoric stresses on a staggered grid in a centered based manner (averages strainrates/old stresses from vertices -> centers).
Take care of the sizes of the input matrixes:
- `Txx`,`Tyy`,`Txy`,`Tii`: 2D matrixes of size (nx,ny) which describe updated deviatoric stress components (x,y,xy, 2nd invariant) at the center point
- `η_vep`: viscoelastoplastic viscosity @ center (nx,ny)
- `Exx`,`Eyy`,`P`: deviatoric strain rate components & pressure @ center 
- `Exy`: shear strain rate @ vertices (nx+1,ny+1)
- `Txx_o`, `Tyy_o`: deviatoric stress normal components of last timestep at center (nx,ny)
- `Txy_o`: deviatoric stress shear components at vertices (nx+1,ny+1)
- `phase_center`: integer array with phase at center points (nx,ny)
- `phase_vertex`: integer array with phase at vertex points (nx+1,ny+1)
- `MatParam`: material parameter array that includes a `CompositeRheology` field which specifies the rheology for different phases
- `dt`: timestep 
"""
function compute_τij_stagcenter!(
    Txx::AbstractArray{_T},
    Tyy::AbstractArray{_T},
    Txy::AbstractArray{_T},
    Tii::AbstractArray{_T},
    η_vep::AbstractArray{_T},
    Exx::AbstractArray{_T},
    Eyy::AbstractArray{_T},
    Exyv::AbstractArray{_T},
    P::AbstractArray{_T},
    Txx_o::AbstractArray{_T},
    Tyy_o::AbstractArray{_T},
    Txyv_o::AbstractArray{_T},
    phase_center::AbstractArray{I},
    phase_vertex::AbstractArray{I},
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    dt::_T,
) where {_T<:Number,I<:Integer,N}
    Threads.@threads for j in axes(Exx,2)
        for i in axes(Exx,1)
            @inbounds Txx[i,j], Tyy[i,j], Txy[i,j], Tii[i,j], η_vep[i,j] = _compute_τij_stagcenter(
                Exx[i,j],
                Eyy[i,j],
                Exyv,
                P[i,j],
                Txx_o[i,j],
                Tyy_o[i,j],
                Txyv_o,
                phase_center,
                phase_vertex,
                MatParam,
                dt,
                i,
                j,
            )
        end
    end
end

function _compute_τij_stagcenter(
    Exx, Eyy, Exyv, Pt, Txx_o, Tyy_o, Txyv_o, phase_center, phase_vertex, MatParam, dt, i, j
)
    args = (; dt=dt, P=Pt, τII_old=0.0)
    # gather strain rate
    εij_v = (Exyv[i, j], Exyv[i + 1, j], Exyv[i, j + 1], Exyv[i + 1, j + 1]) # gather vertices around ij center
    εij = (Exx, Eyy, εij_v)
    # gather old deviatoric stress
    τij_v = (Txyv_o[i, j], Txyv_o[i + 1, j], Txyv_o[i, j + 1], Txyv_o[i + 1, j + 1]) # gather vertices around ij center
    τij_o = (Txx_o, Tyy_o, τij_v)
    # gather material phases
    phases_v = (
        phase_vertex[i, j],
        phase_vertex[i + 1, j],
        phase_vertex[i, j + 1],
        phase_vertex[i + 1, j + 1],
    ) # gather vertices around ij center
    phases = (phase_center[i, j], phase_center[i, j], phases_v)

    # update stress and effective viscosity
    Tij, Tii, η_vep = compute_τij(MatParam, εij, args, τij_o, phases)

    return Tij[1], Tij[2], Tij[3], Tii, η_vep
end

"""
    compute_τij_stagcenter!(Txx, Tyy, Txy, Tii, η_vep,  Exx, Eyy, Exyv, P, Txx_o, Tyy_o, Txyv_o, phase_center, phase_vertex, MatParam, dt)

Updates deviatoric stresses on a staggered grid in a centered based manner (averages strainrates/old stresses from vertices -> centers).
Take care of the sizes of the input matrixes:
- `Txx`,`Tyy`,`Txy`,`Tii`: 2D matrixes of size (nx,ny) which describe updated deviatoric stress components (x,y,xy, 2nd invariant) at the center point
- `η_vep`: viscoelastoplastic viscosity @ center (nx,ny)
- `Exx`,`Eyy`,`P`: deviatoric strain rate components & pressure @ center 
- `Exy`: shear strain rate @ vertices (nx+1,ny+1)
- `Txx_o`, `Tyy_o`: deviatoric stress normal components of last timestep at center (nx,ny)
- `Txy_o`: deviatoric stress shear components at vertices (nx+1,ny+1)
- `phase_center`: integer array with phase at center points (nx,ny)
- `phase_vertex`: integer array with phase at vertex points (nx+1,ny+1)
- `MatParam`: material parameter array that includes a `CompositeRheology` field which specifies the rheology for different phases
- `dt`: timestep 
"""
function compute_p_τij_stagcenter!(
    Txx::AbstractArray{_T},
    Tyy::AbstractArray{_T},
    Txy::AbstractArray{_T},
    Tii::AbstractArray{_T},
    η_vep::AbstractArray{_T},
    P::AbstractArray{_T},
    Exx::AbstractArray{_T},
    Eyy::AbstractArray{_T},
    Exyv::AbstractArray{_T},
    P_o::AbstractArray{_T},
    Txx_o::AbstractArray{_T},
    Tyy_o::AbstractArray{_T},
    Txyv_o::AbstractArray{_T},
    phase_center::AbstractArray{I},
    phase_vertex::AbstractArray{I},
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    dt::_T,
) where {_T<:Number,I<:Integer,N}
    Threads.@threads for j in axes(Exx,2)
        for i in axes(Exx,1)
            @inbounds Txx[i,j], Tyy[i,j], Txy[i,j], Tii[i,j], η_vep[i,j]  = _compute_p_τij_stagcenter(
                Exx[i,j],
                Eyy[i,j],
                Exyv,
                P[i,j],
                Txx_o[i,j],
                Tyy_o[i,j],
                Txyv_o,
                P_o[i,j],
                phase_center,
                phase_vertex,
                MatParam,
                dt,
                i,
                j,
            )
        end
    end
end

function _compute_p_τij_stagcenter(
    Exx,
    Eyy,
    Exyv,
    P,
    Txx_o,
    Tyy_o,
    Txyv_o,
    P_o,
    phase_center,
    phase_vertex,
    MatParam,
    dt,
    i,
    j,
)
    args = (; dt=dt, P=P, τII_old=0.0)
    # gather strain rate
    εij_v = (Exyv[i, j], Exyv[i + 1, j], Exyv[i, j + 1], Exyv[i + 1, j + 1]) # gather vertices around ij center
    εij = (Exx, Eyy, εij_v)
    # gather old deviatoric stress
    τij_v = (Txyv_o[i, j], Txyv_o[i + 1, j], Txyv_o[i, j + 1], Txyv_o[i + 1, j + 1]) # gather vertices around ij center
    τij_o = (Txx_o, Tyy_o, τij_v)
    # gather material phases
    phases_v = (
        phase_vertex[i, j],
        phase_vertex[i + 1, j],
        phase_vertex[i, j + 1],
        phase_vertex[i + 1, j + 1],
    ) # gather vertices around ij center
    phases = (phase_center[i, j], phase_center[i, j], phases_v)

    # update stress and effective viscosity
    P, Tij, Tii, η_vep = compute_p_τij(MatParam, εij, P_o, args, τij_o, phases)

    return P, Tij[1], Tij[2], Tij[3], Tii, η_vep
end

# ----------------------------------------------------------------------------------------

## Helper functions
@inline function staggered_tensor_average(x::NTuple{N,Union{T,NTuple{4,T}}}) where {N,T}
    ntuple(Val(N)) do i
        Base.@_inline_meta
        staggered_tensor_average(x[i])
    end
end

staggered_tensor_average(x::NTuple{N,T}) where {N,T} = sum(x) / N
staggered_tensor_average(x::T) where {T<:Number} = x

@inline function volumetric_strainrate(x::NTuple{3,Union{T,NTuple{4,T}}}) where {T}
    return vol = x[1] + x[2]  #2D
end

@inline function volumetric_strainrate(x::NTuple{6,Union{T,NTuple{4,T}}}) where {T}
    return vol = x[1] + x[2] + x[3] #3D
end
