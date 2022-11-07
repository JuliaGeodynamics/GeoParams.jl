# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

abstract type AbstractElasticity{T} <: AbstractConstitutiveLaw{T} end

export compute_εII,            # calculation routines
    compute_εII!,
    compute_τII,
    compute_τII!,
    dεII_dτII,
    dτII_dεII,
    compute_εvol,
    compute_εvol!,
    compute_p,
    compute_p!,
    dεvol_dp,
    dp_dεvol,
    param_info,
    ConstantElasticity,     # constant
    SetConstantElasticity,  # helper function
    AbstractElasticity,
    isvolumetric, 
    effective_εII, effective_ε

# ConstantElasticity  -------------------------------------------------------

"""
    ConstantElasticity(G, ν, K, E)

Structure that holds parameters for constant, isotropic, linear elasticity.
"""
@with_kw_noshow struct ConstantElasticity{T,U,U1} <: AbstractElasticity{T}
    G::GeoUnit{T,U} = 5e10Pa                                             # Elastic shear modulus
    ν::GeoUnit{T,U1} = 0.5NoUnits                                         # Poisson ratio
    Kb::GeoUnit{T,U} = 2 * G * (1 + ν) / (3 * (1 - 2 * ν))                          # Elastic bulk modulus
    E::GeoUnit{T,U} = 9 * Kb * G / (3 * Kb + G)                                  # Elastic Young's modulus
end

ConstantElasticity(args...) = ConstantElasticity(convert.(GeoUnit, args)...)


# Add multiple dispatch here to allow specifying combinations of 2 elastic parameters (say ν & E), to compute the others
"""
    SetConstantElasticity(; G=nothing, ν=nothing, E=nothing, Kb=nothing)

This allows setting elastic parameters by specifying 2 out of the 4 elastic parameters `G` (Elastic shear modulus), `ν` (Poisson's ratio), `E` (Young's modulus), or `Kb` (bulk modulus).
"""
function SetConstantElasticity(; G=nothing, ν=nothing, E=nothing, Kb=nothing)
    if (!isnothing(G) && !isnothing(ν))
        Kb = 2 * G * (1 + ν) / (3 * (1 - 2 * ν))     # Bulk modulus
        E = 9 * Kb * G / (3 * Kb + G)              # Youngs modulus
    elseif (!isnothing(Kb) && !isnothing(ν))
        G = (3 * Kb * (1 - 2 * ν)) / (2 * (1 + ν))
        E = 9 * Kb * G / (3 * Kb + G)              # Youngs modulus
    elseif (!isnothing(E) && !isnothing(ν))
        Kb = E / (3 * (1 - 2 * ν))              # Bulk modulus
        G = E / (2 * (1 + ν))                # Shear modulus
    elseif (!isnothing(Kb) && !isnothing(G))
        E = 9 * Kb * G / (3 * Kb + G)
        ν = (3 * Kb - 2 * G) / (2 * (3 * Kb + G))  # Poisson's ratio
    end

    if !isa(G, Quantity)
        G = G * Pa
    end
    if !isa(ν, Quantity)
        ν = ν * NoUnits
    end
    if !isa(E, Quantity)
        E = E * Pa
    end
    if !isa(Kb, Quantity)
        Kb = Kb * Pa
    end

    return ConstantElasticity(GeoUnit(G), GeoUnit(ν), GeoUnit(Kb), GeoUnit(E))
end

function param_info(s::ConstantElasticity) # info about the struct
    return MaterialParamsInfo(; Equation=L"Constant elasticity")
end

function isvolumetric(a::ConstantElasticity)
    @unpack_val ν = a
    return ν == 0.5 ? false : true
end

# Calculation routines
"""
    compute_εII(s::ConstantElasticity{_T}, τII; τII_old, dt) 

Computes elastic strainrate given the deviatoric stress at the current (`τII`) and old timestep (`τII_old`), for a timestep `dt`:
```math  
    \\dot{\\varepsilon}^{el} = {1 \\over 2 G} {D \\tau_{II} \\over Dt } ≈ {1 \\over 2 G} {\\tau_{II}- \\tilde{\\tau}_{II}^{old} \\over dt }
```
Note that we here solve the scalar equation, which is sufficient for isotropic cases. In tensor form, it would be

```math  
    {\\dot{\\varepsilon}^{el}}_{ij} = {1 \\over 2 G} { \\tau_{ij} - \\tilde{{\\tau_{ij}}}^{old} \\over dt }
```
here ``\\tilde{{\\tau_{ij}}}^{old}`` is the rotated old deviatoric stress tensor to ensure objectivity (this can be done with Jaumann derivative, or also by using the full rotational formula).

"""
@inline function compute_εII(
    a::ConstantElasticity, τII::_T; τII_old=zero(precision(a)), dt=one(precision(a)), kwargs...
) where {_T}
    @unpack_val G = a
    ε_el = 0.5 * (τII - τII_old) / (G * dt)
    return ε_el
end

@inline function dεII_dτII(a::ConstantElasticity{_T}, τII::_T; τII_old=zero(precision(a)), dt=one(precision(a)), kwargs...
    ) where {_T}
    @unpack_val G = a
    return 0.5 * inv(G * dt)
end

@inline function compute_τII(
    a::ConstantElasticity, εII::_T; τII_old=zero(precision(a)), dt=one(precision(a)), kwargs...
) where {_T}
    @unpack_val G = a
    τII = _T(2) * G * dt * εII + τII_old

    return τII
end

@inline function dτII_dεII(a::ConstantElasticity{_T}, τII_old=zero(precision(a)), dt=one(precision(a)), kwargs...
    ) where {_T}
    @unpack_val G = a
    return _T(2) * G * dt
end

"""
    compute_εII!(ε_el::AbstractArray{_T,N}, s::ConstantElasticity{_T}; τII::AbstractArray{_T,N}, τII_old::AbstractArray{_T,N}, dt::_T, kwargs...) 

In-place computation of the elastic shear strainrate for given deviatoric stress invariants at the previous (`τII_old`) and new (`τII`) timestep, as well as the timestep `dt`  

```math  
    \\dot{\\varepsilon}^{el} = {1 \\over 2 G} {D \\tau_{II} \\over Dt } ≈ {1 \\over 2 G} {\\tau_{II}- \\tau_{II}^{old} \\over dt }
```

"""
function compute_εII!(
    ε_el::AbstractArray{_T,N},
    p::ConstantElasticity{_T},
    τII::AbstractArray{_T,N};
    τII_old::AbstractArray{_T,N},
    dt::_T,
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(τII)
        ε_el[i] = compute_εII(p, τII[i]; τII_old=τII_old[i], dt=dt)
    end
    return nothing
end

"""
    compute_τII!(τII::AbstractArray{_T,N}, s::ConstantElasticity{_T}. ε_el::AbstractArray{_T,N}; τII_old::AbstractArray{_T,N}, dt::_T, kwargs...) 

In-place update of the elastic stress for given deviatoric strainrate invariants and stres invariant at the old (`τII_old`) timestep, as well as the timestep `dt`  

```math  
    \\tau_{II} = 2 G dt \\dot{\\varepsilon}^{el} + \\tau_{II}^{old}
```

"""
function compute_τII!(
    τII::AbstractArray{_T,N},
    p::ConstantElasticity{_T},
    ε_el::AbstractArray{_T,N};
    τII_old::AbstractArray{_T,N},
    dt::_T,
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(ε_el)
        τII[i] = compute_τII(p, ε_el[i]; τII_old=τII_old[i], dt=dt)
    end
    return nothing
end

# Print info 
function show(io::IO, g::ConstantElasticity)
    return print(
        io,
        "Linear elasticity with shear modulus: G = $(UnitValue(g.G)), Poisson's ratio: ν = $(UnitValue(g.ν)), bulk modulus: Kb = $(UnitValue(g.Kb)) and Young's module: E=$(UnitValue(g.E))",
    )
end

"""
    compute_εvol(s::ConstantElasticity{_T}, P; P_old, dt) 

Computes elastic volumetric strainrate given the pressure at the current (`P`) and old timestep (`P_old`), for a timestep `dt`:
```math  
    \\dot{\\vartheta}^{el} = {1 \\over Kb} {D P \\over Dt } ≈ {1 \\over Kb} {P - \\tilde{P^{old} \\over dt }
```

"""
@inline function compute_εvol(
    a::ConstantElasticity, P::_T; P_old=zero(precision(a)), dt=one(precision(a)), kwargs...
) where {_T}
    @unpack_val Kb = a
    εvol_el = - (P - P_old) / (Kb * dt)
    return εvol_el
end

@inline function dεvol_dp(a::ConstantElasticity{_T}, P::_T; P_old=zero(precision(a)), dt=one(precision(a)), kwargs...
    ) where {_T}
    @unpack_val Kb = a
    return - inv(Kb * dt)
end

@inline function compute_p(
    a::ConstantElasticity, εvol::_T; P_old=zero(precision(a)), dt=one(precision(a)), kwargs...
) where {_T}
    @unpack_val Kb = a
    P = - Kb * dt * εvol + P_old

    return P
end

@inline function dp_dεvol(a::ConstantElasticity{_T}, P_old=zero(precision(a)), dt=one(precision(a)), kwargs...
    ) where {_T}
    @unpack_val Kb = a
    return -Kb * dt
end

"""
    compute_εvol!(s::ConstantElasticity{_T}, P; P_old, dt) 

    In-place computation of the elastic volumetric strainrate given the pressure at the current (`P`) and old timestep (`P_old`), for a timestep `dt`:
```math  
    \\dot{\\vartheta}^{el} = {1 \\over Kb} {D P \\over Dt } ≈ {1 \\over Kb} {P - \\tilde{P^{old} \\over dt }
```

"""
function compute_εvol!(
    εvol_el::AbstractArray{_T,N},
    a::ConstantElasticity{_T},
    P::AbstractArray{_T,N};
    P_old::AbstractArray{_T,N},
    dt::_T,
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(P)
        εvol_el[i] = compute_εvol(a, P[i]; P_old=P_old[i], dt=dt)
    end
    return nothing
end

"""
    compute_p!(p::AbstractArray{_T,N}, s::ConstantElasticity{_T}. εvol_el::AbstractArray{_T,N}; P_old::AbstractArray{_T,N}, dt::_T, kwargs...) 

In-place update of the elastic pressure for given volumetric strainrate and pressure at the old (`P_old`) timestep, as well as the timestep `dt`  

```math  
    \\p = Kb dt \\dot{\\vartheta}^{el} + \\p^{old}
```

"""
function compute_p!(
    P::AbstractArray{_T,N},
    a::ConstantElasticity{_T},
    εvol_el::AbstractArray{_T,N};
    P_old::AbstractArray{_T,N},
    dt::_T,
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(εvol_el)
        P[i] = compute_p(a, εvol_el[i]; P_old=P_old[i], dt=dt)
    end
    return nothing
end

#-------------------------------------------------------------------------

# Computational routines needed for computations with the MaterialParams structure 
function compute_εII(s::AbstractMaterialParamsStruct, args)
    if isempty(s.Elasticity)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return s.Elasticity[1](args)
    end
end

function compute_εvol(s::AbstractMaterialParamsStruct, args)
    println("hllo")
    if isempty(s.Elasticity)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return s.Elasticity[1](args)
    end
end

## Kernels for effective strain rate (Eij_eff = Eij + Tij/(2 G dt) ) 

@inline _elastic_ε(v::ConstantElasticity, τij_old, dt) = τij_old / (2 * v.G * dt)
@inline _elastic_ε(v::Vararg{Any, N}) where {N} = 0.0

@inline effective_ε(εij::T, v, τij_old::T, dt) where {T} = εij + elastic_ε(v, τij_old, dt)
@inline effective_ε(εij::T, v, τij_old::T, dt, phase::Int64) where {T} = εij + elastic_ε(v, τij_old, dt, phase)

# By letting ::NTuple{N, Any} it can recursively call itself when 
# strain and stress are nested NamedTuples i.e. in staggered grids
# NOTE: ::NTuple{N, Union{T,Any}} makes type unstable (need fix?)
@inline function effective_ε(
    εij::NTuple{N, Any}, v, τij_old::NTuple{N, Any}, dt
) where N
    return ntuple(i -> effective_ε(εij[i], v, τij_old[i], dt), Val(N))
end

@inline function effective_ε(
    εij::NTuple{N, Any}, v, τij_old::NTuple{N, Any}, dt, phase::Int64
) where N
    return ntuple(i -> effective_ε(εij[i], v, τij_old[i], dt, phase), Val(N))
end

@inline function effective_ε(
    εij::NTuple{N, T}, v, τij_old::NTuple{N,T}, dt
) where {N,T}
    return ntuple(i -> effective_ε(εij[i], v, τij_old[i], dt), Val(N))
end

@inline function effective_ε(
    εij::NTuple{N, T}, v, τij_old::NTuple{N,T}, dt, phase::Int64
) where {N,T}
    return ntuple(i -> effective_ε(εij[i], v, τij_old[i], dt, phase), Val(N))
end

# 2D wrapper 
function effective_ε(εxx, εyy, εxy, v, τxx_old, τyy_old, τxy_old, dt)
    return effective_ε((εxx, εyy, εxy), v, (τxx_old, τyy_old, τxy_old), dt)
end

function effective_ε(εxx, εyy, εxy, v, τxx_old, τyy_old, τxy_old, dt, phase::Int64)
    return effective_ε((εxx, εyy, εxy), v, (τxx_old, τyy_old, τxy_old), dt, phase)
end

function effective_εII(εxx, εyy, εxy, v, τxx_old, τyy_old, τxy_old, dt)
    εxx, εyy, εxy = effective_ε(εxx, εyy, εxy, v, τxx_old, τyy_old, τxy_old, dt)
    εII = second_invariant(εxx, εyy, εxy)
    return εII
end

function effective_εII(εxx, εyy, εxy, v, τxx_old, τyy_old, τxy_old, dt, phase::Int64)
    εxx, εyy, εxy = effective_ε(εxx, εyy, εxy, v, τxx_old, τyy_old, τxy_old, dt, phase)
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


#=
# add methods programmatically
for myType in (:ConstantElasticity,)
    @eval begin
        (s::$(myType))(args)= s(; args...)
        compute_εII(s::$(myType), args) = s(args)
        compute_εII!(ε_el::AbstractArray{_T,N}, s::$(myType){_T}, args) where {_T,N} = compute_εII!(ε_el, s; args...)
        dεII_dτII(s::ConstantElasticity, args) = dεII_dτII(s; args...)
    end
end
=#

#compute_εII(args...) = compute_param(compute_εII, args...)
#compute_εII!(args...) = compute_param!(compute_εII, args...)
#compute_εvol(args...) = compute_param(compute_εvol, args...)
#3ompute_εvol!(args...) = compute_param!(compute_εvol, args...)
#compute_p(args...) = compute_param(compute_p, args...)
#compute_p!(args...) = compute_param!(compute_p, args...)
