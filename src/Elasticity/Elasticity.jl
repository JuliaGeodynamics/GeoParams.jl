# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

abstract type AbstractElasticity{T} <: AbstractConstitutiveLaw{T} end

export compute_εII,            # calculation routines
    compute_εII!,
    compute_τII,
    compute_τII!,
    dεII_dτII,
    dτII_dεII,
    param_info,
    ConstantElasticity,     # constant
    SetConstantElasticity,  # helper function
    AbstractElasticity

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
function compute_εII(
    s::ConstantElasticity, τII::_T; τII_old::_T=zero(_T), dt::_T=1.0, kwargs...
) where {_T}
    @unpack_val G = s
    ε_el = 0.5 * (τII - τII_old) / (G * dt)
    return ε_el
end

function dεII_dτII(s::ConstantElasticity{_T}, τII::_T; dt::_T=1.0, kwargs...) where {_T}
    @unpack_val G = s
    return 0.5 * (G * dt)
end

function compute_τII(
    s::ConstantElasticity, εII::_T; τII_old::_T=zero(_T), dt::_T=1.0, kwargs...
) where {_T}
    @unpack_val G = s
    τII = 2.0 * G * dt * εII + τII_old
    return τII
end

function dτII_dεII(s::ConstantElasticity{_T}, εII::_T; dt::_T=1.0, kwargs...) where {_T}
    @unpack_val G = s
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
#-------------------------------------------------------------------------

# Computational routines needed for computations with the MaterialParams structure 
function compute_εII(s::AbstractMaterialParamsStruct, args)
    if isempty(s.Elasticity)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return s.Elasticity[1](args)
    end
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

compute_εII(args...) = compute_param(compute_εII, args...)
compute_εII!(args...) = compute_param!(compute_εII, args...)
