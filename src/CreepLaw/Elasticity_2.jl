using ..MaterialParameters: MaterialParamsInfo
import GeoParams.param_info

export ConstantElasticity, dεII_dτII, dτII_dεII

# ConstantElasticity  -------------------------------------------------------

"""
    ConstantElasticity(G, ν, K, E)

Structure that holds parameters for constant, isotropic, linear elasticity.
"""
@with_kw_noshow struct ConstantElasticity{T,U,U1} <: AbstractCreepLaw{T}
    G::GeoUnit{T,U} = 5e10Pa          # Elastic shear modulus
    ν::GeoUnit{T,U1} = 0.5NoUnits      # Poisson ratio
    K::GeoUnit{T,U} = 1e11Pa          # Elastic bulk modulus
    E::GeoUnit{T,U} = 1e11Pa          # Elastic Young's modulus
end
ConstantElasticity(args...) = ConstantElasticity(convert.(GeoUnit, args)...)

# Add multiple dispatch here to allow specifying combinations of 2 elastic parameters (say ν & E), to compute the others

function param_info(s::ConstantElasticity) # info about the struct
    return MaterialParamsInfo(; Equation=L"Constant elasticity")
end

# Calculation routines for linear viscous rheologies
# All inputs must be non-dimensionalized (or converted to consitent units) GeoUnits
@inline function computeCreepLaw_EpsII(
    τII::_T, s::ConstantElasticity{_T}; τII_old::_T=zero(_T), dt::_T=one(_T), kwargs...
) where {_T}
    @unpack_val G = s
    ε_el = (τII - τII_old) / (2.0 * G * dt)

    return ε_el
end

@inline function dεII_dτII(τII, s::ConstantElasticity{_T}; dt::_T=1.0, kwargs...) where {_T}
    @unpack_val G = s
    return 1 / (2 * G * dt)
end

@inline function computeCreepLaw_TauII(
    εII::_T, s::ConstantElasticity{_T}; dt::_T=one(_T), kwargs...
) where {_T}
    @unpack_val G = s
    τII = G * dt

    return τII
end

# """
#     compute_elastic_shear_strainrate(s::ConstantElasticity{_T}; τII, τII_old, dt) 

# Computes elastic strainrate given the deviatoric stress at the current (`τII`) and old timestep (`τII_old`), for a timestep `dt`:
# ```math  
#     \\dot{\\varepsilon}^{el} = {1 \\over 2 G} {D \\tau_{II} \\over Dt } ≈ {1 \\over 2 G} {\\tau_{II}- \\tilde{\\tau}_{II}^{old} \\over dt }
# ```
# Note that we here solve the scalar equation, which is sufficient for isotropic cases. In tensor form, it would be

# ```math  
#     {\\dot{\\varepsilon}^{el}}_{ij} = {1 \\over 2 G} { \\tau_{ij} - \\tilde{{\\tau_{ij}}}^{old} \\over dt }
# ```
# here ``\\tilde{{\\tau_{ij}}}^{old}`` is the rotated old deviatoric stress tensor to ensure objectivity (this can be done with Jaumann derivative, or also by using the full rotational formula).

# """

# Print info 
function show(io::IO, g::ConstantElasticity)
    return print(
        io,
        "Linear elasticity with shear modulus: G = $(UnitValue(g.G)), Poison's ratio: ν = $(UnitValue(g.ν)), bulk modulus: K = $(UnitValue(g.K)) and Young's module: E=$(UnitValue(g.E))",
    )
end
