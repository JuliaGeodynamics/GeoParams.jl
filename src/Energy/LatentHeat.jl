module LatentHeat

# This implements latent heat. There are two options:
# 1) Constant latent heat as a source term to the energy equation (usually numerically unstable)
# 2) Latent heat by modifying heat capacity (usually more stable)
# Note that 1) is implemented in this module, but that 2) is added to the HeatCapacity module

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
using ..MaterialParameters: MaterialParamsInfo
using ..HeatCapacity: AbstractHeatCapacity, ConstantHeatCapacity

import Base.show, GeoParams.param_info

abstract type AbstractLatentHeat{T} <: AbstractMaterialParam end

export compute_latent_heat,                  # calculation routines
    compute_latent_heat!,
    param_info,
    ConstantLatentHeat                      # constant (as source)
    
include("../Computations.jl")

# Constant  -------------------------------------------------------
"""
    ConstantLatentHeat(Q_L=400kJ/kg)
    
Set a constant latent heat:
```math  
    Q_L  = cst
```
where ``Q_L`` is the latent heat [``kJ/kg``].
"""
@with_kw_noshow struct ConstantLatentHeat{T,U} <: AbstractLatentHeat{T}
    Q_L::GeoUnit{T,U} = 400kJ / kg                # Latent heat
end
ConstantLatentHeat(args...) = ConstantLatentHeat(convert.(GeoUnit, args)...)

function param_info(s::ConstantLatentHeat) # info about the struct
    return MaterialParamsInfo(; Equation=L"Q_L = cst")
end

# Calculation routine
function (s::ConstantLatentHeat{_T})(; kwargs...) where {_T}
    @unpack_val Q_L = s

    return Q_L
end

compute_latent_heat(s::ConstantLatentHeat; kwargs...) = s()

function (s::ConstantLatentHeat)(I::Integer...)
    @unpack_val Q_L = s

    return fill(Q_L, I...)
end

# Print info 
function show(io::IO, g::ConstantLatentHeat)
    return print(io, "Constant latent heat: Q_L=$(Value(g.Q_L))")
end
#-------------------------------------------------------------------------


# Help info for the calculation routines
"""
    Ql = compute_latent_heat(s:<AbstractLatentHeat)

Returns the latent heat `Q_L`

"""

# Computational routines needed for computations with the MaterialParams structure 
function compute_latent_heat(s::AbstractMaterialParamsStruct, args)
    if isempty(s.LatentHeat)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return s.LatentHeat[1](args)
    end
end

# add methods programmatically
for myType in (:ConstantLatentHeat,)
    @eval begin
        (s::$(myType))(args) = s(; args...)
        compute_latent_heat(s::$(myType), args) = s(args)
    end
end

compute_latent_heat(args::Vararg{Any, N}) where N = compute_param(compute_latent_heat, args...)
compute_latent_heat!(args::Vararg{Any, N}) where N = compute_param!(compute_latent_heat, args...)

end
