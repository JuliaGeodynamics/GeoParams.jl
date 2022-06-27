module LatentHeat

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractLatentHeat{T} <: AbstractMaterialParam end

export compute_latent_heat,                  # calculation routines
    compute_latent_heat!,
    param_info,
    ConstantLatentHeat                  # constant

include("../Computations.jl")
include("../Utils.jl")

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

compute_latent_heat(s::ConstantLatentHeat{_T}; kwargs...) where {_T} = s()

function (s::ConstantLatentHeat{_T})(I::Integer...) where {_T}
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
#compute_latent_heat()

# Computational routines needed for computations with the MaterialParams structure 
function compute_latent_heat(s::AbstractMaterialParamsStruct, args)
    if isempty(s.LatentHeat)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return s.LatentHeat[1](args)
    end
end

# add methods programatically
for myType in (:ConstantLatentHeat,)
    @eval begin
        (s::$(myType))(args) = s(; args...)
        compute_latent_heat(s::$(myType), args) = s(args)
    end
end

compute_latent_heat(args...) = compute_param(compute_latent_heat, args...)
compute_latent_heat!(args...) = compute_param!(compute_latent_heat, args...)

end
