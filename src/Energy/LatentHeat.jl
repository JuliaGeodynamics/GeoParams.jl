module LatentHeat

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractLatentHeat{T} <: AbstractMaterialParam end

export  compute_latent_heat,                  # calculation routines
        param_info,
        ConstantLatentHeat                  # constant
        

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
    Q_L::GeoUnit{T,U}         =   400kJ/kg                # Latent heat
end
ConstantLatentHeat(a...) = ConstantLatentHeat(convert.(GeoUnit,a)...)

function param_info(s::ConstantLatentHeat) # info about the struct
    return MaterialParamsInfo(Equation = L"Q_L = cst")
end

# Calculation routine
function compute_latent_heat(s::ConstantLatentHeat)
    @unpack_val Q_L   = s
    
    return Q_L
end

# Print info 
function show(io::IO, g::ConstantLatentHeat)  
    print(io, "Constant latent heat: Q_L=$(g.Q_L.val)")  
end
#-------------------------------------------------------------------------


# Help info for the calculation routines
"""
    Ql = compute_latent_heat(s:<AbstractLatentHeat)

Returns the latent heat `Q_L`

"""
compute_latent_heat()



end