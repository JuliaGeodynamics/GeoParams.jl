module LatentHeat

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import Base.show

abstract type AbstractLatentHeat <: AbstractMaterialParam end

export  ComputeLatentHeat,                  # calculation routines
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
@with_kw_noshow mutable struct ConstantLatentHeat <: AbstractLatentHeat
    equation::LaTeXString   =   L"Q_L = cst"     
    Q_L::GeoUnit             =   400kJ/kg                # Latent heat
end

# Calculation routine
function ComputeLatentHeat(s::ConstantLatentHeat)
    @unpack Q_L   = s
    
    return Q_L*1.0
end

# Print info 
function show(io::IO, g::ConstantLatentHeat)  
    print(io, "Constant latent heat: Q_L=$(g.Q_L.val)")  
end
#-------------------------------------------------------------------------


# Help info for the calculation routines
"""
    Ql = ComputeLatentHeat(s:<AbstractLatentHeat)

Returns the latent heat `Q_L`

"""
ComputeLatentHeat()



end