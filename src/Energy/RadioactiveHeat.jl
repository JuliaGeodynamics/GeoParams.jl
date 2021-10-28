module RadioactiveHeat

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import Base.show

abstract type AbstractRadioactiveHeat <: AbstractMaterialParam end

export  ComputeRadioactiveHeat,                  # calculation routines
        ConstantRadioactiveHeat                  # constant
        

# Constant  -------------------------------------------------------
"""
    ConstantRadioactiveHeat(H_r=1e-6Watt/m^3)
    
Set a constant radioactive heat:
```math  
    H_r  = cst
```
where ``H_r`` is the radioactive heat source [``Watt/m^3``].
"""
@with_kw_noshow mutable struct ConstantRadioactiveHeat <: AbstractRadioactiveHeat
    equation::LaTeXString   =   L"H_r = cst"     
    H_r::GeoUnit            =   1e-6Watt/m^3             
end

# Calculation routine
function ComputeRadioactiveHeat(s::ConstantRadioactiveHeat)
    @unpack H_r   = s
   
    return Value(H_r)
end

# Print info 
function show(io::IO, g::ConstantRadioactiveHeat)  
    print(io, "Constant radioactive heat: H_r=$(g.H_r.val)")  
end
#-------------------------------------------------------------------------


# Help info for the calculation routines
"""
    H_r = ComputeRadioactiveHeat(s:<AbstractRadioactiveHeat)

Returns the radioactive heat `H_r`

"""
ComputeRadioactiveHeat()



end