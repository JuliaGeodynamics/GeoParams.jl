module RadioactiveHeat

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractRadioactiveHeat{T} <: AbstractMaterialParam end

export  compute_radioactive_heat,                  # calculation routines
        param_info,
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
@with_kw_noshow struct ConstantRadioactiveHeat{T,U} <: AbstractRadioactiveHeat{T}   
    H_r::GeoUnit{T,U}         =   1e-6Watt/m^3             
end
ConstantRadioactiveHeat(a...) = ConstantRadioactiveHeat(convert.(GeoUnit,a)...)

function param_info(s::ConstantRadioactiveHeat) # info about the struct
    return MaterialParamsInfo(Equation = L"H_r = cst")
end

# Calculation routine
function compute_radioactive_heat(s::ConstantRadioactiveHeat)
    @unpack_val H_r   = s
   
    return H_r
end

# Print info 
function show(io::IO, g::ConstantRadioactiveHeat)  
    print(io, "Constant radioactive heat: H_r=$(g.H_r.val)")  
end
#-------------------------------------------------------------------------


# Help info for the calculation routines
"""
    H_r = compute_radioactive_heat(s:<AbstractRadioactiveHeat)

Returns the radioactive heat `H_r`

"""
compute_radioactive_heat()



end