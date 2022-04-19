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
        ConstantRadioactiveHeat,                  # constant
        ExpDepthDependentRadioactiveHeat
        

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

# Exponential Depth-dependent  -------------------------------------------------------
"""
    ExpDepthDependent(H_0=1e-6Watt/m^3, h_r=10e3m, z_0=0m)
    
Sets an exponential depth-dependent radioactive 
```math  
    H_r  = H_0 \\exp \\left( {- {(z - z_0) \\over h_r}} \\right)
```
where ``H_0`` is the radioactive heat source [``Watt/m^3``] at ``z=z_0`` which decays with depth over a characteristic distance ``h_r``.
"""
@with_kw_noshow struct ExpDepthDependentRadioactiveHeat{T,U,U1} <: AbstractRadioactiveHeat{T}   
    H_0::GeoUnit{T,U}         =   1e-6Watt/m^3             
    h_r::GeoUnit{T,U1}        =   10e3m             
    z_0::GeoUnit{T,U1}        =   0m             
end
ExpDepthDependentRadioactiveHeat(a...) = ExpDepthDependentRadioactiveHeat(convert.(GeoUnit,a)...)

function param_info(s::ExpDepthDependentRadioactiveHeat) # info about the struct
    return MaterialParamsInfo(Equation = L"H_r = H_0 \\exp(-(z-z_0)/h_r)")
end

# Calculation routines
function compute_radioactive_heat(s::ExpDepthDependentRadioactiveHeat, z::Quantity)
    @unpack_units H_0, z_0, h_r   = s
    
    H_r = H_0*exp(-(z-z_0)/h_r);

    return H_r
end

function compute_radioactive_heat(s::ExpDepthDependentRadioactiveHeat, z::Number)
    @unpack_val H_0, z_0, h_r   = s
    
    H_r = H_0*exp(-(z-z_0)/h_r);

    return H_r
end

# Calculation routine
function compute_radioactive_heat!(Hr_array::AbstractArray{_T}, s::ExpDepthDependentRadioactiveHeat{_T}, z::_T=zero(_T)) where _T
    @unpack_val H_0, z_0, h_r   = s
        
    Hr_array .= H_0.*exp(-(z .- z_0)./h_r);
    return nothing
end

function compute_radioactive_heat!(Hr_array::AbstractArray{_T,N}, s::ExpDepthDependentRadioactiveHeat{_T}, z::AbstractArray{_T,N}) where {N,_T}
    @unpack_val H_0, z_0, h_r   = s
        
    Hr_array .= H_0.*exp(-(z .- z_0)./h_r);
    return nothing
end

# Print info 
function show(io::IO, g::ExpDepthDependentRadioactiveHeat)  
    print(io, "Exponential depth-dependent radioactive heat: H_r=$(g.H_0.val) exp(-(z-$(g.z_0.val))/$(g.h_r.val))")  
end
#-------------------------------------------------------------------------


# Help info for the calculation routines
"""
    H_r = compute_radioactive_heat(s:<AbstractRadioactiveHeat)

Returns the radioactive heat `H_r`

"""
compute_radioactive_heat()



end