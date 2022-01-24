module GravitationalAcceleration

# This implements the gravitational acceleration

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import Base.show

abstract type AbstractGravity{_T} <: AbstractMaterialParam end

export  ComputeGravity,        # calculation routines
        ConstantGravity        # constant


# Constant Gravity -------------------------------------------------------
"""
    GravityConstant(g=9.81m/s^2)
    
Set a constant value for the gravitational acceleration:
```math  
    g  = 9.81 m s^{-2}
```
"""
@with_kw_noshow struct ConstantGravity{_T,U}   <: AbstractGravity{_T}
 #   equation::LaTeXString   =   L"g = 9.81 m s^{-2}"     
    g::GeoUnit{_T,U}              =   9.81m/s^2               # gravitational acceleration
end
ConstantGravity(args...) = ConstantGravity(convert.(GeoUnit,args)...) 

# Calculation routine
function ComputeGravity(s::ConstantGravity{_T}) where _T
    @unpack g   = s
    
    return g*1.0   # multiply with 1.0, to return Float64
end

# Print info 
function show(io::IO, d::ConstantGravity{_T})  where _T
    print(io, "Gravitational acceleration: g=$(d.g.val)")  
end
#-------------------------------------------------------------------------




# Help info for the calculation routines
"""
ComputeGravity(s:<AbstractGravity)

Returns the gravitational acceleration 

"""
ComputeGravity



end