module GravitationalAcceleration

# This implements the gravitational acceleration

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import Base.show

abstract type AbstractGravity <: AbstractMaterialParam end

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
@with_kw_noshow mutable struct ConstantGravity <: AbstractGravity
    equation::LaTeXString   =   L"g = 9.81 m s^{-2}"     
    g::GeoUnit              =   9.81m/s^2               # gravitational acceleration
end

# Calculation routine
function ComputeGravity(s::ConstantGravity)
    @unpack g   = s
    
    return g*1.0   # multiply with 1.0, to return Float64
end

# Print info 
function show(io::IO, d::ConstantGravity)  
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