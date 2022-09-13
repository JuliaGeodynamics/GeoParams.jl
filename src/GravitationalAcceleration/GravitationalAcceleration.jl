module GravitationalAcceleration

# This implements the gravitational acceleration

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
import Base.show, GeoParams.param_info
using ..MaterialParameters: MaterialParamsInfo

abstract type AbstractGravity{_T} <: AbstractMaterialParam end

export compute_gravity,        # calculation routines
    ConstantGravity,        # constant
    param_info

# Constant Gravity -------------------------------------------------------
"""
    GravityConstant(g=9.81m/s^2)
    
Set a constant value for the gravitational acceleration:
```math  
    g  = 9.81 m s^{-2}
```
"""
@with_kw_noshow struct ConstantGravity{_T,U} <: AbstractGravity{_T}
    g::GeoUnit{_T,U} = 9.81m / s^2               # gravitational acceleration
end
ConstantGravity(args...) = ConstantGravity(convert.(GeoUnit, args)...)

function param_info(s::ConstantGravity) # info about the struct
    return MaterialParamsInfo(; Equation=L"g = 9.81 m s^{-2}")
end

# Calculation routine
function compute_gravity(s::ConstantGravity{_T}) where {_T}
    @unpack_val g = s

    return g
end

# Calculation routine
@inline @generated function compute_gravity(
    MatParam::NTuple{N,AbstractMaterialParamsStruct}, Phase::Integer
) where {N}
    quote
        Base.Cartesian.@nexprs $N i ->
            (MatParam[i].Phase == Phase) && return compute_gravity(MatParam[i].Gravity[1])
    end
end

# Print info 
function show(io::IO, d::ConstantGravity{_T}) where {_T}
    return print(io, "Gravitational acceleration: g=$(UnitValue(d.g))")
end
#-------------------------------------------------------------------------

# Help info for the calculation routines
"""
compute_gravity(s:<AbstractGravity)

Returns the gravitational acceleration 

"""
compute_gravity

end
