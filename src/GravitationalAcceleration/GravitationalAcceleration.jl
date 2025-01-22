module GravitationalAcceleration

# This implements the gravitational acceleration

using Parameters, LaTeXStrings, Unitful, StaticArrays
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
import Base.show, GeoParams.param_info
using ..MaterialParameters: MaterialParamsInfo

abstract type AbstractGravity{_T} <: AbstractMaterialParam end

export compute_gravity, # calculation routines
    ConstantGravity, # constant
    DippingGravity, # constant with dip and strike angles
    param_info

# Constant Gravity -------------------------------------------------------
"""
    GravityConstant(g=9.81m/s^2)
    
Set a constant value for the gravitational acceleration:
```math  
    g  = 9.81 m s^{-2}
```
"""
@with_kw_noshow struct ConstantGravity{_T, U} <: AbstractGravity{_T}
    g::GeoUnit{_T, U} = 9.81m / s^2               # gravitational acceleration
end
ConstantGravity(args...) = ConstantGravity(convert.(GeoUnit, args)...)

function param_info(s::ConstantGravity) # info about the struct
    return MaterialParamsInfo(; Equation = L"g = 9.81 m s^{-2}")
end

# Calculation routine
function compute_gravity(s::ConstantGravity)
    @unpack_val g = s

    return g
end

# Print info
function show(io::IO, d::ConstantGravity)
    return print(io, "Gravitational acceleration: g=$(UnitValue(d.g))")
end
#-------------------------------------------------------------------------

# Constant Gravity -------------------------------------------------------
"""
    DippingGravity(g=9.81m/s^2)
    
Set a constant value for the gravitational acceleration with dip and strike angles:
```math  
    g  = R_z  R_y  9.81 m s^{-2} 
```
"""
@kwdef struct DippingGravity{_T, U} <: AbstractGravity{_T}
    g::GeoUnit{_T, U} = 9.81m / s^2 # gravitational acceleration
    gx::GeoUnit{_T, U} = 0.0e0m / s^2  # gravitational acceleration
    gy::GeoUnit{_T, U} = 0.0e0m / s^2  # gravitational acceleration
    gz::GeoUnit{_T, U} = 9.81m / s^2 # gravitational acceleration
end
DippingGravity(args...) = DippingGravity(convert.(GeoUnit, args)...)

function DippingGravity(α::T1, θ::T2, g::T3) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    sinα, cosα = sincosd(90 - α)
    sinθ, cosθ = sincosd(θ)
    gᵢ = @SVector [zero(T), zero(T), T(g)]

    Ry = @SMatrix [
        cosα    zero(T) sinα
        zero(T)  one(T) zero(T)
        -sinα    zero(T) cosα
    ]

    Rz = @SMatrix [
        cosθ    -sinθ   zero(T)
        sinθ    cosθ    zero(T)
        zero(T) zero(T) one(T)
    ]
    g′ = Rz * (Ry * gᵢ)

    return DippingGravity(; g = g, gx = g′[1], gy = g′[2], gz = g′[3])
end

function param_info(s::DippingGravity) # info about the struct
    return MaterialParamsInfo(; Equation = L"g = ($(s.gx, s.gy, s.gz)) m s^{-2}")
end

# Calculation routine
function compute_gravity(s::DippingGravity)
    @unpack_val gx, gy, gz = s
    return gx, gy, gz
end

# Print info
function show(io::IO, d::DippingGravity)
    return print(io, "Gravitational acceleration: g=$((UnitValue(d.gz), UnitValue(d.gy), UnitValue(d.gz)))")
end
#-------------------------------------------------------------------------


# Calculation routine
@generated function compute_gravity(
        MatParam::NTuple{N, AbstractMaterialParamsStruct}, Phase::Integer
    ) where {N}
    return quote
        @inline
        Base.Cartesian.@nexprs $N i ->
        (MatParam[i].Phase == Phase) && return compute_gravity(MatParam[i].Gravity[1])
    end
end

compute_gravity(MatParam::AbstractMaterialParamsStruct) = compute_gravity(MatParam.Gravity[1])


# Help info for the calculation routines
"""
compute_gravity(s:<AbstractGravity)

Returns the gravitational acceleration 

"""
compute_gravity

end
