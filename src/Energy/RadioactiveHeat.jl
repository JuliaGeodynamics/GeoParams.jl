module RadioactiveHeat

# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, @extractors, add_extractor_functions
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractRadioactiveHeat{T} <: AbstractMaterialParam end

export compute_radioactive_heat, # calculation routines
    compute_radioactive_heat!,
    param_info,
    ConstantRadioactiveHeat, # constant
    ExpDepthDependentRadioactiveHeat

include("../Computations.jl")
#include("../Utils.jl")
# Constant  -------------------------------------------------------
"""
    ConstantRadioactiveHeat(H_r=1e-6Watt/m^3)
    
Set a constant radioactive heat:
```math  
    H_r  = cst
```
where ``H_r`` is the radioactive heat source [``Watt/m^3``].
"""
@with_kw_noshow struct ConstantRadioactiveHeat{T, U} <: AbstractRadioactiveHeat{T}
    H_r::GeoUnit{T, U} = 1.0e-6Watt / m^3
end
ConstantRadioactiveHeat(args...) = ConstantRadioactiveHeat(convert.(GeoUnit, args)...)

function param_info(s::ConstantRadioactiveHeat) # info about the struct
    return MaterialParamsInfo(; Equation = L"H_r = cst")
end

# Calculation routine
function (s::ConstantRadioactiveHeat)(; kwargs...)
    @unpack_val H_r = s

    return H_r
end

compute_radioactive_heat(s::ConstantRadioactiveHeat; kwargs...) = s()

function (s::ConstantRadioactiveHeat)(I::Integer...)
    @unpack_val H_r = s

    return fill(H_r, I...)
end

# Print info
function show(io::IO, g::ConstantRadioactiveHeat)
    return print(io, "Constant radioactive heat: H_r=$(Value(g.H_r))")
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
@with_kw_noshow struct ExpDepthDependentRadioactiveHeat{T, U, U1} <:
    AbstractRadioactiveHeat{T}
    H_0::GeoUnit{T, U} = 1.0e-6Watt / m^3
    h_r::GeoUnit{T, U1} = 10.0e3m
    z_0::GeoUnit{T, U1} = 0m
end
function ExpDepthDependentRadioactiveHeat(args::Vararg{Any, N}) where {N}
    return ExpDepthDependentRadioactiveHeat(convert.(GeoUnit, args)...)
end

function param_info(s::ExpDepthDependentRadioactiveHeat) # info about the struct
    return MaterialParamsInfo(; Equation = L"H_r = H_0 \\exp(-(z-z_0)/h_r)")
end

# Calculation routines
function (s::ExpDepthDependentRadioactiveHeat{_T})(; z::_T = zero(_T), kwargs...) where {_T}
    @unpack_val H_0, z_0, h_r = s

    H_r = H_0 * exp(-(z - z_0) / h_r)

    return H_r
end

function compute_radioactive_heat(
        s::ExpDepthDependentRadioactiveHeat{_T}; z::_T = zero(_T)
    ) where {_T}
    return s(; z = z)
end

# Calculation routine
function (s::ExpDepthDependentRadioactiveHeat)(
        z::AbstractArray{_T, N}; kwargs...
    ) where {_T, N}
    Hr = similar(z)
    @unpack_val H_0, z_0, h_r = s

    @. Hr = H_0 * exp(-(z - z_0) / h_r)
    return Hr
end

function (s::ExpDepthDependentRadioactiveHeat)(
        z::AbstractArray{_T, N}, args...
    ) where {_T, N}
    return s(z; args...)
end
function compute_radioactive_heat(
        s::ExpDepthDependentRadioactiveHeat, z::AbstractArray{_T, N}, args...
    ) where {_T, N}
    return s(z; args...)
end

# Print info
function show(io::IO, g::ExpDepthDependentRadioactiveHeat)
    return print(
        io,
        "Exponential depth-dependent radioactive heat: H_r=$(Value(g.H_0)) exp(-(z-$(Value(g.z_0)))/$(Value(g.h_r)))",
    )
end
#-------------------------------------------------------------------------

# Computational routines needed for computations with the MaterialParams structure
function compute_radioactive_heat(s::AbstractMaterialParamsStruct, args)
    if isempty(s.RadioactiveHeat)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return s.RadioactiveHeat[1](args)
    end
end

# Help info for the calculation routines
"""
    H_r = compute_radioactive_heat(s:<AbstractRadioactiveHeat)

Returns the radioactive heat `H_r`

"""
#compute_radioactive_heat()

"""
    compute_radioactive_heat!(H_r, s:<AbstractRadioactiveHeat, z)

In-place computation of radioactive heat `H_r`

"""
#compute_radioactive_heat!()

# add methods programmatically
for myType in (:ExpDepthDependentRadioactiveHeat, :ConstantRadioactiveHeat)
    @eval begin
        (s::$(myType))(args) = s(; args...)
        compute_radioactive_heat(s::$(myType), args) = s(args)
    end
end

compute_radioactive_heat(args::Vararg{Any, N}) where {N} = compute_param(compute_radioactive_heat, args...)
compute_radioactive_heat!(args::Vararg{Any, N}) where {N} = compute_param!(compute_radioactive_heat, args...)

# extractor methods
for type in (ConstantRadioactiveHeat, ExpDepthDependentRadioactiveHeat)
    @extractors(type, :RadioactiveHeat)
end

end
