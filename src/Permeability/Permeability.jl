module Permeability

using Parameters, Unitful, LaTeXStrings, MuladdMacro
using ..Units
using ..PhaseDiagrams
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, @extractors, add_extractor_functions
import ..Units: isdimensional
using ..MaterialParameters: No_MaterialParam, MaterialParamsInfo
import Base.show, GeoParams.param_info

include("../Computations.jl")

abstract type AbstractPermeability{T} <: AbstractMaterialParam end

export compute_permeability,     # calculation routines
    compute_permeability!,       # in place calculation
    param_info,                  # info about the parameters
    AbstractPermeability,
    ConstantPermeability,        # constant
    HazenPermeability,           # Hazen equation
    PowerLawPermeability,        # Power-law permeability
    CarmanKozenyPermeability     # Carman-Kozeny permeability

# Define "empty" computational routines in case nothing is defined
function compute_permeability!(
    k::_T, s::No_MaterialParam{_T}; ϕ::_T=zero(_T)
) where {_T}
    return zero(_T)
end
function compute_permeability(s::No_MaterialParam{_T}; ϕ::_T=zero(_T)) where {_T}
    return zero(_T)
end

# Constant Permeability
@with_kw_noshow struct ConstantPermeability{_T,U} <: AbstractPermeability{_T}
    k::GeoUnit{_T,U} = 1e-12m^2 # permeability
end
ConstantPermeability(args...) = ConstantPermeability(convert.(GeoUnit, args)...)
isdimensional(s::ConstantPermeability) = isdimensional(s.k)

@inline compute_permeability(s::ConstantPermeability{_T}, args) where {_T} = s(; args...)
@inline compute_permeability(s::ConstantPermeability{_T})       where {_T} = s()

function param_info(s::ConstantPermeability)
    return MaterialParamsInfo(; Equation=L"k = cst")
end

function compute_permeability!(k::AbstractArray, s::ConstantPermeability; kwargs...)
    @unpack_val k_val = s
    k[:] .= k_val
    return nothing
end

function compute_permeability!(k::AbstractArray, s::ConstantPermeability, args)
    return compute_permeability!(k, s; args...)
end

function show(io::IO, g::ConstantPermeability)
    return print(io, "Constant permeability: k=$(UnitValue(g.k))")
end

# Hazen Permeability
@with_kw_noshow struct HazenPermeability{_T,U1,U2} <: AbstractPermeability{_T}
    C::GeoUnit{_T,U1} = 1.0 * NoUnits          # Hazen constant
    D10::GeoUnit{_T,U2} = 1e-4 * m   # Effective grain size
end
HazenPermeability(args...) = HazenPermeability(convert.(GeoUnit, args)...)
isdimensional(s::HazenPermeability) = isdimensional(s.D10)

function param_info(s::HazenPermeability)
    return MaterialParamsInfo(; Equation = L"k = C \cdot D_{10}^2")
end

function (s::HazenPermeability{_T})(; kwargs...) where {_T}
    if kwargs isa Quantity
        @unpack_units C, D10 = s
    else
        @unpack_val   C, D10 = s
    end

    return C * D10^2
end

@inline (s::HazenPermeability)(args)                = s(; args...)
@inline compute_permeability(s::HazenPermeability, args) = s(; args...)


function show(io::IO, g::HazenPermeability)
    return print(io, "Hazen permeability: k = C * D10^2; C=$(g.C); D10=$(g.D10)")
end

# Power-law Permeability
@with_kw_noshow struct PowerLawPermeability{_T,U1,U2,U3, U4} <: AbstractPermeability{_T}
    c::GeoUnit{_T,U1}  = 1.0 * NoUnits      # Power-law constant
    k0::GeoUnit{_T,U2} = 1e-12 * m^2           # reference permeability
    ϕ::GeoUnit{_T,U3}  = 1e-2 * NoUnits     # reference porosity
    n::GeoUnit{_T,U4}  = 3.0 * NoUnits      # exponent
end
PowerLawPermeability(args...) = PowerLawPermeability(convert.(GeoUnit, args)...)
isdimensional(s::PowerLawPermeability) = isdimensional(s.k0)

function param_info(s::PowerLawPermeability)
    return MaterialParamsInfo(; Equation = L"k = c \cdot k_0 \cdot \phi^n")
end

function (s::PowerLawPermeability{_T})(; ϕ=0e0, kwargs...) where {_T}
    if ϕ isa Quantity
        @unpack_units c, k0, n = s
    else
        @unpack_val   c, k0, n = s
    end

    return c * k0 * ϕ^n
end

@inline (s::PowerLawPermeability)(args)                = s(; args...)
@inline compute_permeability(s::PowerLawPermeability, args) = s(args)

function show(io::IO, g::PowerLawPermeability)
    return print(io, "Power-law permeability: k = c* k0 * ϕ^n; c = $(g.c), k0=$(g.k0); n=$(g.n)")
end

# Carman-Kozeny Permeability
@with_kw_noshow struct CarmanKozenyPermeability{_T,U1,U2,U3} <: AbstractPermeability{_T}
    c::GeoUnit{_T,U1} = 1.0 * NoUnits       # Carman-Kozeny constant
    ϕ0::GeoUnit{_T,U2} = 0.01 * NoUnits      # reference porosity
    n::GeoUnit{_T,U3}  = 3.0 * NoUnits     # exponent
end
CarmanKozenyPermeability(args...) = CarmanKozenyPermeability(convert.(GeoUnit, args)...)
# isdimensional(s::CarmanKozenyPermeability) = isdimensional(s.c)

function param_info(s::CarmanKozenyPermeability)
    return MaterialParamsInfo(; Equation = L"k = c \left(\frac{\phi}{\phi_0}\right)^n")
end

function (s::CarmanKozenyPermeability{_T})(; ϕ=1e-2, kwargs...) where {_T}
    if ϕ isa Quantity
        @unpack_units c, ϕ0, n = s
    else
        @unpack_val   c, ϕ0, n = s
    end

    return c * (ϕ / ϕ0)^n
end

@inline (s::CarmanKozenyPermeability)(args)                = s(; args...)
@inline compute_permeability(s::CarmanKozenyPermeability, args) = s(args)

function show(io::IO, g::CarmanKozenyPermeability)
    return print(io, "Carman-Kozeny permeability: k = c * (ϕ / ϕ0)^n; c=$(g.c); ϕ0=$(g.ϕ0); n=$(g.n)")
end

# extractor methods
for type in (ConstantPermeability, HazenPermeability, PowerLawPermeability, CarmanKozenyPermeability)
    @extractors(type, :Permeability)
end

end
