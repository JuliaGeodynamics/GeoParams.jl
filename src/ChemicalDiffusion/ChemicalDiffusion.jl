module ChemicalDiffusion

using GeoParams
using Unitful

using Parameters, LaTeXStrings, Unitful, MuladdMacro
using ..Units
using GeoParams: AbstractMaterialParam, AbstractConstitutiveLaw, AbstractComposite
import GeoParams: param_info, ptr2string
using ..MaterialParameters: MaterialParamsInfo
import Base.show

abstract type AbstractChemicalDiffusion{T} <: AbstractMaterialParam end

@inline diffusion_database(f::F) where {F} = first(f())
@inline diffusion_database_info(f::F) where {F} = last(f())

@inline precision(::AbstractChemicalDiffusion{T}) where T = T

"""
    DiffusionData(; Name, Mineral, Formula, Species, D0, D0_2σ, Ea, Ea_2σ, ΔV, ΔV_2σ, T_range, P0, Orientation, Crystallography, Buffer)

Defines the diffusion data for the chemical diffusion of a given phase and species from an experiment.

The diffusion coefficient `D` [\\mathrm{[m^2/s]}] is given by an Arrhenius equation:
```math
    D = D0 * \\exp\\left(-\\frac{Ea + PΔV} {RT}\\right)
```
where
- ``D0`` is the pre-exponential factor [\\mathrm{[m^2/s]}],
- ``Ea`` is the activation energy [\\mathrm{[J/mol]}],
- ``ΔV`` is the activation volume [\\mathrm{[cm^3/mol]}],
- ``P`` is the pressure [\\mathrm{[Pa]},
- ``T`` is the temperature [\\mathrm{[K]}],
- ``R`` is the gas constant [\\mathrm{[J/(mol K)]}],
"""
struct DiffusionData{T, U1, U2, U3, U4, U5, U6, U7, U8, U9} <: AbstractChemicalDiffusion{T}
    Name::Ptr{UInt8}  # name of the diffusion experiment and paper
    Mineral::Ptr{UInt8}  # name of the mineral
    Formula::Ptr{UInt8}  # chemical formula of the mineral
    Species::Ptr{UInt8}  # element or species being diffused
    D0::GeoUnit{T, U1} # pre-exponential factor
    D0_2σ::GeoUnit{T, U2} # uncertainty at 2σ of the pre-exponential factor
    Ea::GeoUnit{T, U3} # activation energy
    Ea_2σ::GeoUnit{T, U4} # uncertainty at 2σ of the activation energy
    ΔV::GeoUnit{T, U5} # activation volume
    ΔV_2σ::GeoUnit{T, U6} # uncertainty at 2σ of the activation volume
    R::GeoUnit{T, U7}  # gas constant
    T_range_min::GeoUnit{T, U8}  # minimum temperature of the T_range
    T_range_max::GeoUnit{T, U8}  # maximum temperature of the T_range
    P0::GeoUnit{T, U9} # pressure of calibration
    Orientation::Ptr{UInt8}  # Crystal orientation from the diffusion experiment
    Crystallography::Ptr{UInt8}  # Crystallographic system of the mineral
    Buffer::Ptr{UInt8}  # Buffer condition (e.g., NNO) during the experiment
    Fluid::Ptr{UInt8}  # Fluid condition (e.g., anhydrous) during the experiment

    function DiffusionData(;
        Name            = "",  # name of the diffusion experiment and paper
        Mineral         = "",  # name of the mineral
        Formula         = "",  # chemical formula of the mineral
        Species         = "",  # element or species being diffused
        D0              = 0.0m^2/s,  # pre-exponential factor
        D0_2σ           = 0.0m^2/s,  # uncertainty at 2σ of the pre-exponential factor
        Ea              = 0.0J/mol,  # activation energy
        Ea_2σ           = 0.0J/mol,  # uncertainty at 2σ of the activation energy
        ΔV              = 0.0cm^3/mol,  # activation volume
        ΔV_2σ           = 0.0cm^3/mol,  # uncertainty at 2σ of the activation volume
        R               = Unitful.R,  # gas constant
        T_range_min     = 0.0K,  # minimum temperature of the T_range
        T_range_max     = 0.0K,  # maximum temperature of the T_range
        P0              = 0.0Pa,  # pressure of calibration
        Orientation     = "",  # Crystal orientation from the diffusion experiment
        Crystallography = "",  # Crystallographic system of the mineral
        Buffer          = "",  # Buffer condition (e.g., NNO) during the experiment
        Fluid           = ""  # Fluid condition (e.g., anhydrous) during the experiment
    )

        # Convert to GeoUnits
        D0U          = convert(GeoUnit, D0)
        D0_2σU       = convert(GeoUnit, D0_2σ)
        EaU          = convert(GeoUnit, Ea)
        Ea_2σU       = convert(GeoUnit, Ea_2σ)
        ΔVU          = convert(GeoUnit, ΔV)
        ΔV_2σU       = convert(GeoUnit, ΔV_2σ)
        RU           = convert(GeoUnit, R)
        T_range_minU = convert(GeoUnit, T_range_min)
        T_range_maxU = convert(GeoUnit, T_range_max)
        P0U          = convert(GeoUnit, P0)
        # Extract struct types
        T               = typeof(D0U).types[1]
        U1              = typeof(D0U).types[2]
        U2              = typeof(D0_2σU).types[2]
        U3              = typeof(EaU).types[2]
        U4              = typeof(Ea_2σU).types[2]
        U5              = typeof(ΔVU).types[2]
        U6              = typeof(ΔV_2σU).types[2]
        U7              = typeof(RU).types[2]
        U8              = typeof(T_range_minU).types[2]
        U9              = typeof(P0U).types[2]
        name            = pointer(ptr2string(Name))
        mineral         = pointer(ptr2string(Mineral))
        formula         = pointer(ptr2string(Formula))
        species         = pointer(ptr2string(Species))
        orientation     = pointer(ptr2string(Orientation))
        crystallography = pointer(ptr2string(Crystallography))
        buffer          = pointer(ptr2string(Buffer))
        fluid           = pointer(ptr2string(Fluid))

        # Create struct
        return new{T,U1,U2,U3,U4,U5,U6,U7,U8,U9}(
            name, mineral, formula, species, D0U, D0_2σU, EaU, Ea_2σU, ΔVU, ΔV_2σU, RU, T_range_minU, T_range_maxU, P0U, orientation, crystallography, buffer, fluid
        )
    end
end

function param_info(data::DiffusionData) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"D = D0 * \exp\left(-\frac{Ea + PΔV} {RT}\right)",
    )
end


# Calculation routines for the diffusion coefficient without units
@inline function compute_D(data::DiffusionData; T=1K, P=0Pa, kwargs...)

    if P isa Quantity && T isa Quantity
        @unpack_units D0, Ea, ΔV, P0, R = data
    else
        @unpack_val D0, Ea, ΔV, P0, R = data
    end

    D = D0 * exp(- (Ea + (P - P0) * ΔV) / (R * T))

    return D
end


function compute_D!(
    D,
    data::DiffusionData;
    T=ones(size(D))K,
    P=zeros(size(D))Pa,
    kwargs...) 

    @inbounds for i in eachindex(D)
        D[i] = compute_D(data; T=T[i], P=P[i], kwargs...)
    end
end


function show(io::IO, g::DiffusionData)
    return print(
        io,
        "DiffusionData: Mineral = $(unsafe_string(g.Mineral)), Species = $(unsafe_string(g.Species)), D0 = $(Value(g.D0)), Ea = $(Value(g.Ea)), ΔV = $(Value(g.ΔV))",
    )
end


export AbstractChemicalDiffusion,
       DiffusionData,
       compute_D,
       compute_D!


"""
    SetChemicalDiffusion["Name of Chemical Diffusion"]
"""
function SetChemicalDiffusion(
    name:: F;
    D0          = nothing,
    D0_2σ       = nothing,
    Ea          = nothing,
    Ea_2σ       = nothing,
    ΔV          = nothing,
    ΔV_2σ       = nothing,
    T_range_min = nothing,
    T_range_max = nothing,
    P0          = nothing,
    ) where F
    kwargs = (; D0, D0_2σ, Ea, Ea_2σ, ΔV, ΔV_2σ, T_range_min, T_range_max, P0)
    Transform_ChemicalDiffusion(name, kwargs)
end

function SetChemicalDiffusion(
    name::F,
    CharDim::GeoUnits{T};
    D0          = nothing,
    D0_2σ       = nothing,
    Ea          = nothing,
    Ea_2σ       = nothing,
    ΔV          = nothing,
    ΔV_2σ       = nothing,
    T_range_min = nothing,
    T_range_max = nothing,
    P0          = nothing,
    ) where {F, T<:Union{GEO, SI}}
    kwargs = (; D0, D0_2σ, Ea, Ea_2σ, ΔV, ΔV_2σ, T_range_min, T_range_max, P0)
    nondimensionalize(Transform_ChemicalDiffusion(name, kwargs), CharDim)
end

"""
    Transform_ChemicalDiffusion(name)
Transforms units from MPa, kJ etc. to basic units such as Pa, J etc.
"""
Transform_ChemicalDiffusion(name::F) where F = Transform_ChemicalDiffusion(diffusion_database(name))
Transform_ChemicalDiffusion(name::F, kwargs::NamedTuple) where F = Transform_ChemicalDiffusion(diffusion_database(name), kwargs)

function Transform_ChemicalDiffusion(name::F, CharDim::GeoUnits{U}) where {F, U<:Union{GEO,SI}}
    Transform_ChemicalDiffusion(diffusion_database(name), CharDim)
end

function Transform_ChemicalDiffusion(p::AbstractChemicalDiffusion{T}, CharDim::GeoUnits{U}) where {T,U<:Union{GEO,SI}}
    nondimensionalize(Transform_ChemicalDiffusion(p), CharDim)
end

function Transform_ChemicalDiffusion(pp::AbstractChemicalDiffusion{T}) where T

    D0          = Value(pp.D0)
    D0_2σ       = Value(pp.D0_2σ)
    Ea          = Value(pp.Ea)
    Ea_2σ       = Value(pp.Ea_2σ)
    ΔV          = Value(pp.ΔV)
    ΔV_2σ       = Value(pp.ΔV_2σ)
    T_range_min = Value(pp.T_range_min)
    T_range_max = Value(pp.T_range_max)
    P0          = Value(pp.P0)

    D0_SI          = uconvert(m^2 / s, D0 )
    D0_2σ_SI       = uconvert(m^2 / s, D0_2σ)
    Ea_SI          = uconvert(J / mol, pp.Ea)
    Ea_2σ_SI       = uconvert(J / mol, pp.Ea_2σ)
    ΔV_SI          = uconvert(m^3 / mol, pp.ΔV)
    ΔV_2σ_SI       = uconvert(m^3 / mol, pp.ΔV_2σ)
    T_range_min_SI = uconvert(K, T_range_min)
    T_range_max_SI = uconvert(K, T_range_max)
    P0_SI          = uconvert(Pa, P0)


    args       = (Name=unsafe_string(pp.Name),
                  Mineral=unsafe_string(pp.Mineral),
                  Formula=unsafe_string(pp.Formula),
                  Species=unsafe_string(pp.Species),
                  D0=D0_SI, D0_2σ=D0_2σ_SI,
                  Ea=Ea_SI, Ea_2σ=Ea_2σ_SI,
                  ΔV=ΔV_SI, ΔV_2σ=ΔV_2σ_SI,
                  T_range_min_SI=T_range_min_SI,
                  T_range_max_SI=T_range_max_SI,
                  P0=P0_SI,
                  Orientation=unsafe_string(pp.Orientation),
                  Crystallography=unsafe_string(pp.Crystallography),
                  Buffer=unsafe_string(pp.Buffer),
                  Fluid=unsafe_string(pp.Fluid)
                  )

    return DiffusionData(; args...)
end

function Transform_ChemicalDiffusion(pp::AbstractChemicalDiffusion{T}, kwargs::NamedTuple) where T

    (; D0, D0_2σ, Ea, Ea_2σ, ΔV, ΔV_2σ, T_range_min, T_range_max, P0) = kwargs

    D0_new          = isnothing(D0) ? Value(pp.D0) : Value(GeoUnit(D0))
    D0_2σ_new       = isnothing(D0_2σ) ? Value(pp.D0_2σ) : Value(GeoUnit(D0_2σ))
    Ea_new          = isnothing(Ea) ? Value(pp.Ea) : Value(GeoUnit(Ea))
    Ea_2σ_new       = isnothing(Ea_2σ) ? Value(pp.Ea_2σ) : Value(GeoUnit(Ea_2σ))
    ΔV_new          = isnothing(ΔV) ? Value(pp.ΔV) : Value(GeoUnit(ΔV))
    ΔV_2σ_new       = isnothing(ΔV_2σ) ? Value(pp.ΔV_2σ) : Value(GeoUnit(ΔV_2σ))
    T_range_min_new = isnothing(T_range_min) ? Value(pp.T_range_min) : Value(GeoUnit(T_range_min))
    T_range_max_new = isnothing(T_range_max) ? Value(pp.T_range_max) : Value(GeoUnit(T_range_max))
    P0_new          = isnothing(P0) ? Value(pp.P0) : Value(GeoUnit(P0))

    D0_SI          = uconvert(m^2 / s, D0_new)
    D0_2σ_SI       = uconvert(m^2 / s, D0_2σ_new)
    Ea_SI          = uconvert(J / mol, Ea_new)
    Ea_2σ_SI       = uconvert(J / mol, Ea_2σ_new)
    ΔV_SI          = uconvert(m^3 / mol, ΔV_new)
    ΔV_2σ_SI       = uconvert(m^3 / mol, ΔV_2σ_new)
    T_range_min_SI = uconvert(K, T_range_min_new)
    T_range_max_SI = uconvert(K, T_range_max_new)
    P0_SI          = uconvert(Pa, P0_new)

    args       = (Name=unsafe_string(pp.Name),
                  Mineral=unsafe_string(pp.Mineral),
                  Formula=unsafe_string(pp.Formula),
                  Species=unsafe_string(pp.Species),
                  D0=D0_SI, D0_2σ=D0_2σ_SI,
                  Ea=Ea_SI, Ea_2σ=Ea_2σ_SI,
                  ΔV=ΔV_SI, ΔV_2σ=ΔV_2σ_SI,
                  T_range_min=T_range_min_SI,
                  T_range_max=T_range_max_SI,
                  P0=P0_SI,
                  Orientation=unsafe_string(pp.Orientation),
                  Crystallography=unsafe_string(pp.Crystallography),
                  Buffer=unsafe_string(pp.Buffer),
                  Fluid=unsafe_string(pp.Fluid)
                  )

    return DiffusionData(; args...)
end

export SetChemicalDiffusion,
       Transform_ChemicalDiffusion


# load collection of chemical diffusion data
include("Data/Rutile/Rutile.jl")
using .Rutile

export Rutile


end
