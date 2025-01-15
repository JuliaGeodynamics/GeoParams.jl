module ChemicalDiffusion

using GeoParams
using Unitful

using Base: Float64
using Parameters, LaTeXStrings, Unitful, MuladdMacro
using ..Units
using GeoParams: AbstractMaterialParam, AbstractConstitutiveLaw, AbstractComposite
import GeoParams: param_info, fastpow, pow_check, nphase, ntuple_idx, @print, @pow, ptr2string
using ..MaterialParameters: MaterialParamsInfo
import Base.show

abstract type AbstractChemicalDiffusion{T} <: AbstractMaterialParam end


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
struct DiffusionData{T, U1, U2, U3, U4, U5, U6} <: AbstractChemicalDiffusion{T}
    Name::Ptr{UInt8}  # name of the diffusion experiment and paper
    Mineral::Ptr{UInt8}  # name of the mineral
    Formula::Ptr{UInt8}  # chemical formula of the mineral
    Species::Ptr{UInt8}  # element or species being diffused
    D0::GeoUnit{T, U1} # pre-exponential factor
    D0_2σ::GeoUnit{T, U1} # uncertainty at 2σ of the pre-exponential factor
    Ea::GeoUnit{T, U2} # activation energy
    Ea_2σ::GeoUnit{T, U2} # uncertainty at 2σ of the activation energy
    ΔV::GeoUnit{T, U3} # activation volume
    ΔV_2σ::GeoUnit{T, U3} # uncertainty at 2σ of the activation volume
    R::GeoUnit{T, U4}  # gas constant
    T_range::Tuple{GeoUnit{T, U5}, GeoUnit{T, U5}} # temperature range
    P0::GeoUnit{T, U6} # pressure of calibration
    Orientation::Ptr{UInt8}  # Crystal orientation from the diffusion experiment
    Crystallography::Ptr{UInt8}  # Crystallographic system of the mineral
    Buffer::Ptr{UInt8}  # Buffer condition (e.g., NNO) during the experiment
    fluid::Ptr{UInt8}  # Fluid condition (e.g., anhydrous) during the experiment

    function DiffusionData(;
        Name="",  # name of the diffusion experiment and paper
        Mineral="",  # name of the mineral
        Formula="",  # chemical formula of the mineral
        Species="",  # element or species being diffused
        D0=0.0m^2/s,  # pre-exponential factor
        D0_2σ=0.0m^2/s,  # uncertainty at 2σ of the pre-exponential factor
        Ea=0.0J/mol,  # activation energy
        Ea_2σ=0.0J/mol,  # uncertainty at 2σ of the activation energy
        ΔV=0.0cm^3/mol,  # activation volume
        ΔV_2σ=0.0cm^3/mol,  # uncertainty at 2σ of the activation volume
        R=Unitful.R,  # gas constant
        T_range=(0.0K, 0.0K),  # temperature range
        P0=0.0Pa,  # pressure of calibration
        Orientation="",  # Crystal orientation from the diffusion experiment
        Crystallography="",  # Crystallographic system of the mineral
        Buffer="",  # Buffer condition (e.g., NNO) during the experiment
        Fluid=""  # Fluid condition (e.g., anhydrous) during the experiment
    )

        # Convert to GeoUnits
        D0U = convert(GeoUnit, D0)
        D0_2σU = convert(GeoUnit, D0_2σ)
        EaU = convert(GeoUnit, Ea)
        Ea_2σU = convert(GeoUnit, Ea_2σ)
        ΔVU = convert(GeoUnit, ΔV)
        ΔV_2σU = convert(GeoUnit, ΔV_2σ)
        RU = convert(GeoUnit, R)
        T_rangeU = convert.(GeoUnit, T_range)
        P0U = convert(GeoUnit, P0)
        # Extract struct types
        T = typeof(D0U).types[1]
        U1 = typeof(D0U).types[2]
        U2 = typeof(EaU).types[2]
        U3 = typeof(ΔVU).types[2]
        U4 = typeof(RU).types[2]
        U5 = typeof(T_rangeU[1]).types[2]
        U6 = typeof(P0U).types[2]
        name = pointer(ptr2string(Name))
        mineral = pointer(ptr2string(Mineral))
        formula = pointer(ptr2string(Formula))
        species = pointer(ptr2string(Species))
        orientation = pointer(ptr2string(Orientation))
        crystallography = pointer(ptr2string(Crystallography))
        buffer = pointer(ptr2string(Buffer))
        fluid = pointer(ptr2string(Fluid))

        # Create struct
        return new{T,U1,U2,U3,U4,U5,U6}(
            name, mineral, formula, species, D0U, D0_2σU, EaU, Ea_2σU, ΔVU, ΔV_2σU, RU, T_rangeU, P0U, orientation, crystallography, buffer, fluid
        )
    end
end

function param_info(data::DiffusionData) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"D = D0 * \exp\left(-\frac{Ea + PΔV} {RT}\right)",
    )
end

# Calculation routines for the diffusion coefficient
@inline function compute_D(data::DiffusionData; T=1K, P=0Pa, kwargs...)
    @unpack_units D0, Ea, ΔV, P0, R = data

    D = @pow D0 * exp(- (Ea + (P - P0) * ΔV) / (R * T))

    return D
end

function compute_D!(D::AbstractVector, data::DiffusionData; T=ones(size(D)), P=zeros(size(data)), kwargs...)
    @unpack_units D0, Ea, ΔV, P0, R = data

    @inbounds for i in eachindex(D)
        D[i] = @pow  D0 * exp(- (Ea + (P[i] - P0) * ΔV) / (R * T[i]))
    end
end

function show(io::IO, g::DiffusionData)
    return print(
        io,
        "DiffusionData: D0 = $(Value(g.D0)), Ea = $(Value(g.Ea)), ΔV = $(Value(g.ΔV))",
    )
end


export AbstractChemicalDiffusion,
       DiffusionData,
       compute_D,
       compute_D!

# load collection of chemical diffusion data
include("Data/Rutile/Rutile.jl")
using .Rutile

export Rutile


end
