export DiffusionData, AbstractChemicalDiffusion

abstract type AbstractChemicalDiffusion{T} <: AbstractMaterialParam end

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
struct DiffusionData{T, U1, U2, U3, U4, U5} <: AbstractChemicalDiffusion{T}
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
    T_range::Tuple{GeoUnit{T, U4}, GeoUnit{T, U4}} # temperature range
    P0::GeoUnit{T, U5} # pressure of calibration
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
        T_rangeU = convert.(GeoUnit, T_range)
        P0U = convert(GeoUnit, P0)
        # Extract struct types
        T = typeof(D0U).types[1]
        U1 = typeof(D0U).types[2]
        U2 = typeof(EaU).types[2]
        U3 = typeof(ΔVU).types[2]
        U4 = typeof(T_rangeU[1]).types[2]
        U5 = typeof(P0U).types[2]
        name = pointer(ptr2string(Name))
        mineral = pointer(ptr2string(Mineral))
        formula = pointer(ptr2string(Formula))
        species = pointer(ptr2string(Species))
        orientation = pointer(ptr2string(Orientation))
        crystallography = pointer(ptr2string(Crystallography))
        buffer = pointer(ptr2string(Buffer))
        fluid = pointer(ptr2string(Fluid))

        # Create struct
        return new{T,U1,U2,U3,U4,U5}(
            name, mineral, formula, species, D0U, D0_2σU, EaU, Ea_2σU, ΔVU, ΔV_2σU, T_rangeU, P0U, orientation, crystallography, buffer, fluid
        )
    end
end

# load collection of chemical diffusion data
include("Data/Rutile/Rutile.jl")
using .Rutile
