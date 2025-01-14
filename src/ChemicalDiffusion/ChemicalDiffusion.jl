module ChemicalDiffusion


using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, @extractors, add_extractor_functions, ptr2string, GeoUnit
using ..MaterialParameters: MaterialParamsInfo
import Base.show, GeoParams.param_info

abstract type AbstractChemicalDiffusion{T} <: AbstractMaterialParam end

"""
    DiffusionData(; Name, Mineral, Formula, Species, D0, E, V, T_range, P0, Orientation, Crystallography, Buffer)

Defines the diffusion data for a given mineral and species.

The diffusion coefficient `D` [\\mathrm{[m^2/s]}] is given by an Arrhenius equation:
```math
    D = D0 * \\exp\\left(-\\frac{E + PV} {RT}\\right)
```
where
- ``D0`` is the pre-exponential factor [\\mathrm{[m^2/s]}],
- ``E`` is the activation energy [\\mathrm{[J/mol]}],
- ``V`` is the activation volume [\\mathrm{[cm^3/mol]}],
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
    E::GeoUnit{T, U2} # activation energy
    V::GeoUnit{T, U3} # activation volume
    T_range::Tuple{GeoUnit{T, U4}, GeoUnit{T, U4}} # temperature range
    P0::GeoUnit{T, U5} # effective pressure
    Crystal_orientation::Ptr{UInt8}  # Crystal orientation from the diffusion experiment
    Crystallography::Ptr{UInt8}  # Crystallographic system of the mineral
    Buffer::Ptr{UInt8}  # Buffer condition (e.g., NNO) during the experiment

    function DiffusionData(;
        Name="",  # name of the diffusion experiment and paper
        Mineral="",  # name of the mineral
        Formula="",  # chemical formula of the mineral
        Species="",  # element or species being diffused
        D0=0.0m^2/s,  # pre-exponential factor
        E=0.0J/mol,  # activation energy
        V=0.0cm^3/mol,  # activation volume
        T_range=(0.0K, 0.0K),  # temperature range
        P0=0.0Pa,  # effective pressure
        Orientation="",  # Crystal orientation from the diffusion experiment
        Crystallography="",  # Crystallographic system of the mineral
        Buffer=""  # Buffer condition (e.g., NNO) during the experiment
    )

        # Convert to GeoUnits
        D0U = convert(GeoUnit, D0)
        EU = convert(GeoUnit, E)
        VU = convert(GeoUnit, V)
        T_rangeU = convert.(GeoUnit, T_range)
        P0U = convert(GeoUnit, P0)
        # Extract struct types
        T = typeof(D0U).types[1]
        U1 = typeof(D0U).types[2]
        U2 = typeof(EU).types[2]
        U3 = typeof(VU).types[2]
        U4 = typeof(T_rangeU[1]).types[2]
        U5 = typeof(P0U).types[2]
        name = pointer(ptr2string(Name))
        mineral = pointer(ptr2string(Mineral))
        formula = pointer(ptr2string(Formula))
        species = pointer(ptr2string(Species))
        orientation = pointer(ptr2string(Orientation))
        crystallography = pointer(ptr2string(Crystallography))
        buffer = pointer(ptr2string(Buffer))

        # Create struct
        return new{T,U1,U2,U3,U4,U5}(
            name, mineral, formula, species, D0U, EU, VU, T_rangeU, P0U, orientation, crystallography, buffer
        )
    end

end


end