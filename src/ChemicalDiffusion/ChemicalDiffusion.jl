module ChemicalDiffusion

using GeoParams
using Unitful

using Parameters, LaTeXStrings, Unitful, MuladdMacro
using ..Units
using GeoParams: AbstractMaterialParam, AbstractConstitutiveLaw, AbstractComposite
import GeoParams: param_info, ptr2string
using ..MaterialParameters: MaterialParamsInfo
import Base.show

# Exported functions defined in this module
export SetChemicalDiffusion,
    Transform_ChemicalDiffusion,
    AbstractChemicalDiffusion,
    DiffusionData,
    compute_D,
    compute_D!

# load collection of chemical diffusion data
include("Data/Rutile/Rutile.jl")
using .Rutile
include("Data/Garnet/Garnet.jl")
using .Garnet
include("Data/Olivine/Olivine.jl")
using .Olivine
include("Data/Melt/Melt.jl")
using .Melt

# Exported modules of chemical diffusion data
export Rutile, Garnet, Olivine, Melt

abstract type AbstractChemicalDiffusion{T} <: AbstractMaterialParam end

@inline diffusion_database(f::F) where {F} = first(f())
@inline diffusion_database_info(f::F) where {F} = last(f())

@inline precision(::AbstractChemicalDiffusion{T}) where {T} = T

"""
    DiffusionData(; Name, Phase, Formula, Species, Orientation, Crystallography, Buffer, Fluid, Doping, D0, log_D0_1σ, Ea, Ea_1σ, ΔV, ΔV_1σ, afO2, bfO2, nfO2, Charge, T_range, P0, )

Defines the diffusion data for the chemical diffusion of a given phase and species from an experiment.

The diffusion coefficient `D` [\\mathrm{[m^2/s]}] is given by an Arrhenius equation:
```math
    D = D0 * (afO2 * (fO2^nfO2 / bfO2)) * \\exp\\left(-\\frac{Ea + PΔV} {RT}\\right)
```
where
- ``D0`` is the pre-exponential factor [\\mathrm{[m^2/s]}],
- ``fO2`` is the oxygen fugacity,
- ``afO2`` is the prefactor for the fO2 dependency,
- ``bfO2`` is the quotient for the fO2 dependency,
- ``nfO2`` is the exponent for the fO2 dependency,
- ``Ea`` is the activation energy [\\mathrm{[J/mol]}],
- ``ΔV`` is the activation volume [\\mathrm{[cm^3/mol]}],
- ``P`` is the pressure [\\mathrm{[Pa]},
- ``T`` is the temperature [\\mathrm{[K]}],
- ``R`` is the gas constant [\\mathrm{[J/(mol K)]}],
"""
struct DiffusionData{T, U1, U2, U3, U4, U5, U6, U7, U8, U9, U10} <: AbstractChemicalDiffusion{T}
    Name::Ptr{UInt8}  # name of the diffusion experiment and paper
    Phase::Ptr{UInt8}  # name of the phase
    Formula::Ptr{UInt8}  # chemical formula of the mineral
    Species::Ptr{UInt8}  # element or species being diffused
    Orientation::Ptr{UInt8}  # Crystal orientation from the diffusion experiment
    Crystallography::Ptr{UInt8}  # Crystallographic system of the mineral
    Buffer::Ptr{UInt8}  # Buffer condition (e.g., NNO) during the experiment
    Fluid::Ptr{UInt8}  # Fluid condition (e.g., anhydrous) during the experiment
    Doping::Ptr{UInt8}  # Doping condition (e.g., Cr2O3) during the experiment
    D0::GeoUnit{T, U1}  # pre-exponential factor
    log_D0_1σ::GeoUnit{T, U2}  # uncertainty at 1σ of the pre-exponential factor
    Ea::GeoUnit{T, U3}  # activation energy
    Ea_1σ::GeoUnit{T, U4}  # uncertainty at 1σ of the activation energy
    ΔV::GeoUnit{T, U5}  # activation volume
    ΔV_1σ::GeoUnit{T, U6}  # uncertainty at 1σ of the activation volume
    afO2::GeoUnit{T, U2}  # prefactor for f(O2) dependency
    bfO2::GeoUnit{T, U2}  # quotient for f(O2) dependency
    nfO2::GeoUnit{T, U2}  # exponent for f(O2) dependency
    Charge::GeoUnit{T, U7}  # charge of the species
    R::GeoUnit{T, U8}  # gas constant
    T_range_min::GeoUnit{T, U9}  # minimum temperature of the T_range
    T_range_max::GeoUnit{T, U9}  # maximum temperature of the T_range
    P0::GeoUnit{T, U10}  # pressure of calibration

    function DiffusionData(;
            Name = "Unknown",  # name of the diffusion experiment and paper
            Phase = "Unknown",  # name of the mineral
            Formula = "Unknown",  # chemical formula of the mineral
            Species = "Unknown",  # element or species being diffused
            Orientation = "Unknown",  # Crystal orientation from the diffusion experiment
            Crystallography = "Unknown",  # Crystallographic system of the mineral
            Buffer = "Unknown",  # Buffer condition (e.g., NNO) during the experiment
            Fluid = "Unknown",  # Fluid condition (e.g., anhydrous) during the experiment
            Doping = "Pure",  # Doping condition (e.g., Cr) during the experiment
            D0 = 0.0m^2 / s,  # pre-exponential factor
            log_D0_1σ = 0.0NoUnits,  # uncertainty at 1σ of the pre-exponential factor in the log form
            Ea = 0.0J / mol,  # activation energy
            Ea_1σ = 0.0J / mol,  # uncertainty at 1σ of the activation energy
            ΔV = 0.0cm^3 / mol,  # activation volume
            ΔV_1σ = 0.0cm^3 / mol,  # uncertainty at 1σ of the activation volume
            afO2 = 1.0NoUnits,  # prefactor for f(O2) dependency
            bfO2 = 1.0NoUnits,  # quotient for f(O2) dependency
            nfO2 = 1.0NoUnits,  # exponent for f(O2) dependency
            Charge = 0NoUnits,  # charge of the species
            R = Unitful.R,  # gas constant
            T_range_min = 0.0K,  # minimum temperature of the T_range
            T_range_max = 0.0K,  # maximum temperature of the T_range
            P0 = 0.0Pa,  # pressure of calibration
        )

        # Convert to GeoUnits
        D0U = convert(GeoUnit, D0)
        log_D0_1σU = convert(GeoUnit, log_D0_1σ)
        EaU = convert(GeoUnit, Ea)
        Ea_1σU = convert(GeoUnit, Ea_1σ)
        ΔVU = convert(GeoUnit, ΔV)
        ΔV_1σU = convert(GeoUnit, ΔV_1σ)
        RU = convert(GeoUnit, R)
        afO2U = convert(GeoUnit, afO2)
        bfO2U = convert(GeoUnit, bfO2)
        nfO2U = convert(GeoUnit, nfO2)
        ChargeU = convert(GeoUnit, Charge)
        T_range_minU = convert(GeoUnit, T_range_min)
        T_range_maxU = convert(GeoUnit, T_range_max)
        P0U = convert(GeoUnit, P0)
        # Extract struct types
        T = typeof(D0U).types[1]
        U1 = typeof(D0U).types[2]
        U2 = typeof(log_D0_1σU).types[2]
        U3 = typeof(EaU).types[2]
        U4 = typeof(Ea_1σU).types[2]
        U5 = typeof(ΔVU).types[2]
        U6 = typeof(ΔV_1σU).types[2]
        U7 = typeof(ChargeU).types[2]
        U8 = typeof(RU).types[2]
        U9 = typeof(T_range_minU).types[2]
        U10 = typeof(P0U).types[2]
        name = pointer(ptr2string(Name))
        mineral = pointer(ptr2string(Phase))
        formula = pointer(ptr2string(Formula))
        species = pointer(ptr2string(Species))
        orientation = pointer(ptr2string(Orientation))
        crystallography = pointer(ptr2string(Crystallography))
        buffer = pointer(ptr2string(Buffer))
        fluid = pointer(ptr2string(Fluid))
        doping = pointer(ptr2string(Doping))

        # Create struct
        return new{T, U1, U2, U3, U4, U5, U6, U7, U8, U9, U10}(
            name, mineral, formula, species, orientation, crystallography, buffer, fluid, doping, D0U, log_D0_1σU, EaU, Ea_1σU, ΔVU, ΔV_1σU, afO2U, bfO2U, nfO2U, ChargeU, RU, T_range_minU, T_range_maxU, P0U
        )
    end
end

function DiffusionData(
        Name, Phase, Formula, Species, Orientation, Crystallography, Buffer, Fluid, Doping, D0, log_D0_1σ, Ea, Ea_1σ, ΔV, ΔV_1σ, afO2, bfO2, nfO2, Charge, R, T_range_min, T_range_max, P0,
    )
    return DiffusionData(;
        Name = Name,
        Phase = Phase,
        Formula = Formula,
        Species = Species,
        Orientation = Orientation,
        Crystallography = Crystallography,
        Buffer = Buffer,
        Fluid = Fluid,
        Doping = Doping,
        D0 = D0,
        log_D0_1σ = log_D0_1σ,
        Ea = Ea,
        Ea_1σ = Ea_1σ,
        ΔV = ΔV,
        ΔV_1σ = ΔV_1σ,
        afO2 = afO2,
        bfO2 = bfO2,
        nfO2 = nfO2,
        R = R,
        Charge = Charge,
        T_range_min = T_range_min,
        T_range_max = T_range_max,
        P0 = P0,
    )
end

function param_info(data::DiffusionData) # info about the struct
    return MaterialParamsInfo(;
        Equation = L"D = D0 * (afO2 * fO2^(nfO2) / bfO2) * \exp\left(-\frac{Ea + PΔV} {RT}\right)",
    )
end


"""
    compute_D(data::DiffusionData; T=1K, P=1GPa, fO2 = 1NoUnits, kwargs...)

Computes the diffusion coefficient `D` [m^2/s] from the diffusion data `data` at temperature `T` [K], pressure `P` [Pa] and oxygen fugacity `fO2` [NoUnits] from a structure of type `DiffusionData`.
If `T` and `P` are provided without unit, the function assumes the units are in Kelvin and Pascal, respectively, and outputs the diffusion coefficient without unit based on the value in m^2/s.
"""
@inline function compute_D(data::DiffusionData; T = 1K, P = 1GPa, fO2 = 1NoUnits, kwargs...)

    if P isa Quantity && T isa Quantity
        @unpack_units D0, Ea, ΔV, P0, R, afO2, bfO2, nfO2 = data

        # convert to K to prevent affine error with Celsius
        T = uconvert(K, T)
    else
        @unpack_val D0, Ea, ΔV, P0, R, afO2, bfO2, nfO2 = data
    end

    D = @muladd D0 * (afO2 * (fO2^nfO2) / bfO2) * exp(-(Ea + (P - P0) * ΔV) / (R * T))

    return D
end

"""
    compute_D!(D, data::DiffusionData; T=ones(size(D))K, P=zeros(size(D))Pa, fO2 = ones(size(D)), kwargs...)

In-place version of `compute_D(data::DiffusionData; T=1K, P=1GPa, fO2=0NoUnits, kwargs...)`. `D` should be an array of the same size as T, P and fO2.
"""
function compute_D!(
        D::AbstractArray{_T, nDim},
        data::DiffusionData;
        T = ones(size(D))K,
        P = ones(size(D))GPa,
        fO2 = ones(size(D)),
        kwargs...
    ) where {_T, nDim}

    return @inbounds for i in eachindex(D)
        D[i] = compute_D(data; T = T[i], P = P[i], fO2 = fO2[i], kwargs...)
    end
end


function show(io::IO, g::DiffusionData)
    return print(
        io,
        "DiffusionData: Phase = $(unsafe_string(g.Phase)), Species = $(unsafe_string(g.Species)), D0 = $(Value(g.D0)), Ea = $(Value(g.Ea)), ΔV = $(Value(g.ΔV))",
    )
end


"""
    SetChemicalDiffusion["Name of Chemical Diffusion"]
"""
function SetChemicalDiffusion(
        name::F;
        D0 = nothing,
        log_D0_1σ = nothing,
        Ea = nothing,
        Ea_1σ = nothing,
        ΔV = nothing,
        ΔV_1σ = nothing,
        afO2 = nothing,
        bfO2 = nothing,
        nfO2 = nothing,
        Charge = nothing,
        T_range_min = nothing,
        T_range_max = nothing,
        P0 = nothing,
    ) where {F}
    kwargs = (; D0, log_D0_1σ, Ea, Ea_1σ, ΔV, ΔV_1σ, afO2, bfO2, nfO2, Charge, T_range_min, T_range_max, P0)
    return Transform_ChemicalDiffusion(name, kwargs)
end

function SetChemicalDiffusion(
        name::F,
        CharDim::GeoUnits{T};
        D0 = nothing,
        log_D0_1σ = nothing,
        Ea = nothing,
        Ea_1σ = nothing,
        ΔV = nothing,
        ΔV_1σ = nothing,
        afO2 = nothing,
        bfO2 = nothing,
        nfO2 = nothing,
        Charge = nothing,
        T_range_min = nothing,
        T_range_max = nothing,
        P0 = nothing,
    ) where {F, T <: Union{GEO, SI}}
    kwargs = (; D0, log_D0_1σ, Ea, Ea_1σ, ΔV, ΔV_1σ, afO2, bfO2, nfO2, Charge, T_range_min, T_range_max, P0)
    return nondimensionalize(Transform_ChemicalDiffusion(name, kwargs), CharDim)
end

"""
    Transform_ChemicalDiffusion(name)
Transforms units from MPa, kJ etc. to basic units such as Pa, J etc.
"""
Transform_ChemicalDiffusion(name::F) where {F} = Transform_ChemicalDiffusion(diffusion_database(name))
Transform_ChemicalDiffusion(name::F, kwargs::NamedTuple) where {F} = Transform_ChemicalDiffusion(diffusion_database(name), kwargs)

function Transform_ChemicalDiffusion(name::F, CharDim::GeoUnits{U}) where {F, U <: Union{GEO, SI}}
    return Transform_ChemicalDiffusion(diffusion_database(name), CharDim)
end

function Transform_ChemicalDiffusion(p::AbstractChemicalDiffusion{T}, CharDim::GeoUnits{U}) where {T, U <: Union{GEO, SI}}
    return nondimensionalize(Transform_ChemicalDiffusion(p), CharDim)
end

function Transform_ChemicalDiffusion(pp::AbstractChemicalDiffusion)

    D0 = Value(pp.D0)
    log_D0_1σ = Value(pp.log_D0_1σ)
    Ea = Value(pp.Ea)
    Ea_1σ = Value(pp.Ea_1σ)
    ΔV = Value(pp.ΔV)
    ΔV_1σ = Value(pp.ΔV_1σ)
    afO2 = Value(pp.afO2)
    bfO2 = Value(pp.bfO2)
    nfO2 = Value(pp.nfO2)
    Charge = Value(pp.Charge)
    T_range_min = Value(pp.T_range_min)
    T_range_max = Value(pp.T_range_max)
    P0 = Value(pp.P0)

    D0_SI = uconvert(m^2 / s, D0)
    Ea_SI = uconvert(J / mol, pp.Ea)
    Ea_1σ_SI = uconvert(J / mol, pp.Ea_1σ)
    ΔV_SI = uconvert(m^3 / mol, pp.ΔV)
    ΔV_1σ_SI = uconvert(m^3 / mol, pp.ΔV_1σ)
    T_range_min_SI = uconvert(K, T_range_min)
    T_range_max_SI = uconvert(K, T_range_max)
    P0_SI = uconvert(Pa, P0)


    args = (
        Name = unsafe_string(pp.Name),
        Phase = unsafe_string(pp.Phase),
        Formula = unsafe_string(pp.Formula),
        Species = unsafe_string(pp.Species),
        Orientation = unsafe_string(pp.Orientation),
        Crystallography = unsafe_string(pp.Crystallography),
        Buffer = unsafe_string(pp.Buffer),
        Fluid = unsafe_string(pp.Fluid),
        Doping = unsafe_string(pp.Doping),
        D0 = D0_SI, log_D0_1σ = log_D0_1σ,
        Ea = Ea_SI, Ea_1σ = Ea_1σ_SI,
        ΔV = ΔV_SI, ΔV_1σ = ΔV_1σ_SI,
        afO2 = afO2,
        bfO2 = bfO2,
        nfO2 = nfO2,
        Charge = Charge,
        T_range_min_SI = T_range_min_SI,
        T_range_max_SI = T_range_max_SI,
        P0 = P0_SI,
    )

    return DiffusionData(; args...)
end

function Transform_ChemicalDiffusion(pp::AbstractChemicalDiffusion, kwargs::NamedTuple)

    f(a, b) = Value(GeoUnit(a))
    f(::Nothing, b) = Value(b)

    (; D0, log_D0_1σ, Ea, Ea_1σ, ΔV, ΔV_1σ, afO2, bfO2, nfO2, Charge, T_range_min, T_range_max, P0) = kwargs

    D0_new = f(D0, pp.D0)
    log_D0_1σ_new = f(log_D0_1σ, pp.log_D0_1σ)
    Ea_new = f(Ea, pp.Ea)
    Ea_1σ_new = f(Ea_1σ, pp.Ea_1σ)
    ΔV_new = f(ΔV, pp.ΔV)
    ΔV_1σ_new = f(ΔV_1σ, pp.ΔV_1σ)
    afO2_new = f(afO2, pp.afO2)
    bfO2_new = f(bfO2, pp.bfO2)
    nfO2_new = f(nfO2, pp.nfO2)
    Charge_new = f(Charge, pp.Charge)
    T_range_min_new = f(T_range_min, pp.T_range_min)
    T_range_max_new = f(T_range_max, pp.T_range_max)
    P0_new = f(P0, pp.P0)

    D0_SI = uconvert(m^2 / s, D0_new)
    Ea_SI = uconvert(J / mol, Ea_new)
    Ea_1σ_SI = uconvert(J / mol, Ea_1σ_new)
    ΔV_SI = uconvert(m^3 / mol, ΔV_new)
    ΔV_1σ_SI = uconvert(m^3 / mol, ΔV_1σ_new)
    T_range_min_SI = uconvert(K, T_range_min_new)
    T_range_max_SI = uconvert(K, T_range_max_new)
    P0_SI = uconvert(Pa, P0_new)

    args = (
        Name = unsafe_string(pp.Name),
        Phase = unsafe_string(pp.Phase),
        Formula = unsafe_string(pp.Formula),
        Species = unsafe_string(pp.Species),
        Orientation = unsafe_string(pp.Orientation),
        Crystallography = unsafe_string(pp.Crystallography),
        Buffer = unsafe_string(pp.Buffer),
        Fluid = unsafe_string(pp.Fluid),
        Doping = unsafe_string(pp.Doping),
        D0 = D0_SI, log_D0_1σ = log_D0_1σ_new,
        Ea = Ea_SI, Ea_1σ = Ea_1σ_SI,
        ΔV = ΔV_SI, ΔV_1σ = ΔV_1σ_SI,
        afO2 = afO2_new,
        bfO2 = bfO2_new,
        nfO2 = nfO2_new,
        Charge = Charge_new,
        T_range_min = T_range_min_SI,
        T_range_max = T_range_max_SI,
        P0 = P0_SI,
    )

    return DiffusionData(; args...)
end


end
