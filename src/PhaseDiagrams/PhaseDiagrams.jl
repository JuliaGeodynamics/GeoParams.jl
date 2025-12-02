module PhaseDiagrams

# This contains routines to read phase diagrams from disk

using ..Units
using Unitful
using DelimitedFiles
import Base.show
using GeoParams: AbstractMaterialParam, AbstractPhaseDiagramsStruct
import GeoParams.PerpleX_LaMEM_Diagram
import GeoParams.ptr2string
using Adapt
using GeoParams: LinearInterpolator, interpolate

export PhaseDiagram_LookupTable, PerpleX_LaMEM_Diagram, MAGEMin_Diagram, MAGEMin_LookupTable

"""
    Contains data of a Phase Diagram that is regularly spaced in P & T

# Fields
- `Type::Ptr{UInt8}` : String pointer indicating the type of phase diagram (e.g., "Perple_X/MAGEMin/LaMEM")
- `Name::Ptr{UInt8}` : String pointer to the name of the phase diagram file
- `rockRho::Union{T, Nothing}` : Interpolation object for rock density
- `meltRho::Union{T, Nothing}` : Interpolation object for melt density
- `meltFrac::Union{T, Nothing}` : Interpolation object for melt fraction
- `Rho::Union{T, Nothing}` : Interpolation object for total density
- `rockVp::Union{T, Nothing}` : Interpolation object for rock P-wave velocity
- `rockVs::Union{T, Nothing}` : Interpolation object for rock S-wave velocity
- `rockVpVs::Union{T, Nothing}` : Interpolation object for rock Vp/Vs ratio
- `meltVp::Union{T, Nothing}` : Interpolation object for melt P-wave velocity
- `meltVs::Union{T, Nothing}` : Interpolation object for melt S-wave velocity
- `meltVpVs::Union{T, Nothing}` : Interpolation object for melt Vp/Vs ratio
- `Vp::Union{T, Nothing}` : Interpolation object for total P-wave velocity
- `Vs::Union{T, Nothing}` : Interpolation object for total S-wave velocity
- `VpVs::Union{T, Nothing}` : Interpolation object for total Vp/Vs ratio
- `SpecificCp::Union{T, Nothing}` : Interpolation object for specific heat capacity
- `solid_Vp::Union{T, Nothing}` : Interpolation object for solid P-wave velocity
- `solid_Vs::Union{T, Nothing}` : Interpolation object for solid S-wave velocity
- `melt_bulkModulus::Union{T, Nothing}` : Interpolation object for melt bulk modulus
- `solid_bulkModulus::Union{T, Nothing}` : Interpolation object for solid bulk modulus
- `solid_shearModulus::Union{T, Nothing}` : Interpolation object for solid shear modulus
- `Vp_uncorrected::Union{T, Nothing}` : Interpolation object for uncorrected P-wave velocity
- `Vs_uncorrected::Union{T, Nothing}` : Interpolation object for uncorrected S-wave velocity
"""
struct PhaseDiagram_LookupTable{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21} <: AbstractPhaseDiagramsStruct
    Type::Ptr{UInt8}  # using Ptr{UInt8} to avoid issues with String on GPU
    # HeaderText::Ptr{S}
    Name::Ptr{UInt8}
    rockRho::T1
    meltRho::T2
    meltFrac::T3
    Rho::T4
    rockVp::T5
    rockVs::T6
    rockVpVs::T7
    meltVp::T8
    meltVs::T9
    meltVpVs::T10
    Vp::T11
    Vs::T12
    VpVs::T13
    SpecificCp::T14
    solid_Vp::T15
    solid_Vs::T16
    melt_bulkModulus::T17
    solid_bulkModulus::T18
    solid_shearModulus::T19
    Vp_uncorrected::T20           # will hold Vs velocity corrected for pores, fluids, & melt
    Vs_uncorrected::T21
end

# Make PhaseDiagram_LookupTable adaptable for GPU arrays
Adapt.@adapt_structure PhaseDiagram_LookupTable

"""
    PD_Data = PerpleX_LaMEM_Diagram(fname::String; CharDim = nothing, type::AbstractString = "Perple_X/LaMEM")

Reads a precomputed phase diagram in the `Perple_X/MAGEMin/LaMEM` format (which is a phase diagram computed using `Perple_X/MAGEMin`, but formatted in a manner that is readable
using LaMEM or potentially any other geodynamic code).
The data is stored in the `PhaseDiagram_LookupTable` structure.

If the `CharDim` object is specified, the values of all diagrams will be non-dimensionalized.
The type kwarg is a string that indicates the type of phase diagram (default is "Perple_X/LaMEM").

# Example

```julia
julia> PD_Data = PerpleX_LaMEM_Diagram("./test_data/Peridotite.in")
Perple_X/LaMEM Phase Diagram Lookup Table:
                      File    :   ./test_data/Peridotite.in
                      T       :   293.0 - 1573.000039
                      P       :   1.0e7 - 2.9999999944e9
                      fields  :   :meltRho, :meltRho, :meltFrac, :rockRho, :Rho, :rockVp
                                  :rockVs, :rockVpVs, :meltVp, :meltVs, :meltVpVs
                                  :Vp, :Vs, :VpVs
```
Once imported, the properties on the diagram can be interpolated in a simple manner:
```
julia> PD_Data.Rho(1500,1e7)
3042.836820256982
```
This also works for vectors or arrays:
```julia
julia> T = [1500 1800; 1233 1300]
julia> P = [1e8 1e9; 1e7 1e7]
julia> rho = PD_Data.Rho.(T,P)
```
(Note the dot `.` in front of the bracket while evaluating arrays).

The fields that are available depend on what is listed in the diagram file.
The units of the fields are automatically evaluated, and employed to non-dimensionalize the parameters if `CharDim` is specified.

# Algorithm

Internally, we employ a linear 2D interpolation scheme for evaluating the phase diagram values at arbitrary (T,P) points using bilinear
interpolation. Values outside the range of the diagram are set to the boundary of the diagram. The interpolation object is directly encoded in the `PhaseDiagram_LookupTable`` object.
"""
function PerpleX_LaMEM_Diagram(fname::String; CharDim = nothing, type::AbstractString = "Perple_X/LaMEM")

    # Read header:
    #  the first 50 lines are comments (freely useable), followed by data
    n = 55
    header = open(readlines, fname)[1:55]

    name = pointer(ptr2string(fname))
    type_ptr = pointer(ptr2string(type))

    # Parse the names of the columns in the data file
    fields = split(header[49], "[")[2:end]      # this line should contain the names and units of the columns in the file

    # Throw error message if fields are not specified
    if length(fields) == 0
        error("Line 49 in the phase diagram file should contain the names of the columns")
    end

    fields_keys = Vector{Union{Missing, Symbol}}(missing, length(fields))
    fields_units = Vector{Union{Missing, Any}}(missing, length(fields))
    for i in 1:length(fields)
        # error catch
        if length(split(fields[i], ",")) != 3
            error(
                "Error parsing column: $(fields[i]). Should have 3 fields (name, symbol, units)",
            )
        end

        fields_keys[i] = Symbol(split(fields[i], ",")[1])
        unit_string = split(split(fields[i], ",")[3], "]")[1]
        fields_units[i] = uparse(unit_string)
    end

    Punit = fields_units[findfirst(i -> fields_keys[i] in [:Pressure], 1:length(fields_keys))]
    Tunit = fields_units[findfirst(i -> fields_keys[i] in [:Temperature], 1:length(fields_keys))]
    # Determine the range
    # T0 = parse(Float64, header[50]) * u"K"      # in K
    # dT = parse(Float64, header[51]) * u"K"
    T0 = parse(Float64, header[50]) * Tunit      # in K
    dT = parse(Float64, header[51]) * Tunit
    numT = parse(Int64, header[52])

    # convert to GeoUnits if necessary
    T0 = uconvert(u"K", T0)
    dT = uconvert(u"K", dT)

    P0 = parse(Float64, header[53]) * Punit    # in bar or Pa (will be convert to Pa later)
    dP = parse(Float64, header[54]) * Punit
    numP = parse(Int64, header[55])

    # convert to GeoUnits if necessary
    P0 = uconvert(u"Pa", P0)
    dP = uconvert(u"Pa", dP)
    Pmax = P0 + dP * (numP - 1)
    Tmax = T0 + dT * (numT - 1)
    # Tvec = T0:dT:(T0 + dT * (numT - 1))              # 1D vector
    # Pvec = P0:dP:(P0 + dP * (numP - 1))

    # In the LaMEM/Perple_X file format, the first 50 lines are comments
    data = readdlm(fname; skipstart = 55, header = false)        # read numerical data

    # Shape of 2D arrays:
    siz = (numT, numP)

    # Initialize fields in the order they are defined in the PhaseDiagram_LookupTable structure
    Struct_Fieldnames = fieldnames(PhaseDiagram_LookupTable)[3:end] # fieldnames from structure

    # Process all fields that are present in the phase diagram (and non-dimensionalize if requested)
    Struct_Fields = Vector{Union{Nothing, LinearInterpolator}}(
        nothing, length(Struct_Fieldnames)
    )

    # Loop through all fields and create an interpolation object for each (which is optionally non-dimensionalized)
    for (i, field) in enumerate(fields_keys)
        ind = findall(Struct_Fieldnames .== field)
        if length(ind) > 0
            Struct_Fields[ind[1]] = CreateInterpolationObject_PhaseDiagram(
                data[:, i], T0, dT, numT, Tmax, P0, dP, numP, Pmax, siz, fields_units[i], CharDim
            )
        end
    end

    # Some fields have melt and solid part; we can reconstruct the total part as an arithmetic average:
    Struct_Fields = ComputeTotalField_withMeltFraction(
        :Rho, :meltRho, :rockRho, :meltFrac, Struct_Fields, Struct_Fieldnames
    )
    Struct_Fields = ComputeTotalField_withMeltFraction(
        :Vp, :meltVp, :rockVp, :meltFrac, Struct_Fields, Struct_Fieldnames
    )
    Struct_Fields = ComputeTotalField_withMeltFraction(
        :Vs, :meltVs, :rockVs, :meltFrac, Struct_Fields, Struct_Fieldnames
    )
    Struct_Fields = ComputeTotalField_withMeltFraction(
        :VpVs, :meltVpVs, :rockVpVs, :meltFrac, Struct_Fields, Struct_Fieldnames
    )

    # Store in phase diagram structure
    PD_data = PhaseDiagram_LookupTable(
        type_ptr, name, Struct_Fields...
    )

    return PD_data
end

# Print info
function show(io::IO, d::AbstractPhaseDiagramsStruct)
    # T = d.rockRho.knots[1]
    # P = d.rockRho.knots[2]

    println(io, "$(ptr2string(d.Type)) Phase Diagram Lookup Table: ")
    println(io, "                      File    :   $(ptr2string(d.Name))")
    println(io, "                      T       :   $(d.rockRho.T0) - $(d.rockRho.Tmax)")
    println(io, "                      P       :   $(d.rockRho.P0) - $(d.rockRho.Pmax)")

    lst = fieldnames(typeof(d))
    str = ""
    for i in 4:length(lst)
        if !isnothing(getfield(d, lst[i]))
            if i == 4
                str = "                      fields  :   :$(lst[i]),"
            else
                str = str * " :$(lst[i]),"
            end
            if mod((i - 3), 5) == 0
                str = str[1:(end - 1)]
                println(io, str)
                str = "                                 "
            end
        end
    end
    str = str[1:(end - 1)]
    return println(io, str)
end

# Internal routine that creates an interpolation object from a column of the data
function CreateInterpolationObject_PhaseDiagram(
        data_vec::AbstractArray{Float64, 1}, T0, dT, numT, Tmax, P0, dP, numP, Pmax, siz::Tuple{Int64, Int64}, units, CharDim
    )
    data_units = reshape(data_vec, siz) * units      # Create 2D array

    # # Convert to Pa & K
    # Pvec_Pa = Float64.(uconvert.(u"Pa", Pvec))
    # Tvec_K = Float64.(uconvert.(u"K", Tvec))

    # Optional: nondimensionalize values as well as T and P
    if CharDim == nothing
        T0 = ustrip(T0)
        dT = ustrip(dT)
        P0 = ustrip(P0)
        dP = ustrip(dP)
        Tmax = ustrip(Tmax)
        Pmax = ustrip(Pmax)
        data = ustrip.(data_units)
    else
        P0 = nondimensionalize(P0, CharDim)
        dP = nondimensionalize(dP, CharDim)
        T0 = nondimensionalize(T0, CharDim)
        dT = nondimensionalize(dT, CharDim)
        Tmax = nondimensionalize(Tmax, CharDim)
        Pmax = nondimensionalize(Pmax, CharDim)

        # data_tmp = data_units isa Matrix{<:GeoUnit} ? data_units : GeoUnit.(data_units)
        # data = nondimensionalize(data_tmp, CharDim)
        data = if data_units isa Matrix{<:Quantity}
            nondimensionalize(data_units, CharDim)
        else
            data_units
        end
    end

    # Create interpolation object
    intp_data = interpolate(T0, dT, numT, Tmax, P0, dP, numP, Pmax, data)

    return intp_data
end

# Internal routine, which combines solid & melt properties with melt fraction
function ComputeTotalField_withMeltFraction(
        totalData::Symbol,
        meltData::Symbol,
        solidData::Symbol,
        meltFrac::Symbol,
        Struct_Fields,
        Struct_Fieldnames,
    )
    ind_totalData = findall(Struct_Fieldnames .== totalData)
    ind_meltData = findall(Struct_Fieldnames .== meltData)
    ind_solidData = findall(Struct_Fieldnames .== solidData)
    ind_meltFrac = findall(Struct_Fieldnames .== meltFrac)
    # extract arrays from interpolation objects
    if (length(ind_totalData) > 0) &
            (length(ind_meltData) > 0) &
            (length(ind_solidData) > 0) &
            (length(ind_meltFrac) > 0) &
            !(isnothing(Struct_Fields[ind_meltFrac[1]])) &
            !(isnothing(Struct_Fields[ind_solidData[1]])) &
            !(isnothing(Struct_Fields[ind_meltData[1]]))
        if Struct_Fields[ind_meltFrac[1]] != nothing
            ϕ = Struct_Fields[ind_meltFrac[1]].coefs      # melt fraction
            S = Struct_Fields[ind_solidData[1]].coefs      # solid property
            M = Struct_Fields[ind_meltData[1]].coefs      # melt property

            # Compute average
            Result = (1.0 .- ϕ) .* S .+ ϕ .* M

            # Create interpolation object of average
            # Tvec = Struct_Fields[ind_meltFrac[1]].knots[1]       # temperature data
            # Pvec = Struct_Fields[ind_meltFrac[1]].knots[2]       # pressure data
            T0_new = Struct_Fields[ind_meltFrac[1]].T0
            dT_new = Struct_Fields[ind_meltFrac[1]].dT
            numT_new = Struct_Fields[ind_meltFrac[1]].numT
            Tmax_new = Struct_Fields[ind_meltFrac[1]].Tmax
            P0_new = Struct_Fields[ind_meltFrac[1]].P0
            dP_new = Struct_Fields[ind_meltFrac[1]].dP
            numP_new = Struct_Fields[ind_meltFrac[1]].numP
            Pmax_new = Struct_Fields[ind_meltFrac[1]].Pmax
            # intp_Result = interpolate((Tvec, Pvec), Result)
            intp_Result = interpolate(T0_new, dT_new, numT_new, Tmax_new, P0_new, dP_new, numP_new, Pmax_new, Result)

            # assign
            Struct_Fields[ind_totalData[1]] = intp_Result
        else
            # use solid fraction
            S = Struct_Fields[ind_solidData[1]]      # solid property

            # assign
            Struct_Fields[ind_totalData[1]] = S
        end
    end

    return Struct_Fields
end

struct MAGEMin_LookupTable{T} <: AbstractPhaseDiagramsStruct
    Type::Ptr{UInt8}  # using Ptr{UInt8} to avoid issues with String on GPU
    Name::Ptr{UInt8}
    meltRho::T
    meltFrac::T
    rockRho::T
    Vp::T
    Vs::T
    SpecificCp::T
    Rho::T
end

Adapt.@adapt_structure MAGEMin_LookupTable

function MAGEMin_Diagram(fname::String; CharDim = nothing, type::AbstractString = "MAGEMin")

    # Read header:
    #  the first 50 lines are comments (freely useable), followed by data
    n = 55
    header = open(readlines, fname)[1:55]

    name = pointer(ptr2string(fname))
    type_ptr = pointer(ptr2string(type))

    # Parse the names of the columns in the data file
    fields = split(header[49], "[")[2:end]      # this line should contain the names and units of the columns in the file

    # Throw error message if fields are not specified
    if length(fields) == 0
        error("Line 49 in the phase diagram file should contain the names of the columns")
    end

    fields_keys = Vector{Symbol}(undef, length(fields))
    fields_units = Vector{Unitful.Units}(undef, length(fields))

    for i in 1:length(fields)
        # error catch
        if length(split(fields[i], ",")) != 3
            error(
                "Error parsing column: $(fields[i]). Should have 3 fields (name, symbol, units)",
            )
        end

        fields_keys[i] = Symbol(split(fields[i], ",")[1])
        unit_string = split(split(fields[i], ",")[3], "]")[1]
        fields_units[i] = uparse(unit_string)
    end

    Punit = fields_units[findfirst(i -> fields_keys[i] in [:Pressure], 1:length(fields_keys))]
    Tunit = fields_units[findfirst(i -> fields_keys[i] in [:Temperature], 1:length(fields_keys))]
    # Determine the range
    T0 = parse(Float64, header[50]) * Tunit      # in K
    dT = parse(Float64, header[51]) * Tunit
    numT = parse(Int64, header[52])

    # convert to GeoUnits if necessary
    T0 = uconvert(u"K", T0)
    dT = uconvert(u"K", dT)

    P0 = parse(Float64, header[53]) * Punit    # in bar or Pa (will be convert to Pa later)
    dP = parse(Float64, header[54]) * Punit
    numP = parse(Int64, header[55])

    # convert to GeoUnits if necessary
    P0 = uconvert(u"Pa", P0)
    dP = uconvert(u"Pa", dP)
    Pmax = P0 + dP * (numP - 1)
    Tmax = T0 + dT * (numT - 1)

    data = readdlm(fname; skipstart = 55, header = false)        # read numerical data

    # Shape of 2D arrays:
    siz = (numT, numP)

    # # Initialize fields in the order they are defined in the PhaseDiagram_LookupTable structure
    Struct_Fieldnames = fieldnames(MAGEMin_LookupTable)[3:end] # fieldnames from structure

    Struct_Fields = Vector{LinearInterpolator{Float64, Int64, Array{Float64,2}}}(
        undef, length(Struct_Fieldnames)
    )

    # for i in 1:length(Struct_Fieldnames)
    #     field = fields_keys[i]
    #     ind = findall(Struct_Fieldnames .== field)
    #     if length(ind) > 0
    #         Struct_Fields[ind[1]] = CreateInterpolationObject_PhaseDiagram(
    #             data[:, i], T0, dT, numT, Tmax, P0, dP, numP, Pmax, siz, fields_units[i], CharDim
    #         )
    #     end
    # end
    for i in 1:6
        Struct_Fields[i] = CreateInterpolationObject_PhaseDiagram(
            data[:, i], T0, dT, numT, Tmax, P0, dP, numP, Pmax, siz, fields_units[i], CharDim
        )
    end

    # Some fields have melt and solid part; we can reconstruct the total part as an arithmetic average:
    Struct_Fields = ComputeTotalField_withMeltFraction(
        :Rho, :meltRho, :rockRho, :meltFrac, Struct_Fields, Struct_Fieldnames
    )

    PD_data = MAGEMin_LookupTable(
        type_ptr, name, Struct_Fields...
    )

    return PD_data
end

end
