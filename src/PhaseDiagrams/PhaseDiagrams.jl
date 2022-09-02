module PhaseDiagrams

# This contains routines to read phase diagrams from disk 

using ..Units
using Unitful
using DelimitedFiles, Interpolations
import Base.show
using GeoParams: AbstractMaterialParam, AbstractPhaseDiagramsStruct
import GeoParams.PerpleX_LaMEM_Diagram

export PhaseDiagram_LookupTable, PerpleX_LaMEM_Diagram, ComputeDensity

"""
    Contains data of a Phase Diagram that is regularly spaced in P & T
"""
struct PhaseDiagram_LookupTable{S,T,nothing} <: AbstractPhaseDiagramsStruct
    Type::S
    HeaderText::Vector{S}
    Name::S
    rockRho::Union{T,nothing}
    meltRho::Union{T,nothing}
    meltFrac::Union{T,nothing}
    Rho::Union{T,nothing}
    rockVp::Union{T,nothing}
    rockVs::Union{T,nothing}
    rockVpVs::Union{T,nothing}
    meltVp::Union{T,nothing}
    meltVs::Union{T,nothing}
    meltVpVs::Union{T,nothing}
    Vp::Union{T,nothing}
    Vs::Union{T,nothing}
    VpVs::Union{T,nothing}
    cpxFrac::Union{T,nothing}
    solid_Vp::Union{T,nothing}
    solid_Vs::Union{T,nothing}
    melt_bulkModulus::Union{T,nothing}
    solid_bulkModulus::Union{T,nothing}
    solid_shearModulus::Union{T,nothing}
    Vp_uncorrected::Union{T,nothing}           # will hold Vs velocity corrected for pores, fluids, & melt 
    Vs_uncorrected::Union{T,nothing}
end

"""
    PD_Data = PerpleX_LaMEM_Diagram(fname::String; CharDim = nothing)

Reads a precomputed phase diagram in the `LaMEM/Perple_X` format (which is a phase diagram computed using `Perple_X`, but formatted in a manner that is readable using LaMEM).
The data is stored in the `PhaseDiagram_LookupTable` structure.

If the `CharDim` object is specified, the values of all diagrams will be non-dimensionalized.

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

Internally, we employ linear interpolation, as provided by the [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) package.
Values outside the range of the diagram are set to the boundary of the diagram. The interpolation object is directly encoded in the `PhaseDiagram_LookupTable`` object.  

"""
function PerpleX_LaMEM_Diagram(fname::String; CharDim=nothing)

    # Read header: 
    #  the first 50 lines are comments (freely useable), followed by data 
    n = 55
    header = open(readlines, `head -n $(n) $(fname)`)
    header_text = header[1:49]

    # Parse the names of the collumns in the data file 
    fields = split(header[49], "[")[2:end]      # this line should contain the names and units of the collumns in the file 

    # Throw error message if fields are not specified
    if length(fields) == 0
        error("Line 49 in the phase diagram file should contain the names of the collums")
    end

    fields_keys = Vector{Union{Missing,Symbol}}(missing, length(fields))
    fields_units = Vector{Union{Missing,Any}}(missing, length(fields))
    for i in 1:length(fields)
        # error catch
        if length(split(fields[i], ",")) != 3
            error(
                "Error parsing collumn: $(fields[i]). Should have 3 fields (name, symbol, units)",
            )
        end

        fields_keys[i] = Symbol(split(fields[i], ",")[1])
        unit_string = split(split(fields[i], ",")[3], "]")[1]
        fields_units[i] = uparse(unit_string)
    end

    # Determine the range 
    T0 = parse(Float64, header[50]) * u"K"      # in K
    dT = parse(Float64, header[51]) * u"K"
    numT = parse(Int64, header[52])

    P0 = parse(Float64, header[53]) * u"bar"    # in bar (will be convert to Pa later)
    dP = parse(Float64, header[54]) * u"bar"
    numP = parse(Int64, header[55])

    Tvec = T0:dT:(T0 + dT * (numT - 1))              # 1D vector
    Pvec = P0:dP:(P0 + dP * (numP - 1))

    # In the LaMEM/Perple_X file format, the first 50 lines are comments
    data = readdlm(fname; skipstart=55, header=false)        # read numerical data

    # Shape of 2D arrays:
    siz = (numT, numP)

    # Initialize fields in the order they are defined in the PhaseDiagram_LookupTable structure 
    Struct_Fieldnames = fieldnames(PhaseDiagram_LookupTable)[4:end] # fieldnames from structure

    # Process all fields that are present in the phase diagram (and non-dimensionalize if requested)
    Struct_Fields = Vector{Union{Nothing,Interpolations.Extrapolation}}(
        nothing, length(Struct_Fieldnames)
    )

    # Loop through all fields and create an interpolation object for each (which is optionally non-dimensionalized)
    for (i, field) in enumerate(fields_keys)
        ind = findall(Struct_Fieldnames .== field)
        if length(ind) > 0
            Struct_Fields[ind[1]] = CreateInterpolationObject_PhaseDiagram(
                data[:, i], Tvec, Pvec, siz, fields_units[i], CharDim
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
        "Perple_X/MAGEMin/LaMEM", header_text, fname, Struct_Fields...
    )

    return PD_data
end

# Print info 
function show(io::IO, d::PhaseDiagram_LookupTable)
    T = d.rockRho.itp.knots[1]
    P = d.rockRho.itp.knots[2]

    println(io, "$(d.Type) Phase Diagram Lookup Table: ")
    println(io, "                      File    :   $(d.Name)")
    println(io, "                      T       :   $(minimum(T)) - $(maximum(T))")
    println(io, "                      P       :   $(minimum(P)) - $(maximum(P))")

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

# Internal routine that creates an interpolation object from a collumn of the data
function CreateInterpolationObject_PhaseDiagram(
    data_vec::Vector{Float64}, Tvec, Pvec, siz::Tuple{Int64,Int64}, units, CharDim
)
    data_units = reshape(data_vec, siz) * units      # Create 2D array 

    # Convert to Pa & K
    Pvec_Pa = Float64.(uconvert.(u"Pa", Pvec))
    Tvec_K = Float64.(uconvert.(u"K", Tvec))

    # Optional: nondimensionalize values as well as T and P
    if CharDim == nothing
        Pvec = ustrip.(Pvec_Pa)
        Tvec = ustrip.(Tvec_K)
        data = ustrip.(data_units)
    else
        Pvec = nondimensionalize(Pvec_Pa, CharDim)
        Tvec = nondimensionalize(Tvec_K, CharDim)
        data = nondimensionalize(data_units, CharDim)
    end

    # Create interpolation object
    intp_data = LinearInterpolation((Tvec, Pvec), data; extrapolation_bc=Flat())

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
            ϕ = Struct_Fields[ind_meltFrac[1]].itp.coefs      # melt fraction
            S = Struct_Fields[ind_solidData[1]].itp.coefs      # solid property
            M = Struct_Fields[ind_meltData[1]].itp.coefs      # melt property

            # Compute average
            Result = (1.0 .- ϕ) .* S .+ ϕ .* M

            # Create interpolation object of average
            Tvec = Struct_Fields[ind_meltFrac[1]].itp.knots[1]       # temperature data
            Pvec = Struct_Fields[ind_meltFrac[1]].itp.knots[2]       # pressure data
            intp_Result = LinearInterpolation((Tvec, Pvec), Result; extrapolation_bc=Flat())

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

end
