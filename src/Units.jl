"""
    This provides units and creates a non-dimensionalization object
"""
module Units
using Unitful
import Unitful: superscript
using Parameters
using Setfield # allows modifying fields in immutable struct

import Base:
    show, isapprox, isequal, convert, length, size, getindex, setindex!, getproperty
import Base.Broadcast: broadcasted
using GeoParams:
    AbstractMaterialParam,
    AbstractMaterialParamsStruct,
    AbstractPhaseDiagramsStruct,
    PerpleX_LaMEM_Diagram

# Define additional units that are useful in geodynamics 
@unit Myrs "Myrs" MillionYears 1000000u"yr" false

function __init__()
    return Unitful.register(Units)
end

Unitful.register(Units)

# Define a number of useful units
const km = u"km"
const m = u"m"
const cm = u"cm"
const mm = u"mm"
const Myrs = u"Myrs"
const yr = u"yr"
const s = u"s"
const kg = u"kg"
const g = u"g"
const Pa = u"Pa"
const MPa = u"MPa"
const kbar = u"kbar"
const Pas = u"Pa*s"
const K = u"K"
const C = u"°C"
const mol = u"mol"
const kJ = u"kJ"
const J = u"J"
const Watt = u"W"
const μW = u"μW"
const μm = u"μm"

export km,
    m,
    cm,
    mm,
    μm,
    Myrs,
    yr,
    s,
    MPa,
    Pa,
    kbar,
    Pas,
    K,
    C,
    g,
    kg,
    mol,
    J,
    kJ,
    Watt,
    μW,
    GeoUnit,
    GeoUnits,
    GEO_units,
    SI_units,
    NO_units,
    AbstractGeoUnit,
    nondimensionalize,
    dimensionalize,
    superscript,
    upreferred,
    GEO,
    SI,
    NONE,
    isDimensional,
    Value,
    NumValue,
    Unit,
    UnitValue,
    isdimensional,
    compute_units

include("unpack.jl")    # adds macros for unpacking GeoUnit variables with or w/out units

"""
AbstractGeoUnit

Abstract supertype for geo units.
"""
abstract type AbstractGeoUnit{TYPE,DIMENSIONAL} <: Number end

abstract type AbstractUnitType end

struct GEO <: AbstractUnitType end
struct SI <: AbstractUnitType end
struct NONE <: AbstractUnitType end

# The GeoUnit struct encodes dimensional info in the type info
"""
    Structure that holds a GeoUnit parameter, and encodes the units and whether it is dimensional or not in the name.

    Having that is useful, as non-dimensionalization removes the units from a number and we thus no longer know how to transfer it back to the correct units.
    With the GeoUnit struct, this information is retained, and we can thus redimensionalize it at a later StepRange


"""
struct GeoUnit{T,U}
    val::T
    unit::U
    isdimensional::Bool
end

# Different ways of specifying the GeoUnit:
function GeoUnit(val)
    return GeoUnit{typeof(ustrip.(val)),typeof(unit(val[1]))}(
        ustrip.(val), unit(val[1]), isa(val[1], Union{Unitful.FreeUnits,Unitful.Quantity})
    )
end

function GeoUnit{T}(val) where {T,U}
    return GeoUnit{T,typeof(unit(val[1]))}(
        T.(ustrip.(val)),
        unit(val[1]),
        isa(val[1], Union{Unitful.FreeUnits,Unitful.Quantity}),
    )
end

function GeoUnit{T,U}(val) where {T,U}
    return GeoUnit{T,U}(
        T.(ustrip.(val)),
        unit(val[1]),
        isa(val[1], Union{Unitful.FreeUnits,Unitful.Quantity}),
    )
end

Base.convert(t::Type{GeoUnit{T,U}}, x::Quantity{T,U,X}) where {T,U,X} = GeoUnit(x)
Base.convert(t::Type{GeoUnit{T,U}}, x::T) where {T,U} = GeoUnit(x)

# Define methods to deal with cases when the input has integers 
GeoUnit(val::Union{Int64,AbstractArray{Int64}}) = GeoUnit(Float64.(val))
GeoUnit(val::Union{Int32,AbstractArray{Int32}}) = GeoUnit(Float32.(val))
function GeoUnit(val::Union{Quantity{Int64},AbstractArray{<:Quantity{<:Int64}}})
    return GeoUnit(Float64.(val))
end
function GeoUnit(val::Union{Quantity{Int32},AbstractArray{<:Quantity{<:Int32}}})
    return GeoUnit(Float32.(val))
end

# helper functions
Unit(v::GeoUnit{T,U}) where {T,U} = Unitful.unit(v.unit * 1)
isdimensional(v::GeoUnit{T,U}) where {T,U} = v.isdimensional                 # is it a nondimensional number or not?
isdimensional(v::Number) = false                           # nope
NumValue(v::GeoUnit) = v.val                           # numeric value, with no units
NumValue(v::Number) = v                               # numeric value
NumValue(v::AbstractArray) = v                               # numeric value
Value(v::GeoUnit) = Unitful.Quantity.(v.val, v.unit)  # value, with units

function UnitValue(v::GeoUnit{T,U}) where {T,U}
    if v.isdimensional
        return Value(v)             # returns value with units
    else
        return NumValue(v.val)
    end
end

Base.isapprox(x::GeoUnit, y::Number; kwargs...) = Base.isapprox(x.val, y; kwargs...)
Base.isequal(x::GeoUnit, y::Number) = Base.isequal(x.val, y)
Base.isequal(x::GeoUnit, y::AbstractArray) = Base.isequal(x.val, y)
Base.isequal(x::GeoUnit, y::GeoUnit) = Base.isequal(x.val, y.val)

Base.convert(::Type{<:AbstractArray}, v::GeoUnit) = v.val
Base.convert(::Type{<:Number}, v::GeoUnit) = v.val
Base.convert(::Type{GeoUnit}, v::Number) = GeoUnit(v)
Base.convert(::Type{GeoUnit}, v::Int32) = GeoUnit(Float32(v))
Base.convert(::Type{GeoUnit}, v::Int64) = GeoUnit(Float64(v))
Base.convert(::Type{GeoUnit}, v::Quantity) = GeoUnit(v)
Base.convert(::Type{GeoUnit}, v::AbstractArray) = GeoUnit(v)
function Base.convert(::Type{GeoUnit}, v::AbstractArray{Int32})
    return GeoUnit(ustrip.(Float32.(v)) * unit(v[1]))
end
function Base.convert(::Type{GeoUnit}, v::AbstractArray{Int64})
    return GeoUnit(ustrip.(Float64.(v)) * unit(v[1]))
end

Base.convert(::Type{GeoUnit{T}}, v::Quantity) where {T} = GeoUnit(T.(v))
Base.convert(::Type{GeoUnit{T}}, v::Number) where {T} = GeoUnit(T(v))
Base.convert(::Type{GeoUnit{T}}, v::AbstractArray) where {T} = GeoUnit(T.(v))

#Base.convert(::Type{GeoUnit{T,U}},  v::T)    where {T,U}    =   GeoUnit{T,typeof(unit(v[1]))}(v) 

Base.promote_rule(::Type{GeoUnit}, ::Type{Quantity}) = GeoUnit

function Base.show(io::IO, x::GeoUnit{T,U}) where {T,U} # output
    val = x.val
    if x.isdimensional == true
        println("GeoUnit{dimensional, $(x.unit)}, ")
    else
        println("GeoUnit{nondimensional, $(x.unit)}, ")
    end
    return show(io, MIME("text/plain"), val)
end

# define a few basic routines so we can easily operate with GeoUnits
Base.length(v::GeoUnit) = length(v.val)
Base.size(v::GeoUnit) = size(v.val)
Base.getindex(A::GeoUnit{T,U}, inds::Vararg{Int,N}) where {T,U,N} = A.val[inds...]

for op in (:+, :-, :*, :/)

    # Multiply with number
    @eval Base.$op(x::GeoUnit, y::Number) = $(op)(x.val, y)
    @eval Base.$op(x::Number, y::GeoUnit) = $(op)(x, y.val)

    # Multiplying a GeoUnit with another one, returns a GeoUnit
    @eval function Base.$op(x::GeoUnit{T1,U1}, y::GeoUnit{T2,U2}) where {T1,T2,U1,U2}
        return GeoUnit($(op)(UnitValue(x), UnitValue(y)))
    end
    @eval Base.$op(x::GeoUnit, y::Quantity) = $(op).(UnitValue(x), y)
    @eval Base.$op(x::Quantity, y::GeoUnit) = $(op).(x, UnitValue(y))

    # If we multiply a GeoUnit with an abstract array, we only return values, not units (use GeoUnits for that)
    @eval Base.$op(x::GeoUnit, y::AbstractArray) = broadcast($op, NumValue(x), y)
    @eval Base.$op(x::AbstractArray, y::GeoUnit) = broadcast($op, x, NumValue(y))
    @eval Base.$op(x::GeoUnit, y::AbstractArray{<:Quantity}) = broadcast($op, Value(x), y)
    @eval Base.$op(x::AbstractArray{<:Quantity}, y::GeoUnit) = broadcast($op, x, Value(y))
    # Broadcasting
    @eval function Base.broadcasted(::typeof($(op)), A::GeoUnit, B::AbstractArray)
        return broadcast($(op), NumValue(A), B)
    end
    @eval function Base.broadcasted(::typeof($(op)), A::AbstractArray, B::GeoUnit)
        return broadcast($(op), A, NumValue(B))
    end
    @eval function Base.broadcasted(
        ::typeof($(op)), A::GeoUnit, B::AbstractArray{<:Quantity}
    )
        return broadcast($(op), Value(A), B)
    end
    @eval function Base.broadcasted(
        ::typeof($(op)), A::AbstractArray{<:Quantity}, B::GeoUnit
    )
        return broadcast($(op), A, Value(B))
    end
end

# Broadcasting
Base.broadcasted(::typeof(^), A::GeoUnit, B::AbstractArray) = broadcast(^, NumValue(A), B)
Base.broadcasted(::typeof(^), A::AbstractArray, B::GeoUnit) = broadcast(^, A, NumValue(B))
function Base.broadcasted(::typeof(^), A::GeoUnit, B::AbstractArray{<:Quantity})
    return broadcast(^, Value(A), B)
end
function Base.broadcasted(::typeof(^), A::AbstractArray{<:Quantity}, B::GeoUnit)
    return broadcast(^, A, Value(B))
end

Base.getindex(x::GeoUnit, i::Int64, j::Int64, k::Int64) = GeoUnit(x.val[i, j, k] * x.unit)
Base.getindex(x::GeoUnit, i::Int64, j::Int64) = GeoUnit(x.val[i, j] * x.unit)
Base.getindex(x::GeoUnit, i::Int64) = GeoUnit(x.val[i] * x.unit)
function Base.setindex!(A::GeoUnit{T,U}, val, inds::Vararg{Int,N}) where {T,U,N}
    return A.val[inds...] = val
end

"""
    GeoUnits

    Structure that holds parameters used for non-dimensionalization
"""
@with_kw_noshow struct GeoUnits{TYPE}
    # Selectable input parameters
    temperature = 1               #   Characteristic temperature  [C or K]
    length = 1               #   Characteristic length unit  [km, m or -]
    stress = 1               #   Characteristic stress unit  [MPa, Pa or -]
    time = 1               #   Characteristic time unit    [Myrs, s pr -]
    viscosity = 1

    # primary characteristic units
    K = 1                      # temperature in SI units for material parameter scaling
    s = 1                      # time in SI units for material parameter scaling
    m = 1                      # length in SI units for material parameter scaling
    Pa = 1                      # stress in SI units 
    kg = upreferred(Pa * m * s^2)    # compute mass from pascal. Note: this may result in very large values

    # Main SI units (used later for nondimensionalization)
    Length = m
    Mass = kg
    Time = s
    Temperature = K
    Amount = 1mol
    Second = s

    # Not defined, as they are not common in geodynamics:
    #Current
    #Luminosity

    # Derived units
    N = kg * m / s^2        # Newton
    J = N * m             # Joule
    W = J / s             # Watt
    area = m^2             # area   in SI units for material parameter scaling
    volume = m^3             # volume
    velocity = m / s
    density = kg / m^3
    acceleration = m / s^2
    force = kg * m / s^2
    strainrate = 1 / s
    heatcapacity = J / kg / K
    conductivity = W / m / K

    # Helpful
    SecYear = 3600 * 24 * 365.25
    Myrs = 1e6
    cmYear = SecYear * 100     # to transfer m/s -> cm/yr
end

"""
    GEO_units(;length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)

Creates a non-dimensionalization object using GEO units.

GEO units implies that upon dimensionalization, `time` will be in `Myrs`, `length` in `km`, stress in `MPa`, etc.
which is more convenient for typical geodynamic simulations than SI units
The characteristic values given as input can be in arbitrary units (`km` or `m`), provided the unit is specified.

# Examples:
```julia-repl
julia> CharUnits = GEO_units()
Employing GEO units 
Characteristic values: 
         length:      1000 km
         time:        0.3169 Myrs
         stress:      10 MPa
         temperature: 1000.0 °C
julia> CharUnits.velocity
1.0e-7 m s⁻¹
```
If we instead have a crustal-scale simulation, it is likely more appropriate to use a different characteristic `length`:
```julia-repl
julia> CharUnits = GEO_units(length=10km)
Employing GEO units 
Characteristic values: 
         length:      10 km
         time:        0.3169 Myrs
         stress:      10 MPa
         temperature: 1000.0 °C
```
"""
function GEO_units(; length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
    if unit(temperature) == NoUnits
        temperature = temperature * C
    end
    if unit(length) == NoUnits
        length = length * u"km"
    end
    if unit(stress) == NoUnits
        stress = stress * u"MPa"
    end
    if unit(viscosity) == NoUnits
        viscosity = viscosity * u"Pa*s"
    end

    T = uconvert(C, temperature)
    Le = uconvert(km, length)
    Sigma = uconvert(MPa, stress)
    Eta = uconvert(Pas, viscosity)

    T_SI = uconvert(K, T)
    Le_SI = uconvert(m, Le)
    Sigma_SI = uconvert(Pa, Sigma)
    Time_SI = Eta / Sigma_SI
    t = uconvert(Myrs, Time_SI)

    return GeoUnits{GEO}(;
        length=Le,
        temperature=T,
        stress=Sigma,
        viscosity=Eta,
        time=t,
        m=Le_SI,
        K=T_SI,
        Pa=Sigma_SI,
        s=Time_SI,
    )
end

"""
    CharUnits = SI_units(length=1000m, temperature=1000K, stress=10Pa, viscosity=1e20)

Specify the characteristic values using SI units 

# Examples:
```julia-repl
julia> CharUnits = SI_units(length=1000m)
Employing SI units 
Characteristic values: 
         length:      1000 m
         time:        1.0e19 s
         stress:      10 Pa
         temperature: 1000.0 K
```
Note that the same can be achieved if the input is given in `km`:
```julia-repl
julia> CharUnits = SI_units(length=1km)
```
"""
function SI_units(; length=1000m, temperature=1000K, stress=10Pa, viscosity=1e20Pas)
    if unit(temperature) == NoUnits
        temperature = temperature * K
    end
    if unit(length) == NoUnits
        length = length * u"m"
    end
    if unit(stress) == NoUnits
        stress = stress * u"Pa"
    end
    if unit(viscosity) == NoUnits
        viscosity = viscosity * u"Pa*s"
    end

    T = uconvert(K, temperature)
    Le = uconvert(m, length)
    Sigma = uconvert(Pa, stress)
    Eta = uconvert(Pas, viscosity)

    T_SI = uconvert(K, T)
    Le_SI = uconvert(m, Le)
    Sigma_SI = uconvert(Pa, Sigma)
    Time_SI = Eta / Sigma_SI
    t = uconvert(s, Time_SI)

    return GeoUnits{SI}(;
        length=Le,
        temperature=T,
        stress=Sigma,
        viscosity=Eta,
        time=t,
        m=Le_SI,
        K=T_SI,
        Pa=Sigma_SI,
        s=Time_SI,
    )
end

"""
    CharUnits = NO_units(length=1, temperature=1, stress=1, viscosity=1)
   
Specify the characteristic values in non-dimensional units

# Examples:
```julia-repl
julia> using GeoParams;
julia> CharUnits = NO_units()
Employing NONE units 
Characteristic values: 
         length:      1
         time:        1.0 
         stress:      1
         temperature: 1.0
```
"""
function NO_units(; length=1, temperature=1, stress=1, viscosity=1)
    if unit(temperature) != NoUnits
        error("temperature should not have units")
    end
    if unit(length) != NoUnits
        error("length should not have units")
    end
    if unit(stress) != NoUnits
        error("stress should not have units")
    end
    if unit(viscosity) != NoUnits
        error("viscosity should not have units")
    end

    T = temperature
    Le = length
    Sigma = stress
    Eta = viscosity
    Time = Eta / Sigma

    return GeoUnits{NONE}(;
        length=Le,
        temperature=T,
        stress=Sigma,
        viscosity=Eta,
        time=Time,
        m=Le,
        K=T,
        Pa=Sigma,
        s=Time,
    )
end

"""
    nondimensionalize(param, CharUnits::GeoUnits{TYPE})

Nondimensionalizes `param` using the characteristic values specified in `CharUnits`

# Example 1
```julia-repl
julia> using GeoParams;
julia> CharUnits =   GEO_units();
julia> v         =   3cm/yr
3 cm yr⁻¹ 
julia> v_ND      =   nondimensionalize(v, CharUnits) 
0.009506426344208684
```
# Example 2
In geodynamics one sometimes encounters more funky units
```julia-repl
julia> CharUnits =   GEO_units();
julia> A         =   6.3e-2MPa^-3.05*s^-1
0.063 MPa⁻³·⁰⁵ s⁻¹
julia> A_ND      =   nondimensionalize(A, CharUnits) 
7.068716262102384e14
```

In case you are interested to see how the units of `A` look like in different units, use this function from the [Unitful](https://github.com/PainterQubits/Unitful.jl) package:
```julia-repl
julia> uconvert(u"Pa^-3.05*s^-1",A) 
3.157479571851836e-20 Pa⁻³·⁰⁵
```
and to see it decomposed in the basic `SI` units of length, mass and time:
```julia-repl
julia> upreferred(A)
3.1574795718518295e-20 m³·⁰⁵ s⁵·¹ kg⁻³·⁰⁵
```
"""
function nondimensionalize(param::GeoUnit{T,U}, g::GeoUnits{TYPE}) where {T,U,TYPE}
    if param.isdimensional
        char_val = compute_units(param, g)
        #    @show char_val, typeof(char_val)
        val_ND = upreferred.(param.val * param.unit) / char_val
        val_ND_number = convert(T, val_ND)
        param_ND = GeoUnit{T,U}(val_ND_number, param.unit, false)          # store new value, but keep original dimensions

    else
        param_ND = param
    end
    return param_ND
end

# in case the parameter is already non-dimensional:
nondimensionalize(param::String, g::GeoUnits{TYPE}) where {TYPE} = param

# in case it is a unitful quantity
function nondimensionalize(
    param::Union{Unitful.Quantity{T,K,M},AbstractArray{<:Quantity{T,K,M}}},
    g::GeoUnits{TYPE},
) where {TYPE,T,K,M}
    param_Geo = GeoUnit(param)
    result = nondimensionalize(param_Geo, g)

    return UnitValue(result)
end

# If it is an array, but has no units we cannot know how to nondimensionalize it
nondimensionalize(param::AbstractArray{<:Number}, g::GeoUnits{TYPE}) where {TYPE} = param

# This computes the characteristic value
function compute_units(
    param::GeoUnit{<:Union{T,AbstractArray{T}},U}, g::GeoUnits{TYPE}
) where {T,U,TYPE}
    dim = Unitful.dimension(param.unit)              # Basic SI units

    # Note: `char_val`` is type unstable within the routine
    # Yet, the dimensions at the end are known as they must be the same as `U` 
    #  (perhaps converted to preferred units such as m instead of km)
    char_val = 1.0
    foreach((typeof(dim).parameters[1])) do y
        val = upreferred(getproperty(g, Unitful.name(y)))       # Retrieve the characteristic value from structure g
        pow = Float64(y.power)                                  # power by which it should be multiplied   
        char_val *= val^pow                                     # multiply characteristic value
    end

    value::T = ustrip(char_val)                                 # numerical value
    char_val_out = upreferred(param.unit) * value                # this is done for type-stability

    return char_val_out
end

"""
    MatParam_ND = nondimensionalize(MatParam::AbstractMaterialParam, CharUnits::GeoUnits{TYPE})

Non-dimensionalizes a material parameter structure (e.g., Density, CreepLaw)

"""
function nondimensionalize(MatParam::AbstractMaterialParam, g::GeoUnits{TYPE}) where {TYPE}
    for param in fieldnames(typeof(MatParam))
        if isa(getfield(MatParam, param), GeoUnit)
            z = getfield(MatParam, param)

            # non-dimensionalize:
            z = nondimensionalize(z, g)

            # Replace field (using Setfield package):
            MatParam = set(MatParam, Setfield.PropertyLens{param}(), z)
        elseif isa(getfield(MatParam, param), AbstractMaterialParam)
            # The field contains another AbstractMaterialParam
            z = getfield(MatParam, param)
            z = nondimensionalize(z, g)
            MatParam = set(MatParam, Setfield.PropertyLens{param}(), z)
        end
    end
    return MatParam
end

"""
    nondimensionalize(phase_mat::MaterialParams, g::GeoUnits{TYPE})

nondimensionalizes all fields within the Material Parameters structure that contain material parameters
"""
function nondimensionalize(
    phase_mat::AbstractMaterialParamsStruct, g::GeoUnits{TYPE}
) where {TYPE}
    for param in fieldnames(typeof(phase_mat))
        fld = getfield(phase_mat, param)
        if length(fld) > 0
            if typeof(fld[1]) <: AbstractPhaseDiagramsStruct

                # in case we employ a phase diagram 
                temp = PerpleX_LaMEM_Diagram(fld[1].Name; CharDim=g)
                fld_new = (temp,)
                #setfield!(phase_mat, param, fld_new)
                phase_mat = set(phase_mat, Setfield.PropertyLens{param}(), fld_new)

            else
                # otherwise non-dimensionalize all fields and create a new tuple
                id = findall(isa.(fld, AbstractMaterialParam))
                if length(id) > 0
                    # Create a new tuple with non-dimensionalized fields:
                    fld_new = ntuple(
                        i -> Units.nondimensionalize(fld[id[i]], g), length(id)
                    )
                    #                    setfield!(phase_mat, param, fld_new)        # to be changed for immutable struct
                    phase_mat = set(phase_mat, Setfield.PropertyLens{param}(), fld_new)
                end
            end
        end
    end
    # phase_mat.Nondimensional = true
    phase_mat = set(phase_mat, Setfield.PropertyLens{:Nondimensional}(), true)
    return phase_mat
end

"""
    dimensionalize(param, param_dim::Unitful.FreeUnits, CharUnits::GeoUnits{TYPE})

Dimensionalizes `param` into the dimensions `param_dim` using the characteristic values specified in `CharUnits`.  

# Example
```julia-repl
julia> CharUnits =   GEO_units();
julia> v_ND      =   nondimensionalize(3cm/yr, CharUnits) 
0.031688087814028945
julia> v_dim     =   dimensionalize(v_ND, cm/yr, CharUnits) 
3.0 cm yr⁻¹
```

"""
function dimensionalize(
    param_ND, param_dim::Unitful.FreeUnits, g::GeoUnits{TYPE}
) where {TYPE}
    char_val = compute_units(GeoUnit(1.0 * param_dim), g)         # Determine characteristic units
    param = uconvert.(param_dim, param_ND * char_val)

    return param
end

function dimensionalize(param_ND::GeoUnit{T,U}, g::GeoUnits{TYPE}) where {T,U,TYPE}
    if isdimensional(param_ND) == false
        char_val = compute_units(param_ND, g)                       # Determine characteristic units
        val = uconvert.(param_ND.unit, param_ND.val * char_val)   # dimensionalize
        param = GeoUnit{T,U}(ustrip.(val), param_ND.unit, true)   # store new value, but keep original dimensions

        return param
    else
        return param_ND
    end
end

"""
    dimensionalize(MatParam::AbstractMaterialParam, CharUnits::GeoUnits{TYPE})

Dimensionalizes a material parameter structure (e.g., Density, CreepLaw)

"""
function dimensionalize(MatParam::AbstractMaterialParam, g::GeoUnits{TYPE}) where {TYPE}
    for param in fieldnames(typeof(MatParam))
        if isa(getfield(MatParam, param), GeoUnit)
            z = getfield(MatParam, param)
            z_dim = dimensionalize(z, g)

            # Replace field (using Setfield package):
            MatParam = set(MatParam, Setfield.PropertyLens{param}(), z_dim)

        elseif isa(getfield(MatParam, param), AbstractMaterialParam)
            # The field contains another AbstractMaterialParam
            z = getfield(MatParam, param)
            z_dim = dimensionalize(z, g)
            MatParam = set(MatParam, Setfield.PropertyLens{param}(), z_dim)
        end
    end
    return MatParam
end

"""
    Dimensionalize(phase_mat::MaterialParams, g::GeoUnits{TYPE})

Dimensionalizes all fields within the Material Parameters structure that contain material parameters
"""
function dimensionalize(
    phase_mat::AbstractMaterialParamsStruct, g::GeoUnits{TYPE}
) where {TYPE}
    for param in fieldnames(typeof(phase_mat))
        fld = getfield(phase_mat, param)
        if length(fld) > 0
            id = findall(isa.(fld, AbstractMaterialParam))
            if length(id) > 0
                # Create a new tuple with non-dimensionalized fields:
                fld_new = ntuple(i -> Units.dimensionalize(fld[id[i]], g), length(id))

                phase_mat = set(phase_mat, Setfield.PropertyLens{param}(), fld_new)
            end
        end
    end
    phase_mat = set(phase_mat, Setfield.PropertyLens{:Nondimensional}(), false)

    return phase_mat
end

"""
    isDimensional(MatParam::AbstractMaterialParam)

`true` if MatParam is in dimensional units.    
"""
function isDimensional(MatParam::AbstractMaterialParam)
    isDim = false
    for param in fieldnames(typeof(MatParam))
        if isa(getfield(MatParam, param), GeoUnit)
            z = getfield(MatParam, param)
            if z.isdimensional == true
                isDim = true
            end
        end
    end
    return isDim
end

# Define a view for the GEO_Units structure
function show(io::IO, g::GeoUnits{TYPE}) where {TYPE}
    return print(
        io,
        "Employing $TYPE units \n",
        "Characteristic values: \n",
        "         length:      $(g.length)\n",
        "         time:        $(round(ustrip(g.time),digits=4)) $(unit(g.time))\n",
        "         stress:      $(g.stress)\n",
        "         temperature: $(Float64(g.temperature))\n",
    )
end

# This replaces the viewer of the Unitful package, such that the printing of units is done as floats (better)
function Unitful.superscript(i::Rational{Int64}; io=nothing)
    string(superscript(float(i)))
    if io === nothing
        iocontext_value = nothing
    else
        iocontext_value = get(io, :fancy_exponent, nothing)
    end
    if iocontext_value isa Bool
        fancy_exponent = iocontext_value
    else
        v = get(ENV, "UNITFUL_FANCY_EXPONENTS", Sys.isapple() ? "true" : "false")
        t = tryparse(Bool, lowercase(v))
        fancy_exponent = (t === nothing) ? false : t
    end
    if fancy_exponent
        return superscript(float(i))
    else
        return "^" * string(float(i))
    end
end

function Unitful.superscript(i::Float64)
    map(repr(i)) do c
        if c == '-'
            '\u207b'
        elseif c == '1'
            '\u00b9'
        elseif c == '2'
            '\u00b2'
        elseif c == '3'
            '\u00b3'
        elseif c == '4'
            '\u2074'
        elseif c == '5'
            '\u2075'
        elseif c == '6'
            '\u2076'
        elseif c == '7'
            '\u2077'
        elseif c == '8'
            '\u2078'
        elseif c == '9'
            '\u2079'
        elseif c == '0'
            '\u2070'
        elseif c == '.'
            '\u0387'
        else
            error("unexpected character")
        end
    end
end

end
