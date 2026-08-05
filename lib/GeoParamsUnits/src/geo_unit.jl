"""
    AbstractGeoUnit{T, U}

Abstract supertype for `GeoUnit` values with numeric storage type `T` and unit type `U`.
"""
abstract type AbstractGeoUnit{T, U} end

abstract type AbstractUnitType end

struct GEO <: AbstractUnitType end
struct SI <: AbstractUnitType end
struct NONE <: AbstractUnitType end

"""
    GeoUnit{T, U}

Store a numeric value, its units, and whether the value is currently dimensional.
Nondimensionalization retains `unit`, allowing the value to be dimensionalized later.
The dimensional state is stored in `isdimensional` and is not encoded in the type.
"""
struct GeoUnit{T, U} <: AbstractGeoUnit{T, U}
    val::T
    unit::U
    isdimensional::Bool
end

@inline function _representative_value(val::AbstractArray)
    isempty(val) && throw(ArgumentError("cannot construct a GeoUnit from an empty array"))
    return first(val)
end

@inline _representative_value(val) = val

# Different ways of specifying the GeoUnit:
function GeoUnit(val)
    representative = _representative_value(val)
    return GeoUnit{typeof(ustrip.(val)), typeof(unit(representative))}(
        ustrip.(val),
        unit(representative),
        isa(representative, Union{Unitful.FreeUnits, Unitful.Quantity}),
    )
end

function GeoUnit(val::Nothing)
    return GeoUnit{Nothing, typeof(NoUnits)}(val, NoUnits, false)
end

function GeoUnit(fun::F) where {F <: Function}
    return GeoUnit{typeof(fun), typeof(NoUnits)}(fun, NoUnits, false)
end

function GeoUnit{T}(val) where {T}
    representative = _representative_value(val)
    return GeoUnit{T, typeof(unit(representative))}(
        T.(ustrip.(val)),
        unit(representative),
        isa(representative, Union{Unitful.FreeUnits, Unitful.Quantity}),
    )
end

function GeoUnit{T, U}(val) where {T, U}
    representative = _representative_value(val)
    return GeoUnit{T, U}(
        T.(ustrip.(val)),
        unit(representative),
        isa(representative, Union{Unitful.FreeUnits, Unitful.Quantity}),
    )
end

Base.convert(t::Type{GeoUnit{T, U}}, x::Quantity{T, U, X}) where {T, U, X} = GeoUnit(x)
Base.convert(t::Type{GeoUnit{T, U}}, x::T) where {T, U} = GeoUnit(x)
# Handle cross-type conversions (e.g. GeoUnit{Float32} -> GeoUnit{Float64})
Base.convert(::Type{GeoUnit{T, U}}, x::GeoUnit{S, U}) where {T, U, S} = GeoUnit{T, U}(T(x.val), x.unit, x.isdimensional)

# Define methods to deal with cases when the input has integers
GeoUnit(val::Union{Int64, AbstractArray{Int64}}) = GeoUnit(Float64.(val))
GeoUnit(val::Union{Int32, AbstractArray{Int32}}) = GeoUnit(Float32.(val))
function GeoUnit(val::Union{Quantity{Int64}, AbstractArray{<:Quantity{<:Int64}}})
    return GeoUnit(Float64.(val))
end
function GeoUnit(val::Union{Quantity{Int32}, AbstractArray{<:Quantity{<:Int32}}})
    return GeoUnit(Float32.(val))
end

# helper functions
Unit(v::GeoUnit{T, U}) where {T, U} = Unitful.unit(v.unit * 1)
isdimensional(v::GeoUnit{T, U}) where {T, U} = v.isdimensional                 # is it a nondimensional number or not?
isdimensional(v::Number) = false                            # nope
isDimensional(v::GeoUnit) = isdimensional(v)
isDimensional(v::Number) = false
NumValue(v::GeoUnit) = v.val                                # numeric value, with no units
NumValue(v::Number) = v                                     # numeric value
NumValue(v::AbstractArray) = v                              # numeric value
Value(v::GeoUnit) = Unitful.Quantity.(v.val, v.unit)        # value, with units
Fun(v::GeoUnit) = v.val
unpack_units(x::NTuple{N, GeoUnit}) where {N} = ntuple(i -> x[i].unit * x[i].val, Val(N))
unpack_vals(x::NTuple{N, GeoUnit}) where {N} = ntuple(i -> x[i].val, Val(N))

function UnitValue(v::GeoUnit{T, U}) where {T, U}
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
Base.convert(::Type{<:Real}, v::GeoUnit) = v.val
Base.convert(::Type{GeoUnit}, v::Number) = GeoUnit(v)
Base.convert(::Type{GeoUnit}, v::Int32) = GeoUnit(Float32(v))
Base.convert(::Type{GeoUnit}, v::Int64) = GeoUnit(Float64(v))
Base.convert(::Type{GeoUnit}, v::Quantity) = GeoUnit(v)
Base.convert(::Type{GeoUnit}, v::AbstractArray) = GeoUnit(v)
Base.convert(::Type{GeoUnit}, v::AbstractArray{Int32}) = GeoUnit(v)
Base.convert(::Type{GeoUnit}, v::AbstractArray{Int64}) = GeoUnit(v)

Base.convert(::Type{GeoUnit{T}}, v::Quantity) where {T} = GeoUnit(T.(v))
Base.convert(::Type{GeoUnit{T}}, v::Number) where {T} = GeoUnit(T(v))
Base.convert(::Type{GeoUnit{T}}, v::AbstractArray) where {T} = GeoUnit(T.(v))

#Base.convert(::Type{GeoUnit{T,U}},  v::T)    where {T,U}    =   GeoUnit{T,typeof(unit(v[1]))}(v)

Base.promote_rule(::Type{GeoUnit}, ::Type{Quantity}) = GeoUnit

function Base.show(io::IO, x::GeoUnit{T, U}) where {T, U} # output
    val = x.val
    if x.isdimensional == true
        println("GeoUnit{dimensional, $(x.unit)}, ")
    else
        println("GeoUnit{nondimensional, $(x.unit)}, ")
    end
    return show(io, MIME("text/plain"), val)
end
