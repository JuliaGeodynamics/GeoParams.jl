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
"""
struct GeoUnit{T, U} <: AbstractGeoUnit{T, U}
    val::T
    unit::U
    isdimensional::Bool
end

@inline function _representative_value(value::AbstractArray)
    isempty(value) && throw(ArgumentError("cannot construct a GeoUnit from an empty array"))
    return first(value)
end

@inline _representative_value(value) = value

@inline _promote_storage(value::Int32) = Float32(value)
@inline _promote_storage(value::Int64) = Float64(value)
@inline _promote_storage(value::AbstractArray{Int32}) = Float32.(value)
@inline _promote_storage(value::AbstractArray{Int64}) = Float64.(value)
@inline _promote_storage(value) = value

@inline _storage_type(::Type{Int32}) = Float32
@inline _storage_type(::Type{Int64}) = Float64
@inline _storage_type(::Type{T}) where {T} = T

function GeoUnit(value::DynamicQuantities.AbstractQuantity{T}) where {T}
    S = _storage_type(T)
    stored_value = convert(S, ustrip(value))
    stored_unit = unit(value)
    return GeoUnit{S, typeof(stored_unit)}(
        stored_value,
        stored_unit,
        true,
    )
end

function GeoUnit(value::T) where {T <: Number}
    S = _storage_type(T)
    return GeoUnit{S, typeof(NoUnits)}(convert(S, value), NoUnits, false)
end

function GeoUnit(
        value::AbstractArray{<:DynamicQuantities.AbstractQuantity{T}},
    ) where {T}
    representative = _representative_value(value)
    S = _storage_type(T)
    stored_value = S.(ustrip.(value))
    stored_unit = unit(representative)
    return GeoUnit{typeof(stored_value), typeof(stored_unit)}(
        stored_value,
        stored_unit,
        true,
    )
end

function GeoUnit(value::AbstractArray{T}) where {T <: Number}
    _representative_value(value)
    S = _storage_type(T)
    stored_value = S.(value)
    return GeoUnit{typeof(stored_value), typeof(NoUnits)}(stored_value, NoUnits, false)
end

GeoUnit(::Nothing) = GeoUnit{Nothing, typeof(NoUnits)}(nothing, NoUnits, false)
GeoUnit(fun::F) where {F <: Function} = GeoUnit{F, typeof(NoUnits)}(fun, NoUnits, false)

function GeoUnit{T}(value) where {T}
    representative = _representative_value(value)
    stored_unit = unit(representative)
    return GeoUnit{T, typeof(stored_unit)}(
        T.(ustrip.(value)),
        stored_unit,
        _isquantity(representative),
    )
end

function GeoUnit{T, U}(value) where {T, U}
    representative = _representative_value(value)
    return GeoUnit{T, U}(
        T.(ustrip.(value)),
        unit(representative),
        _isquantity(representative),
    )
end

Base.convert(::Type{GeoUnit{T, U}}, value::GeoUnit{S, U}) where {T, U, S} =
    GeoUnit{T, U}(T(value.val), value.unit, value.isdimensional)

Unit(value::GeoUnit) = value.unit
isdimensional(value::GeoUnit) = value.isdimensional
isdimensional(::Number) = false
isDimensional(value::GeoUnit) = isdimensional(value)
isDimensional(::Number) = false
NumValue(value::GeoUnit) = value.val
NumValue(value::Number) = value
NumValue(value::AbstractArray) = value
Value(value::GeoUnit) = value.val .* value.unit
Fun(value::GeoUnit) = value.val
unpack_units(values::NTuple{N, GeoUnit}) where {N} =
    ntuple(i -> values[i].unit .* values[i].val, Val(N))
unpack_vals(values::NTuple{N, GeoUnit}) where {N} = ntuple(i -> values[i].val, Val(N))

UnitValue(value::GeoUnit) = value.isdimensional ? Value(value) : NumValue(value)

Base.isapprox(x::GeoUnit, y::Number; kwargs...) = Base.isapprox(x.val, y; kwargs...)
Base.isequal(x::GeoUnit, y::Number) = Base.isequal(x.val, y)
Base.isequal(x::GeoUnit, y::AbstractArray) = Base.isequal(x.val, y)
Base.isequal(x::GeoUnit, y::GeoUnit) = Base.isequal(x.val, y.val)

Base.convert(::Type{<:AbstractArray}, value::GeoUnit) = value.val
Base.convert(::Type{<:Real}, value::GeoUnit) = value.val
Base.convert(::Type{GeoUnit}, value::Number) = GeoUnit(value)
Base.convert(::Type{GeoUnit}, value::DynamicQuantity) = GeoUnit(value)
Base.convert(::Type{GeoUnit}, value::AbstractArray) = GeoUnit(value)
Base.convert(::Type{GeoUnit{T}}, value::DynamicQuantity) where {T} = GeoUnit(T(ustrip(value)) * unit(value))
Base.convert(::Type{GeoUnit{T}}, value::Number) where {T} = GeoUnit(T(value))
Base.convert(::Type{GeoUnit{T}}, value::AbstractArray) where {T} = GeoUnit(T.(value))

Base.promote_rule(::Type{GeoUnit}, ::Type{<:DynamicQuantity}) = GeoUnit

function Base.show(io::IO, value::GeoUnit)
    state = value.isdimensional ? "dimensional" : "nondimensional"
    println(io, "GeoUnit{$state, $(value.unit)}, ")
    return show(io, MIME("text/plain"), value.val)
end
