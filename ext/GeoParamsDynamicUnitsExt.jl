module GeoParamsDynamicUnitsExt

import GeoParams
import GeoParamsDynamicUnits
import GeoParamsUnits
import Unitful

const GPU = GeoParamsUnits
const GPD = GeoParamsDynamicUnits
const DQ = GeoParamsDynamicUnits.DynamicQuantities

import Base: convert
import GeoParams: compute_units, dimensionalize, dimensionalize_and_strip, nondimensionalize

@inline _rational_exponent(value) = Rational(value)

function _unitful_unit(dimensions::DQ.AbstractDimensions)
    return Unitful.NoUnits *
        Unitful.u"m"^_rational_exponent(dimensions.length) *
        Unitful.u"kg"^_rational_exponent(dimensions.mass) *
        Unitful.u"s"^_rational_exponent(dimensions.time) *
        Unitful.u"A"^_rational_exponent(dimensions.current) *
        Unitful.u"K"^_rational_exponent(dimensions.temperature) *
        Unitful.u"cd"^_rational_exponent(dimensions.luminosity) *
        Unitful.u"mol"^_rational_exponent(dimensions.amount)
end

@inline _to_unitful(value::Number) = value

function _to_unitful(value::DQ.AbstractQuantity)
    expanded = DQ.uexpand(value)
    return DQ.ustrip(expanded) * _unitful_unit(DQ.dimension(expanded))
end

@inline _to_unitful(value::AbstractArray{<:DQ.AbstractQuantity}) = _to_unitful.(value)

@inline function _dimension_power(::Val{Name}, ::Unitful.Dimensions{D}) where {Name, D}
    for component in D
        Unitful.name(component) === Name && return Unitful.power(component)
    end
    return 0
end

function _dynamic_dimensions(dimensions::Unitful.Dimensions)
    return DQ.Dimensions(;
        length = _dimension_power(Val(:Length), dimensions),
        mass = _dimension_power(Val(:Mass), dimensions),
        time = _dimension_power(Val(:Time), dimensions),
        current = _dimension_power(Val(:Current), dimensions),
        temperature = _dimension_power(Val(:Temperature), dimensions),
        luminosity = _dimension_power(Val(:Luminosity), dimensions),
        amount = _dimension_power(Val(:Amount), dimensions),
    )
end

@inline _to_dynamic(value::Number) = value

function _to_dynamic(value::Unitful.AbstractQuantity)
    preferred = Unitful.upreferred(value)
    dimensions = _dynamic_dimensions(Unitful.dimension(preferred))
    return DQ.Quantity(Unitful.ustrip(preferred), dimensions)
end

@inline _to_dynamic(value::AbstractArray{<:Unitful.AbstractQuantity}) = _to_dynamic.(value)

function convert(::Type{GPU.GeoUnit}, value::GPD.GeoUnit)
    value.isdimensional && return GPU.GeoUnit(_to_unitful(GPD.Value(value)))
    unit_value = _to_unitful(value.unit)
    stored_unit = Unitful.unit(unit_value)
    return GPU.GeoUnit{typeof(value.val), typeof(stored_unit)}(
        value.val, stored_unit, value.isdimensional
    )
end

convert(::Type{GPU.GeoUnit}, value::DQ.AbstractQuantity) = GPU.GeoUnit(_to_unitful(value))
convert(::Type{GPU.GeoUnit}, value::AbstractArray{<:DQ.AbstractQuantity}) =
    GPU.GeoUnit(_to_unitful(value))
convert(::Type{GPU.GeoUnit{T}}, value::DQ.AbstractQuantity) where {T} =
    GPU.GeoUnit{T}(_to_unitful(value))

function convert(::Type{GPD.GeoUnit}, value::GPU.GeoUnit)
    dynamic_unit = _to_dynamic(1.0 * value.unit)
    return GPD.GeoUnit{typeof(value.val), typeof(dynamic_unit)}(
        value.val, dynamic_unit, value.isdimensional
    )
end

convert(::Type{GPD.GeoUnit}, value::Unitful.AbstractQuantity) = GPD.GeoUnit(_to_dynamic(value))
convert(::Type{GPD.GeoUnit}, value::AbstractArray{<:Unitful.AbstractQuantity}) =
    GPD.GeoUnit(_to_dynamic(value))

@inline _converted_data(converter, units) = map(converter, getfield(units, :data))

convert(::Type{GPU.GeoUnits}, units::GPD.GeoUnits{GPD.GEO}) =
    GPU.GeoUnits{GPU.GEO}(; _converted_data(_to_unitful, units)...)
convert(::Type{GPU.GeoUnits}, units::GPD.GeoUnits{GPD.SI}) =
    GPU.GeoUnits{GPU.SI}(; _converted_data(_to_unitful, units)...)
convert(::Type{GPU.GeoUnits}, units::GPD.GeoUnits{GPD.NONE}) =
    GPU.GeoUnits{GPU.NONE}(; _converted_data(_to_unitful, units)...)

convert(::Type{GPD.GeoUnits}, units::GPU.GeoUnits{GPU.GEO}) =
    GPD.GeoUnits{GPD.GEO}(; _converted_data(_to_dynamic, units)...)
convert(::Type{GPD.GeoUnits}, units::GPU.GeoUnits{GPU.SI}) =
    GPD.GeoUnits{GPD.SI}(; _converted_data(_to_dynamic, units)...)
convert(::Type{GPD.GeoUnits}, units::GPU.GeoUnits{GPU.NONE}) =
    GPD.GeoUnits{GPD.NONE}(; _converted_data(_to_dynamic, units)...)

const DynamicInput = Union{
    DQ.AbstractQuantity,
    GPD.GeoUnit,
    AbstractArray{<:DQ.AbstractQuantity},
    NTuple{N, Union{DQ.AbstractQuantity, GPD.GeoUnit}} where {N},
}

nondimensionalize(value::DynamicInput, units::GPD.GeoUnits) =
    GPD.nondimensionalize(value, units)
compute_units(value::GPD.GeoUnit, units::GPD.GeoUnits) = GPD.compute_units(value, units)
dimensionalize(value, target::DQ.AbstractQuantity, units::GPD.GeoUnits) =
    GPD.dimensionalize(value, target, units)
dimensionalize(value, target::DQ.AffineUnit, units::GPD.GeoUnits) =
    GPD.dimensionalize(value, target, units)
dimensionalize(value, target::GPD.DynamicUnitSymbol, units::GPD.GeoUnits) =
    GPD.dimensionalize(value, target, units)
dimensionalize(value::GPD.GeoUnit, units::GPD.GeoUnits) = GPD.dimensionalize(value, units)
dimensionalize_and_strip(value, target::DQ.AbstractQuantity, units::GPD.GeoUnits) =
    GPD.ustrip(GPD.dimensionalize(value, target, units))
dimensionalize_and_strip(value, target::DQ.AffineUnit, units::GPD.GeoUnits) =
    GPD.ustrip(GPD.dimensionalize(value, target, units))
dimensionalize_and_strip(value, target::GPD.DynamicUnitSymbol, units::GPD.GeoUnits) =
    GPD.ustrip(GPD.dimensionalize(value, target, units))

nondimensionalize(material::GeoParams.AbstractMaterialParam, units::GPD.GeoUnits) =
    GeoParams.nondimensionalize(material, convert(GPU.GeoUnits, units))
nondimensionalize(material::GeoParams.AbstractMaterialParamsStruct, units::GPD.GeoUnits) =
    GeoParams.nondimensionalize(material, convert(GPU.GeoUnits, units))
nondimensionalize(
    materials::NTuple{N, GeoParams.AbstractMaterialParamsStruct}, units::GPD.GeoUnits
) where {N} = GeoParams.nondimensionalize(materials, convert(GPU.GeoUnits, units))

dimensionalize(material::GeoParams.AbstractMaterialParam, units::GPD.GeoUnits) =
    GeoParams.dimensionalize(material, convert(GPU.GeoUnits, units))
dimensionalize(material::GeoParams.AbstractMaterialParamsStruct, units::GPD.GeoUnits) =
    GeoParams.dimensionalize(material, convert(GPU.GeoUnits, units))
dimensionalize(
    materials::NTuple{N, GeoParams.AbstractMaterialParamsStruct}, units::GPD.GeoUnits
) where {N} = GeoParams.dimensionalize(materials, convert(GPU.GeoUnits, units))

end
