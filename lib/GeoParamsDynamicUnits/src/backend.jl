const DynamicQuantity = DynamicQuantities.AbstractQuantity
const DynamicUnit = Union{DynamicQuantity, DynamicQuantities.AffineUnit, DynamicUnitSymbol}

@inline _isquantity(value) = value isa DynamicQuantity
@inline _isunit(value) = value isa DynamicUnit
@inline _isdimensionless(value) = iszero(dimension(value))

@inline unit(value::DynamicQuantity) =
    uexpand(Quantity(one(ustrip(value)), dimension(value)))
@inline unit(values::AbstractArray{<:DynamicQuantity}) = unit(_representative_value(values))
@inline unit(::Number) = NoUnits

@inline upreferred(value::DynamicQuantity) = uexpand(value)
@inline upreferred(value) = value

@inline function _normalized_unit(target::DynamicQuantity)
    return Quantity(one(ustrip(target)), dimension(target))
end

@inline _normalized_unit(target::DynamicUnitSymbol) = _materialize(target)

function _convert_unit(target::DynamicQuantity, value::DynamicQuantity)
    if dimension(target) isa DynamicQuantities.AbstractSymbolicDimensions
        return uconvert(_normalized_unit(target), uexpand(value))
    end
    dimension(target) == dimension(value) ||
        throw(DynamicQuantities.DimensionError(target, value))
    return uexpand(value)
end

function _convert_unit(target::DynamicQuantities.AffineUnit, value::DynamicQuantity)
    return ustrip(target, uexpand(value)) * target
end


@inline _convert_unit(target::DynamicUnitSymbol, value::DynamicQuantity) =
    uconvert(_materialize(target), uexpand(value))

@inline _convert_unit(::Number, value::DynamicQuantity) = ustrip(value)

@inline _expanded_value(value::DynamicQuantity) = ustrip(uexpand(value))

function _quantity_with_dimensions(value, dimensions::DynamicQuantities.AbstractDimensions)
    return Quantity(value, dimensions)
end

function _unit_for_dimensions(target::DynamicQuantity)
    return _normalized_unit(target)
end

function _unit_for_dimensions(target::DynamicQuantities.AffineUnit)
    return Quantity(1.0, target.basedim)
end


@inline _unit_for_dimensions(target::DynamicUnitSymbol) = _materialize(target)


@inline _unit_for_dimensions(::Number) = NoUnits
