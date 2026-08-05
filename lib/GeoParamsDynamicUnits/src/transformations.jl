"""
    nondimensionalize(parameter, units)

Nondimensionalize a quantity or `GeoUnit` using the supplied characteristic units.
"""
function nondimensionalize(
        parameter::GeoUnit{T, U}, units::GeoUnits
    ) where {T <: Number, U}
    parameter.isdimensional || return parameter
    characteristic = compute_units(parameter, units)
    ratio = upreferred(Value(parameter)) / characteristic
    value = convert(T, ustrip(ratio))
    return GeoUnit{T, U}(value, parameter.unit, false)
end

function nondimensionalize(
        parameter::GeoUnit{A, U}, units::GeoUnits
    ) where {A <: AbstractArray, U}
    parameter.isdimensional || return parameter
    characteristic = compute_units(parameter, units)
    ratio = upreferred.(Value(parameter)) ./ characteristic
    value = convert(A, ustrip.(ratio))
    return GeoUnit{A, U}(value, parameter.unit, false)
end

nondimensionalize(parameter::String, ::GeoUnits) = parameter
nondimensionalize(parameter::String, ::Nothing) = parameter
nondimensionalize(parameter::Ptr{UInt8}, ::GeoUnits) = unsafe_string(parameter)
nondimensionalize(parameter::Ptr{UInt8}, ::Nothing) = unsafe_string(parameter)

function nondimensionalize(parameter::DynamicQuantity, units::GeoUnits)
    return NumValue(nondimensionalize(GeoUnit(parameter), units))
end

function nondimensionalize(
        parameter::Union{DynamicQuantity, AbstractArray{<:DynamicQuantity}},
        ::Nothing,
    )
    @warn "The input parameter is not being nondimensionalized, as no characteristic units are given"
    return ustrip.(parameter)
end

function nondimensionalize(parameter::GeoUnit, ::Nothing)
    @warn "The input parameter is not being nondimensionalized, as no characteristic units are given"
    return parameter
end

function nondimensionalize(
        ::Union{
            Integer,
            AbstractFloat,
            AbstractArray{<:Integer},
            AbstractArray{<:AbstractFloat},
            AbstractArray{Any},
        },
        ::GeoUnits,
    )
    throw(ArgumentError("The input parameter should have units"))
end

function nondimensionalize(
        ::NTuple{
            N,
            Union{
                Integer,
                AbstractFloat,
                AbstractArray{<:Integer},
                AbstractArray{<:AbstractFloat},
                AbstractArray{Any},
            },
        },
        ::GeoUnits,
    ) where {N}
    throw(ArgumentError("The input parameter should have units"))
end

nondimensionalize(parameter::AbstractArray{<:Number}, ::GeoUnits) = parameter

function nondimensionalize(parameter::AbstractArray{<:DynamicQuantity}, units::GeoUnits)
    return map(value -> nondimensionalize(value, units), parameter)
end

function nondimensionalize(parameter::AbstractArray{<:Number}, ::Nothing)
    @warn "The input parameter is not being nondimensionalized, as no characteristic units are given"
    return parameter
end

function nondimensionalize(
        parameter::NTuple{N, Union{DynamicQuantity, GeoUnit}},
        units::Union{GeoUnits, Nothing},
    ) where {N}
    return ntuple(i -> nondimensionalize(parameter[i], units), Val(N))
end

nondimensionalize(::Tuple{}, ::GeoUnits) = ()
nondimensionalize(args...) = nondimensionalize(Tuple(args[1:(end - 1)]), args[end])

const _DIMENSION_SCALES = (
    (:length, :m),
    (:mass, :kg),
    (:time, :s),
    (:temperature, :K),
    (:amount, :Amount),
)

@inline _storage_number_type(::Type{T}) where {T <: Number} = T
@inline _storage_number_type(::Type{T}) where {T <: AbstractArray} = eltype(T)

"""
    compute_units(parameter::GeoUnit, units::GeoUnits)

Compute the characteristic DynamicQuantities quantity for `parameter`.
"""
function compute_units(parameter::GeoUnit{T}, units::GeoUnits) where {T}
    return _compute_units(parameter, units, _storage_number_type(T))
end

function _compute_units(parameter::GeoUnit, units::GeoUnits, ::Type{T}) where {T}
    dimensions = dimension(upreferred(_unit_for_dimensions(parameter.unit)))
    value = one(T)

    for (dimension_name, scale_name) in _DIMENSION_SCALES
        power = getproperty(dimensions, dimension_name)
        iszero(power) && continue
        scale = convert(T, _expanded_value(getproperty(units, scale_name)))
        value = convert(T, value * scale^power)
    end

    unsupported = (:current, :luminosity)
    any(name -> !iszero(getproperty(dimensions, name)), unsupported) &&
        throw(ArgumentError("current and luminosity dimensions are not supported"))

    return _quantity_with_dimensions(value::T, dimensions)
end

"""
    dimensionalize(value, target_unit, units)

Dimensionalize a value into `target_unit` using the supplied characteristic units.
"""
function dimensionalize(value, target_unit::DynamicUnit, units::GeoUnits)
    characteristic = compute_units(GeoUnit(_unit_for_dimensions(target_unit)), units)
    dimensional = value .* characteristic
    return _convert_unit.(Ref(target_unit), dimensional)
end

function dimensionalize(value, ::DynamicUnit, ::Nothing)
    @warn "The input parameter is not being dimensionalized, as no characteristic units are given"
    return ustrip.(value)
end

function dimensionalize(parameter::GeoUnit{T, U}, units::GeoUnits) where {T, U}
    parameter.isdimensional && return parameter
    characteristic = compute_units(parameter, units)
    dimensional = parameter.val .* characteristic
    unit_scale = ustrip(parameter.unit)
    values = convert(T, ustrip.(dimensional) ./ unit_scale)
    return GeoUnit{T, U}(values, parameter.unit, true)
end

@inline udim(args::Vararg{Any, N}) where {N} = ustrip.(dimensionalize(args...))

"""
    dimensionalize_and_strip(args...)

Dimensionalize and return plain numeric values.
"""
dimensionalize_and_strip(args::Vararg{Any, N}) where {N} =
    ustrip.(dimensionalize(args...))

macro dimstrip(args...)
    return quote
        dimensionalize_and_strip($(esc.(args)...))
    end
end
