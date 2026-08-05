Base.length(value::GeoUnit) = length(value.val)
Base.size(value::GeoUnit) = size(value.val)
Base.eltype(::Type{GeoUnit{T, U}}) where {T, U} = GeoUnit{eltype(T), U}
Base.IteratorSize(::Type{<:GeoUnit}) = Base.HasLength()
Base.getindex(value::GeoUnit, indices::Vararg{Int, N}) where {N} =
    GeoUnit(value.val[indices...], value.unit, value.isdimensional)
Base.iterate(value::GeoUnit, state = 1) =
    state > length(value) ? nothing : (GeoUnit(value.val[state], value.unit, value.isdimensional), state + 1)

for op in (:+, :-, :*, :/)
    @eval Base.$op(x::GeoUnit, y::Number) = $(op)(x.val, y)
    @eval Base.$op(x::Number, y::GeoUnit) = $(op)(x, y.val)

    @eval function Base.$op(x::GeoUnit, y::GeoUnit)
        dimensional = x.isdimensional && y.isdimensional
        new_value = $(op)(upreferred.(UnitValue(x)), upreferred.(UnitValue(y)))
        new_unit = if dimensional
            unit(new_value)
        elseif ($(op) == +) || ($(op) == -)
            Unit(x)
        else
            $(op)(Unit(x), Unit(y))
        end
        return GeoUnit(ustrip.(new_value), new_unit, dimensional)
    end

    @eval Base.$op(x::GeoUnit, y::DynamicQuantity) =
        $(op).(upreferred.(UnitValue(x)), upreferred.(y))
    @eval Base.$op(x::DynamicQuantity, y::GeoUnit) =
        $(op).(upreferred.(x), upreferred.(UnitValue(y)))

    @eval Base.$op(x::GeoUnit, y::AbstractArray) = broadcast($op, NumValue(x), y)
    @eval Base.$op(x::AbstractArray, y::GeoUnit) = broadcast($op, x, NumValue(y))
    @eval Base.$op(x::GeoUnit, y::AbstractArray{<:DynamicQuantity}) =
        broadcast($op, NumValue(x), y)
    @eval Base.$op(x::AbstractArray{<:DynamicQuantity}, y::GeoUnit) =
        broadcast($op, x, Value(y))

    @eval Base.broadcasted(::typeof($(op)), x::GeoUnit, y::AbstractArray) =
        broadcast($(op), NumValue(x), y)
    @eval Base.broadcasted(::typeof($(op)), x::AbstractArray, y::GeoUnit) =
        broadcast($(op), x, NumValue(y))
    @eval Base.broadcasted(
        ::typeof($(op)), x::GeoUnit, y::AbstractArray{<:DynamicQuantity}
    ) = broadcast($(op), Value(x), y)
    @eval Base.broadcasted(
        ::typeof($(op)), x::AbstractArray{<:DynamicQuantity}, y::GeoUnit
    ) = broadcast($(op), x, Value(y))

    @eval function Base.broadcasted(::typeof($(op)), x::GeoUnit, y::GeoUnit)
        dimensional = x.isdimensional && y.isdimensional
        new_value = broadcast($(op), upreferred.(UnitValue(x)), upreferred.(UnitValue(y)))
        new_unit = if dimensional
            unit(new_value)
        elseif ($(op) == +) || ($(op) == -)
            Unit(x)
        else
            $(op)(Unit(x), Unit(y))
        end
        return GeoUnit(ustrip.(new_value), new_unit, dimensional)
    end
end

Base.broadcasted(::typeof(^), x::GeoUnit, y::AbstractArray) = broadcast(^, NumValue(x), y)
Base.broadcasted(::typeof(^), x::AbstractArray, y::GeoUnit) = broadcast(^, x, NumValue(y))
Base.broadcasted(::typeof(^), x::GeoUnit, y::AbstractArray{<:DynamicQuantity}) =
    broadcast(^, Value(x), y)
Base.broadcasted(::typeof(^), x::AbstractArray{<:DynamicQuantity}, y::GeoUnit) =
    broadcast(^, x, Value(y))

function Base.setindex!(value::GeoUnit, replacement, indices::Vararg{Int, N}) where {N}
    return value.val[indices...] = replacement
end
