# define a few basic routines so we can easily operate with GeoUnits
Base.length(v::GeoUnit) = length(v.val)
Base.size(v::GeoUnit) = size(v.val)
Base.eltype(::Type{GeoUnit{T, U}}) where {T, U} = GeoUnit{eltype(T), U}
Base.IteratorSize(::Type{<:GeoUnit}) = Base.HasLength()
Base.getindex(A::GeoUnit{T, U}, inds::Vararg{Int, N}) where {T, U, N} = GeoUnit(A.val[inds...], A.unit, A.isdimensional)
Base.iterate(S::GeoUnit, state = 1) = state > length(S) ? nothing : (GeoUnit(S.val[state], S.unit, S.isdimensional), state + 1)

for op in (:+, :-, :*, :/)

    # Multiply with number
    @eval Base.$op(x::GeoUnit, y::Number) = $(op)(x.val, y)
    @eval Base.$op(x::Number, y::GeoUnit) = $(op)(x, y.val)

    # Multiplying/dividing/adding/subtracting a GeoUnit with another one, returns a GeoUnit
    @eval function Base.$op(x::GeoUnit, y::GeoUnit)
        isdimensional_new_unit = x.isdimensional * y.isdimensional
        new_value = $(op)(upreferred.(UnitValue(x)), upreferred.(UnitValue(y)))
        if isdimensional_new_unit
            new_unit = unit(new_value)
        else
            if ($(op) == +) ||  ($(op) == -)
                new_unit = Unit(x)
            else
                new_unit = $(op)(Unit(x), Unit(y))
            end
        end

        return GeoUnit(ustrip.(new_value), new_unit, isdimensional_new_unit)
    end
    @eval Base.$op(x::GeoUnit, y::Quantity) = $(op).(UnitValue(x), y)
    @eval Base.$op(x::Quantity, y::GeoUnit) = $(op).(x, UnitValue(y))

    # If we multiply a GeoUnit with an abstract array, we only return values, not units (use GeoUnits for that)
    @eval Base.$op(x::GeoUnit, y::AbstractArray) = broadcast($op, NumValue(x), y)
    @eval Base.$op(x::AbstractArray, y::GeoUnit) = broadcast($op, x, NumValue(y))
    @eval Base.$op(x::GeoUnit, y::AbstractArray{<:Quantity}) = broadcast($op, NumValue(x), y)

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

    @eval function Base.broadcasted(
            ::typeof($(op)), A::GeoUnit, B::GeoUnit
        )
        isdimensional_new = A.isdimensional * B.isdimensional
        new_value = broadcast($(op), upreferred.(UnitValue(A)), upreferred.(UnitValue(B)))
        if isdimensional_new
            new_unit = new_value isa Quantity ? unit(new_value) : unit(first(new_value))
        else
            if ($(op) == +) || ($(op) == -)
                new_unit = Unit(A)
            else
                new_unit = $(op)(Unit(A), Unit(B))
            end
        end
        return GeoUnit(ustrip.(new_value), new_unit, isdimensional_new)
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
function Base.setindex!(A::GeoUnit{T, U}, val, inds::Vararg{Int, N}) where {T, U, N}
    return A.val[inds...] = val
end
