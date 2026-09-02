using UnPack
export @unpack_val, @unpack_units
export precision_of, convert_precision

"""
    precision_of(x)

Scalar type in which a calculation seeded by `x` should be evaluated. Units are
stripped and arrays and tuples report their element type.

Only a float selects its own type. Anything else — an integer, or a number that
wraps a value for automatic differentiation — expresses no precision preference
and selects `Float64`, which leaves the stored parameters alone and lets
ordinary promotion carry the result back into the caller's type. Forcing
parameters into a tracked number instead would put constants on the AD tape.
"""
@inline precision_of(x) = _precision_of(ustrip(x))
@inline precision_of(x::AbstractArray) = precision_of(zero(eltype(x)))
@inline precision_of(x::Tuple) = precision_of(zero(eltype(x)))

@inline _precision_of(::T) where {T <: AbstractFloat} = T
@inline _precision_of(::Any) = Float64
# A dual number is differentiated at the precision of the value it carries.
@inline _precision_of(d::Dual) = _precision_of(value(d))

"""
    convert_precision(T, x)

Convert `x` to evaluation type `T`, retaining units. Only plain real numbers are
converted: arrays, and numbers such as dual numbers that carry more than a value,
pass through so that they keep promoting the result to their own type.
"""
@inline convert_precision(::Type{T}, x) where {T} = x
@inline convert_precision(::Type{T}, x::AbstractFloat) where {T} = convert(T, x)
@inline convert_precision(::Type{T}, x::Integer) where {T} = convert(T, x)
@inline convert_precision(::Type{T}, x::Quantity) where {T} =
    convert_precision(T, ustrip(x)) * unit(x)
@inline convert_precision(::Type{T}, x::Tuple) where {T} =
    map(y -> convert_precision(T, y), x)

# Element-wise conversion of a `GeoUnit` payload, which may be an array.
@inline _val_precision(::Type{T}, v::Number) where {T} = convert(T, v)
@inline _val_precision(::Type{T}, v) where {T} = T.(v)
@inline _val_precision(::Type{T}, v::AbstractArray{T}) where {T} = v

# Builds the body shared by all four macro forms. `withunits` multiplies the
# value back by its unit; `T` is `nothing` for the untyped forms.
function _unpack_geounit(args, withunits::Bool, T)
    args.head != :(=) && error("Expression needs to be of form `a, b = c`")
    items, suitecase = args.args
    items = isa(items, Symbol) ? [items] : items.args
    suitecase_instance = gensym()

    kd = map(items) do key
        field = :($UnPack.unpack($suitecase_instance, Val{$(Expr(:quote, key))}()))
        val = T === nothing ? :($field.val) : :($_val_precision($T, $field.val))
        rhs = withunits ? :($val .* $field.unit) : val
        return :($key = $rhs)
    end

    return quote
        local $suitecase_instance = $suitecase # handles if suitecase is an expression
        $(Expr(:block, kd...))
        $suitecase_instance # return RHS of `=` as standard in Julia
    end
end

"""
    @unpack_val ρ, α = r
    @unpack_val T ρ, α = r

Unpack the numerical values of `GeoUnit` fields, without their units. All
requested variables must be `GeoUnit`s.

The second form converts each value to evaluation type `T`, so a formula runs in
the precision of its solver state rather than that of the stored parameters.

This is a modification of the `@unpack` macro from `UnPack.jl`, which retrieves
the full variables.

# Example
```julia
julia> struct Density{T}
        ρ::GeoUnit{T}
        α::GeoUnit{T}
       end
julia> r = Density(GeoUnit(100kg/m^3),GeoUnit(4e-5/K));
julia> @unpack_val ρ,α = r
julia> α
4.0e-5
julia> @unpack_val Float32 ρ,α = r
julia> α
4.0f-5
```
"""
macro unpack_val(args)
    return esc(_unpack_geounit(args, false, nothing))
end

macro unpack_val(T, args)
    return esc(_unpack_geounit(args, false, T))
end

"""
    @unpack_units ρ, α = r
    @unpack_units T ρ, α = r

Unpack `GeoUnit` fields as `Quantity`s, retaining their units. All requested
variables must be `GeoUnit`s.

The second form gives each quantity numerical type `T`.

# Example
```julia
julia> struct Density{T}
        ρ::GeoUnit{T}
        α::GeoUnit{T}
       end
julia> r = Density(GeoUnit(100kg/m^3),GeoUnit(4e-5/K));
julia> @unpack_units ρ,α = r
julia> α
4.0e-5 K⁻¹·⁰
```
"""
macro unpack_units(args)
    return esc(_unpack_geounit(args, true, nothing))
end

macro unpack_units(T, args)
    return esc(_unpack_geounit(args, true, T))
end
