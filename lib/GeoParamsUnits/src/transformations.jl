"""
    nondimensionalize(param, CharUnits::GeoUnits{TYPE})

Nondimensionalizes `param` using the characteristic values specified in `CharUnits`

# Example 1
```julia-repl
julia> using GeoParamsUnits;
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
function nondimensionalize(param::GeoUnit{T, U}, g::GeoUnits{TYPE}) where {T, U, TYPE}

    if param.isdimensional
        char_val = compute_units(param, g)
        val_ND = upreferred.(param.val * param.unit) / char_val
        val_ND_number = convert(T, val_ND)
        param_ND = GeoUnit{T, U}(val_ND_number, param.unit, false)          # store new value, but keep original dimensions

    else
        param_ND = param
    end
    return param_ND
end

# in case the parameter is already non-dimensional:
nondimensionalize(param::String, ::GeoUnits{TYPE}) where {TYPE} = param
nondimensionalize(param::String, ::Nothing) = param
nondimensionalize(param::Ptr{UInt8}, ::GeoUnits{TYPE}) where {TYPE} = unsafe_string(param)
nondimensionalize(param::Ptr{UInt8}, ::Nothing) = unsafe_string(param)

# in case it is a unitful quantity
function nondimensionalize(param::Unitful.Quantity, g::GeoUnits)
    param_Geo = GeoUnit(param)
    result = nondimensionalize(param_Geo, g)
    return NumValue(result)
end

function nondimensionalize(
        param::Union{Unitful.Quantity, AbstractArray{<:Quantity}},
        g::Nothing,
    )
    @warn("The input parameter is not being nondimensionalized, as no characteristic units are given")
    return ustrip.(param)
end

function nondimensionalize(param::GeoUnit, ::Nothing)
    @warn("The input parameter is not being nondimensionalized, as no characteristic units are given")
    return param
end

function nondimensionalize(
        ::Union{Integer, AbstractFloat, AbstractArray{<:Integer}, AbstractArray{<:AbstractFloat}, AbstractArray{Any}},
        ::GeoUnits,
    )
    throw(ArgumentError("The input parameter should have units"))
end

function nondimensionalize(
        ::NTuple{N, Union{Integer, AbstractFloat, AbstractArray{<:Integer}, AbstractArray{<:AbstractFloat}, AbstractArray{Any}}},
        ::GeoUnits,
    ) where {N}
    throw(ArgumentError("The input parameter should have units"))
end

# If it is an array, but has no units we cannot know how to nondimensionalize it
nondimensionalize(param::AbstractArray{<:Number}, g::GeoUnits) = param

function nondimensionalize(
        param::AbstractArray{<:Quantity},
        g::GeoUnits
    )
    param_Geo = GeoUnit.(param)
    result = map(p -> nondimensionalize(p, g), param_Geo)
    return NumValue.(result)
end

function nondimensionalize(param::AbstractArray{<:Number}, ::Nothing)
    @warn("The input parameter is not being nondimensionalized, as no characteristic units are given")
    return param
end

"""
    param = nondimensionalize(param::NTuple{N,Quantity}, g::GeoUnits{TYPE})

nondimensionalizes a tuple of parameters
"""
function nondimensionalize(
        param::NTuple{N, Union{Quantity, GeoUnit}}, g::Union{GeoUnits, Nothing}
    ) where {N}
    return ntuple(Val(N)) do i
        Base.@_inline_meta
        nondimensionalize(param[i], g)
    end
end

nondimensionalize(::Tuple{}, ::GeoUnits) = ()

nondimensionalize(args...) = nondimensionalize(Tuple(args[1:(end - 1)]), args[end])

@inline dimension_types(::Unitful.Dimensions{T}) where {T} = T

# This computes the characteristic value
function compute_units(
        param::GeoUnit{<:Union{T, AbstractArray{T}}, U}, g::GeoUnits
    ) where {T, U}
    dim = Unitful.dimension(param.unit)              # Basic SI units
    value::T = if dim == NoDims
        one(T)
    else
        prod(dimension_types(dim)) do y
            val::T = upreferred(getproperty(g, Unitful.name(y))).val      # Retrieve the characteristic value from structure g
            pow = Float64(y.power)                                  # power by which it should be multiplied
            val^pow                                     # multiply characteristic value
        end
    end
    char_val_out = upreferred(param.unit) * value               # this is done for type-stability

    return char_val_out
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
        param_ND, param_dim::Unitful.FreeUnits, g::GeoUnits
    )
    char_val = compute_units(GeoUnit(1.0 * param_dim), g)         # Determine characteristic units
    param = uconvert.(param_dim, param_ND * char_val)
    return param
end

function dimensionalize(
        param_ND, ::Unitful.FreeUnits, ::Nothing
    )
    @warn("The input parameter is not being dimensionalized, as no characteristic units are given")
    return ustrip.(param_ND)
end

function dimensionalize(param_ND::GeoUnit{T, U}, g::GeoUnits) where {T, U}
    param = if isdimensional(param_ND) == false
        char_val = compute_units(param_ND, g)                   # Determine characteristic units
        val = uconvert.(param_ND.unit, param_ND.val * char_val) # dimensionalize
        GeoUnit{T, U}(ustrip.(val), param_ND.unit, true)        # store new value, but keep original dimensions
    else
        param_ND
    end
    return param
end

@inline udim(args::Vararg{Any, N}) where {N} = ustrip(dimensionalize(args...))

"""
    dimensionalize_and_strip(args...)

Converts the input arguments to their dimensional (unitful) form using `dimensionalize`,
and then strips the units using `ustrip`. This function is useful when you want to
apply units to values and immediately retrieve their plain numeric values.
"""
dimensionalize_and_strip(args::Vararg{Any, N}) where {N} = ustrip(dimensionalize(args...))

macro dimstrip(args...)
    return quote
        dimensionalize_and_strip($(esc.(args)...))
    end
end
