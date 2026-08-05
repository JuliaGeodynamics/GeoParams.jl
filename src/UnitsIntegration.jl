using Setfield

import .Units: dimensionalize, isDimensional, nondimensionalize

"""
    nondimensionalize(param::AbstractMaterialParam, units)

Nondimensionalize every `GeoUnit` field in a GeoParams material parameter.
"""
@generated function nondimensionalize(
        material::T, units::Union{GeoUnits, Nothing}
    ) where {T <: AbstractMaterialParam}
    fields = fieldnames(material)
    nfields = length(fields)
    return quote
        values = Base.@ntuple $nfields i ->
        _nondimensionalize(getfield(material, $fields[i]), units)
        strip_type_parameters(material)(values...)
    end
end

@inline strip_type_parameters(::T) where {T} = Base.typename(T).wrapper

@inline function _nondimensionalize(
        value::Union{GeoUnit, AbstractMaterialParam}, units::Union{GeoUnits, Nothing}
    )
    return nondimensionalize(value, units)
end

@inline function _nondimensionalize(
        value::NTuple{N, GeoUnit}, units::Union{GeoUnits, Nothing}
    ) where {N}
    return ntuple(i -> nondimensionalize(value[i], units), Val(N))
end

@inline _nondimensionalize(value, ::Any) = value

"""
    nondimensionalize(material::AbstractMaterialParamsStruct, units)

Nondimensionalize all material-parameter fields in a phase definition.
"""
@generated function nondimensionalize(
        material::AbstractMaterialParamsStruct, units::Union{GeoUnits, Nothing}
    )
    fields = fieldnames(material)
    nfields = length(fields)
    return quote
        Base.@_inline_meta
        Base.@nexprs $nfields i ->
        material = begin
            fieldname = $(fields)[i]
            value = getfield(material, fieldname)
            nondimensionalize_MatParam(value, material, fieldname, units)
        end
        material = set(material, Setfield.PropertyLens{:Nondimensional}(), true)
        return material
    end
end

@inline nondimensionalize_MatParam(::Any, material, ::Vararg{Any, N}) where {N} = material
@inline nondimensionalize_MatParam(::Tuple{}, material, ::Vararg{Any, N}) where {N} = material

@inline _nondimensionalize_MatParam(value::AbstractMaterialParam, units) =
    nondimensionalize(value, units)
@inline _nondimensionalize_MatParam(value::AbstractPhaseDiagramsStruct, units) =
    PerpleX_LaMEM_Diagram(unsafe_string(value.Name); CharDim = units)
@inline _nondimensionalize_MatParam(value, units) = value

function nondimensionalize_MatParam(
        values::NTuple{N, AbstractMaterialParam}, material, fieldname, units
    ) where {N}
    new_values = ntuple(Val(N)) do i
        @inline
        _nondimensionalize_MatParam(values[i], units)
    end
    return set(material, Setfield.PropertyLens{fieldname}(), new_values)
end

function nondimensionalize(
        materials::NTuple{N, AbstractMaterialParamsStruct}, units::GeoUnits
    ) where {N}
    return ntuple(i -> nondimensionalize(materials[i], units), Val(N))
end

"""
    dimensionalize(param::AbstractMaterialParam, units)

Dimensionalize every `GeoUnit` field in a GeoParams material parameter.
"""
@generated function dimensionalize(
        material::AbstractMaterialParam, units::Union{GeoUnits, Nothing}
    )
    fields = fieldnames(material)
    nfields = length(fields)
    return quote
        Base.@_inline_meta
        Base.@nexprs $nfields i ->
        material = begin
            fieldname = $(fields)[i]
            _dimensionalize(material, getfield(material, fieldname), units, fieldname)
        end
        return material
    end
end

@inline function _dimensionalize(
        material,
        value::Union{GeoUnit, AbstractMaterialParam},
        units::Union{GeoUnits, Nothing},
        fieldname,
    )
    new_value = dimensionalize(value, units)
    return set(material, Setfield.PropertyLens{fieldname}(), new_value)
end

@inline _dimensionalize(material, ::Any, ::Any, ::Any) = material

"""
    dimensionalize(material::AbstractMaterialParamsStruct, units)

Dimensionalize all material-parameter fields in a phase definition.
"""
@generated function dimensionalize(
        material::AbstractMaterialParamsStruct, units::Union{GeoUnits, Nothing}
    )
    fields = fieldnames(material)
    nfields = length(fields)
    return quote
        Base.@_inline_meta
        Base.@nexprs $nfields i ->
        material = begin
            fieldname = $(fields)[i]
            value = getfield(material, fieldname)
            dimensionalize_MatParam(value, material, fieldname, units)
        end
        material = set(material, Setfield.PropertyLens{:Nondimensional}(), false)
        return material
    end
end

@inline dimensionalize_MatParam(::Tuple{}, material, ::Vararg{Any, N}) where {N} = material
@inline dimensionalize_MatParam(::Any, material, ::Vararg{Any, N}) where {N} = material

@inline _dimensionalize_MatParam(value::AbstractMaterialParam, units) =
    dimensionalize(value, units)
@inline _dimensionalize_MatParam(value, ::Any) = value

@inline function dimensionalize_MatParam(
        values::NTuple{N, AbstractMaterialParam}, material, fieldname, units
    ) where {N}
    new_values = ntuple(Val(N)) do i
        @inline
        _dimensionalize_MatParam(values[i], units)
    end
    return set(material, Setfield.PropertyLens{fieldname}(), new_values)
end

@inline function dimensionalize_MatParam(
        value::AbstractPhaseDiagramsStruct, material, fieldname, units
    )
    phase_diagram = PerpleX_LaMEM_Diagram(unsafe_string(value.Name); CharDim = units)
    return set(material, Setfield.PropertyLens{fieldname}(), (phase_diagram,))
end

function dimensionalize(
        materials::NTuple{N, AbstractMaterialParamsStruct},
        units::Union{GeoUnits, Nothing},
    ) where {N}
    return ntuple(i -> dimensionalize(materials[i], units), Val(N))
end

"""
    isDimensional(param::AbstractMaterialParam)

Return `true` when at least one direct `GeoUnit` field is dimensional.
"""
function isDimensional(material::AbstractMaterialParam)
    for fieldname in fieldnames(typeof(material))
        value = getfield(material, fieldname)
        value isa GeoUnit && value.isdimensional && return true
    end
    return false
end
