"""
DynamicQuantities-backed unit handling and nondimensionalization for GeoParams.
"""
module GeoParamsDynamicUnits

using DynamicQuantities

import Base:
    show, isapprox, isequal, convert, length, size, getindex, setindex!, getproperty, iterate
import Base.Broadcast: broadcasted

export km,
    m,
    cm,
    mm,
    μm,
    Myr,
    yr,
    s,
    GPa,
    MPa,
    Pa,
    bar,
    kbar,
    Pas,
    K,
    C,
    g,
    kg,
    mol,
    J,
    kJ,
    Watt,
    μW,
    GeoUnit,
    GeoUnits,
    GEO_units,
    SI_units,
    NO_units,
    AbstractGeoUnit,
    nondimensionalize,
    dimensionalize,
    dimensionalize_and_strip,
    @dimstrip,
    superscript,
    upreferred,
    GEO,
    SI,
    NONE,
    isDimensional,
    Value,
    Fun,
    NumValue,
    unpack_units,
    unpack_vals,
    Unit,
    UnitValue,
    isdimensional,
    compute_units,
    udim,
    upgrade_GeoUnits

include("unit_constants.jl")
include("backend.jl")
include("geo_unit.jl")
include("arithmetic.jl")
include("characteristic_units.jl")
include("transformations.jl")
include("unpack.jl")
include("display.jl")
include("compat.jl")

function _register_dynamic_unit(name::Symbol, value::DynamicQuantity)
    name in DynamicQuantities.ALL_SYMBOLS && return nothing
    push!(DynamicQuantities.Units.UNIT_SYMBOLS, name)
    push!(DynamicQuantities.Units.UNIT_VALUES, value)
    DynamicQuantities.update_all_values(name, value)
    return nothing
end

function __init__()
    _register_dynamic_unit(:MPa, 1.0e6 * DynamicQuantities.Units.Pa)
    _register_dynamic_unit(:GPa, 1.0e9 * DynamicQuantities.Units.Pa)
    _register_dynamic_unit(:kbar, 1.0e3 * DynamicQuantities.Units.bar)
    _register_dynamic_unit(:μW, 1.0e-6 * DynamicQuantities.Units.W)
    _MPa_UNIT[] = Quantity(1.0, SymbolicDimensions(; MPa = 1))
    _GPa_UNIT[] = Quantity(1.0, SymbolicDimensions(; GPa = 1))
    _KBAR_UNIT[] = Quantity(1.0, SymbolicDimensions(; kbar = 1))
    _MICROWATT_UNIT[] = Quantity(1.0, SymbolicDimensions(; μW = 1))
    return nothing
end

end
