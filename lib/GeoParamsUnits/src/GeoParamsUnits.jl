"""
    This provides units and creates a non-dimensionalization object
"""
module GeoParamsUnits

using Unitful

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
include("geo_unit.jl")
include("arithmetic.jl")
include("characteristic_units.jl")
include("transformations.jl")
include("unpack.jl")
include("display.jl")
include("compat.jl")

end
