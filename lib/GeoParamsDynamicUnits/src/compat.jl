"""
    upgrade_GeoUnits(units)

Recreate a legacy characteristic-unit object using the current representation.
"""
function upgrade_GeoUnits(units)
    return GEO_units(;
        length = units.length,
        temperature = units.temperature,
        stress = units.stress,
        viscosity = units.viscosity,
    )
end
