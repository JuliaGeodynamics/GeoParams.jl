"""
    upgrade_GeoUnits(a)

Function which recreates an outdated `GeoUnits` object in line with the current format
"""
function upgrade_GeoUnits(a)
    return GEO_units(length = a.length, temperature = a.temperature, stress = a.stress, viscosity = a.viscosity)
end
