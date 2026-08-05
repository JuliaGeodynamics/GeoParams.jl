"""
    GeoUnits{TYPE}

Characteristic values used for nondimensionalization.
"""
struct GeoUnits{TYPE, F <: NamedTuple}
    data::F
end

@inline function Base.getproperty(units::GeoUnits, name::Symbol)
    name === :data && return getfield(units, :data)
    return getproperty(getfield(units, :data), name)
end

@inline Base.propertynames(units::GeoUnits, private::Bool = false) =
    propertynames(getfield(units, :data), private)

function GeoUnits{TYPE}(;
        temperature = 1,
        length = 1,
        stress = 1,
        time = 1,
        viscosity = 1,
        K = 1,
        s = 1,
        m = 1,
        Pa = 1,
        kg = upreferred(Pa * m * s^2),
        Length = m,
        Mass = kg,
        Time = s,
        Temperature = K,
        Amount = 1mol,
        Second = s,
        N = kg * m / s^2,
        J = N * m,
        W = J / s,
        area = m^2,
        volume = m^3,
        velocity = m / s,
        density = kg / m^3,
        acceleration = m / s^2,
        force = kg * m / s^2,
        strainrate = 1 / s,
        heatcapacity = J / kg / K,
        conductivity = W / m / K,
        chemicaldiffusion = m^2 / s,
        SecYear = 3600 * 24 * 365.25,
        Myr = 1.0e6,
        cmYear = SecYear * 100,
    ) where {TYPE}
    data = (;
        temperature,
        length,
        stress,
        time,
        viscosity,
        K,
        s,
        m,
        Pa,
        kg,
        Length,
        Mass,
        Time,
        Temperature,
        Amount,
        Second,
        N,
        J,
        W,
        area,
        volume,
        velocity,
        density,
        acceleration,
        force,
        strainrate,
        heatcapacity,
        conductivity,
        chemicaldiffusion,
        SecYear,
        Myr,
        cmYear,
    )
    return GeoUnits{TYPE, typeof(data)}(data)
end

@inline _with_default_unit(value::DynamicQuantity, ::DynamicQuantity) = value
@inline _with_default_unit(value::DynamicQuantity, ::DynamicQuantities.AffineUnit) = value
@inline _with_default_unit(value::DynamicQuantity, ::DynamicUnitSymbol) = value
@inline _with_default_unit(value::Number, default_unit::DynamicQuantity) = value * default_unit
@inline _with_default_unit(value::Number, default_unit::DynamicQuantities.AffineUnit) =
    value * default_unit
@inline _with_default_unit(value::Number, default_unit::DynamicUnitSymbol) =
    value * default_unit

"""
    GEO_units(; length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)

Create geodynamics-oriented characteristic units.
"""
function GEO_units(;
        length = 1000km,
        temperature = 1000C,
        stress = 10MPa,
        viscosity = 1.0e20Pas,
    )
    temperature = _with_default_unit(temperature, C)
    length = _with_default_unit(length, km)
    stress = _with_default_unit(stress, MPa)
    viscosity = _with_default_unit(viscosity, Pas)

    T = _convert_unit(C, temperature)
    Le = _convert_unit(km, length)
    Sigma = _convert_unit(MPa, stress)
    Eta = _convert_unit(Pas, viscosity)

    T_SI = upreferred(T)
    Le_SI = upreferred(Le)
    Sigma_SI = upreferred(Sigma)
    Time_SI = upreferred(Eta) / Sigma_SI
    time = _convert_unit(Myr, Time_SI)

    return GeoUnits{GEO}(;
        length = Le,
        temperature = T,
        stress = Sigma,
        viscosity = Eta,
        time,
        m = Le_SI,
        K = T_SI,
        Pa = Sigma_SI,
        s = Time_SI,
    )
end

"""
    SI_units(; length=1000m, temperature=1000K, stress=10Pa, viscosity=1e20Pas)

Create SI-oriented characteristic units.
"""
function SI_units(;
        length = 1000m,
        temperature = 1000K,
        stress = 10Pa,
        viscosity = 1.0e20Pas,
    )
    temperature = _with_default_unit(temperature, K)
    length = _with_default_unit(length, m)
    stress = _with_default_unit(stress, Pa)
    viscosity = _with_default_unit(viscosity, Pas)

    T = _convert_unit(K, temperature)
    Le = _convert_unit(m, length)
    Sigma = _convert_unit(Pa, stress)
    Eta = _convert_unit(Pas, viscosity)

    T_SI = upreferred(T)
    Le_SI = upreferred(Le)
    Sigma_SI = upreferred(Sigma)
    Time_SI = upreferred(Eta) / Sigma_SI
    time = _convert_unit(s, Time_SI)

    return GeoUnits{SI}(;
        length = Le,
        temperature = T,
        stress = Sigma,
        viscosity = Eta,
        time,
        m = Le_SI,
        K = T_SI,
        Pa = Sigma_SI,
        s = Time_SI,
    )
end

"""
    NO_units(; length=1, temperature=1, stress=1, viscosity=1)

Create dimensionless characteristic values.
"""
function NO_units(; length = 1, temperature = 1, stress = 1, viscosity = 1)
    for (name, value) in pairs((; temperature, length, stress, viscosity))
        _isquantity(value) && throw(ArgumentError("$name should not have units"))
    end

    time = viscosity / stress
    return GeoUnits{NONE}(;
        length,
        temperature,
        stress,
        viscosity,
        time,
        m = length,
        K = temperature,
        Pa = stress,
        s = time,
    )
end
