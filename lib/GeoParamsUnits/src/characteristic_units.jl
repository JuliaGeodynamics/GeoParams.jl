"""
    GeoUnits

    Structure that holds parameters used for non-dimensionalization
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

"""
    GEO_units(;length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)

Creates a non-dimensionalization object using GEO units.

GEO units implies that upon dimensionalization, `time` will be in `Myr`, `length` in `km`, stress in `MPa`, etc.
which is more convenient for typical geodynamic simulations than SI units
The characteristic values given as input can be in arbitrary units (`km` or `m`), provided the unit is specified.

# Examples:
```julia-repl
julia> CharUnits = GEO_units()
Employing GEO units
Characteristic values:
         length:      1000 km
         time:        0.3169 Myr
         stress:      10 MPa
         temperature: 1000.0 °C
julia> CharUnits.velocity
1.0e-7 m s⁻¹
```
If we instead have a crustal-scale simulation, it is likely more appropriate to use a different characteristic `length`:
```julia-repl
julia> CharUnits = GEO_units(length=10km)
Employing GEO units
Characteristic values:
         length:      10 km
         time:        0.3169 Myr
         stress:      10 MPa
         temperature: 1000.0 °C
```
"""
function GEO_units(; length = 1000km, temperature = 1000C, stress = 10MPa, viscosity = 1.0e20Pas)
    if unit(temperature) == NoUnits
        temperature = temperature * C
    end
    if unit(length) == NoUnits
        length = length * u"km"
    end
    if unit(stress) == NoUnits
        stress = stress * u"MPa"
    end
    if unit(viscosity) == NoUnits
        viscosity = viscosity * u"Pa*s"
    end

    T = uconvert(C, temperature)
    Le = uconvert(km, length)
    Sigma = uconvert(MPa, stress)
    Eta = uconvert(Pas, viscosity)

    T_SI = uconvert(K, T)
    Le_SI = uconvert(m, Le)
    Sigma_SI = uconvert(Pa, Sigma)
    Time_SI = Eta / Sigma_SI
    t = uconvert(Myr, Time_SI)

    return GeoUnits{GEO}(;
        length = Le,
        temperature = T,
        stress = Sigma,
        viscosity = Eta,
        time = t,
        m = Le_SI,
        K = T_SI,
        Pa = Sigma_SI,
        s = Time_SI,
    )
end

"""
    CharUnits = SI_units(length=1000m, temperature=1000K, stress=10Pa, viscosity=1e20)

Specify the characteristic values using SI units

# Examples:
```julia-repl
julia> CharUnits = SI_units(length=1000m)
Employing SI units
Characteristic values:
         length:      1000 m
         time:        1.0e19 s
         stress:      10 Pa
         temperature: 1000.0 K
```
Note that the same can be achieved if the input is given in `km`:
```julia-repl
julia> CharUnits = SI_units(length=1km)
```
"""
function SI_units(; length = 1000m, temperature = 1000K, stress = 10Pa, viscosity = 1.0e20Pas)
    if unit(temperature) == NoUnits
        temperature = temperature * K
    end
    if unit(length) == NoUnits
        length = length * u"m"
    end
    if unit(stress) == NoUnits
        stress = stress * u"Pa"
    end
    if unit(viscosity) == NoUnits
        viscosity = viscosity * u"Pa*s"
    end

    T = uconvert(K, temperature)
    Le = uconvert(m, length)
    Sigma = uconvert(Pa, stress)
    Eta = uconvert(Pas, viscosity)

    T_SI = uconvert(K, T)
    Le_SI = uconvert(m, Le)
    Sigma_SI = uconvert(Pa, Sigma)
    Time_SI = Eta / Sigma_SI
    t = uconvert(s, Time_SI)

    return GeoUnits{SI}(;
        length = Le,
        temperature = T,
        stress = Sigma,
        viscosity = Eta,
        time = t,
        m = Le_SI,
        K = T_SI,
        Pa = Sigma_SI,
        s = Time_SI,
    )
end

"""
    CharUnits = NO_units(length=1, temperature=1, stress=1, viscosity=1)

Specify the characteristic values in non-dimensional units

# Examples:
```julia-repl
julia> using GeoParamsUnits;
julia> CharUnits = NO_units()
Employing NONE units
Characteristic values:
         length:      1
         time:        1.0
         stress:      1
         temperature: 1.0
```
"""
function NO_units(; length = 1, temperature = 1, stress = 1, viscosity = 1)
    if unit(temperature) != NoUnits
        error("temperature should not have units")
    end
    if unit(length) != NoUnits
        error("length should not have units")
    end
    if unit(stress) != NoUnits
        error("stress should not have units")
    end
    if unit(viscosity) != NoUnits
        error("viscosity should not have units")
    end

    T = temperature
    Le = length
    Sigma = stress
    Eta = viscosity
    Time = Eta / Sigma

    return GeoUnits{NONE}(;
        length = Le,
        temperature = T,
        stress = Sigma,
        viscosity = Eta,
        time = Time,
        m = Le,
        K = T,
        Pa = Sigma,
        s = Time,
    )
end
