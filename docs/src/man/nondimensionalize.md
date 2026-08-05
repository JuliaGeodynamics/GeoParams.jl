# Nondimensionalization

Create a nondimensionalization object in which we specify characteristic values, which are later employed to non-dimensionalize (or dimensionalize) all model parameters. Choose between `GEO`, `SI` or `NO` units:
- `SI` units assume all input and output is in `SI` units. Very general, but for typical geodynamic simulations often not so useful (as a million years has many seconds, resulting in large numbers).
- `GEO` units uses `SI` units throughout but assumes that input/output are in a format that is more convenient for typical geodynamic use cases, such as `Myr`,`km` or `MPa`
- `NO` units are nondimensional units. Note that for parameters to be correctly non-dimensionalized in this case, you still have to indicate units (such as that `velocity` is given in `m/s`).

A dimensional parameter can be transformed into a non-dimensional one with `nondimensionalize`.

The units API is also available independently through `GeoParamsUnits` and
`GeoParamsDynamicUnits`. Both packages are installed with GeoParams, while `GeoParamsUnits`
remains the default backend, preserving `GeoParams.Units` and the existing Unitful-based API.
Loading `GeoParamsDynamicUnits` activates an extension that accepts DynamicQuantities
quantities, `GeoUnit` values, and `GeoUnits` characteristic scales at the GeoParams boundary.

```julia
using GeoParams
import GeoParamsDynamicUnits as DUnits

density = ConstantDensity(; ρ = 2900 * DUnits.kg / DUnits.m^3)
scales = DUnits.GEO_units()
density_nd = nondimensionalize(density, scales)
```

Material-parameter storage continues to use `GeoParamsUnits.GeoUnit` internally. This keeps the
existing public API and serialized material types stable while allowing DynamicQuantities input
and characteristic units through the extension. Explicit `convert` methods are provided between
the two packages' `GeoUnit` and `GeoUnits` types.

See [Unit backends](@ref) for standalone loading, backend comparison, conversion examples, and
the extension's storage contract.

# Specify characteristic values
Characteristic values can be defined in 3 ways.

```@docs
AbstractGeoUnit
GeoUnit
GeoUnits
GEO_units
SI_units
NO_units
@unpack_val
@unpack_units
```

# (Non)-dimensionalize parameters
Once characteristic values have been defined, you can use them to non-dimensionalize or dimensionalize any parameter.
```@docs
nondimensionalize
dimensionalize
isDimensional
```
