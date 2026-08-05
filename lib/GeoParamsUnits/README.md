# GeoParamsUnits.jl

`GeoParamsUnits` is the Unitful-backed units and nondimensionalization package used by
[`GeoParams.jl`](https://github.com/JuliaGeodynamics/GeoParams.jl). It can also be used on its
own, without loading GeoParams' material-parameter implementations.

## Quick start

```julia
import GeoParamsUnits as U

scales = U.GEO_units(;
    length = 1000 * U.km,
    temperature = 1000 * U.C,
    stress = 10 * U.MPa,
    viscosity = 1.0e20 * U.Pas,
)

velocity = 3 * U.cm / U.yr
velocity_nd = U.nondimensionalize(velocity, scales)
velocity_again = U.dimensionalize(velocity_nd, U.cm / U.yr, scales)
```

`GEO_units`, `SI_units`, and `NO_units` construct characteristic scales. `GeoUnit` retains a
value's physical dimensions while it is nondimensional, allowing `dimensionalize` to restore it
later. Scalar, array, tuple, arithmetic, display, and `@unpack_val`/`@unpack_units` workflows are
supported.

## Relationship to GeoParams

This is GeoParams' default backend. Existing names such as `GeoParams.GeoUnit`,
`GeoParams.GEO_units`, and `GeoParams.Units` continue to refer to this package.

For a DynamicQuantities-backed implementation with the same exported API, see
[`GeoParamsDynamicUnits`](../GeoParamsDynamicUnits). Loading that package alongside GeoParams
activates GeoParams' conversion extension; GeoParams material structures continue to store this
package's `GeoUnit` internally for compatibility.

## Development

From the repository root:

```sh
julia --project=lib/GeoParamsUnits --startup-file=no -e 'using Pkg; Pkg.test()'
julia --project=lib/GeoParamsUnits --startup-file=no lib/GeoParamsUnits/benchmark/benchmarks.jl
```

The package is independently versioned under `lib/GeoParamsUnits`, requires Julia 1.10 or later,
and is tested independently in the repository CI.
