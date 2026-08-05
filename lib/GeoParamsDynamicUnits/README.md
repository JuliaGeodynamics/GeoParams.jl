# GeoParamsDynamicUnits.jl

`GeoParamsDynamicUnits` implements the GeoParams units and nondimensionalization API using
[`DynamicQuantities.jl`](https://github.com/SymbolicML/DynamicQuantities.jl). Its exported names
and numerical behavior match `GeoParamsUnits`, while its concrete quantity and dimension types
come from DynamicQuantities.

## Quick start

```julia
import GeoParamsDynamicUnits as U

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

The package supports `GEO_units`, `SI_units`, `NO_units`, `GeoUnit`, scalar and array
transformations, tuples, arithmetic, fractional dimensions, affine Celsius, custom geoscience
units, display, and `@unpack_val`/`@unpack_units`.

## Using it with GeoParams

Loading this package activates `GeoParamsDynamicUnitsExt`:

```julia
using GeoParams
import GeoParamsDynamicUnits as DU

density = ConstantDensity(; ρ = 2900 * DU.kg / DU.m^3)
scales = DU.GEO_units()
density_nd = nondimensionalize(density, scales)

# Explicit conversion is also available after the extension loads.
unitful_scales = convert(GeoParams.Units.GeoUnits, scales)
```

Dynamic quantities, `DU.GeoUnit`, and `DU.GeoUnits` can cross the GeoParams API boundary. To
preserve existing public and serialized material types, GeoParams converts them to its default
`GeoParamsUnits.GeoUnit` storage. This is interoperability, not a runtime replacement of every
internal material field with a DynamicQuantities type.

`using GeoParams` alone does not load this package, DynamicQuantities, or the extension.

## Development

From the repository root:

```sh
julia --project=lib/GeoParamsDynamicUnits --startup-file=no -e 'using Pkg; Pkg.test()'
julia --project=lib/GeoParamsDynamicUnits --startup-file=no lib/GeoParamsDynamicUnits/benchmark/benchmarks.jl
```

The package is independently versioned under `lib/GeoParamsDynamicUnits`, requires Julia 1.10
or later, and is tested both independently and through the GeoParams extension in CI.
