# Unit backends

GeoParams ships two independently usable unit packages. They expose the same GeoParams-specific
API but use different underlying quantity representations.

| Package | Quantity backend | Intended use |
|:--|:--|:--|
| `GeoParamsUnits` | Unitful | Default GeoParams behavior and existing applications |
| `GeoParamsDynamicUnits` | DynamicQuantities | Runtime dimensions and DynamicQuantities-based applications |

Both packages provide the exported unit constants, `GeoUnit`, `GeoUnits`, `GEO_units`,
`SI_units`, `NO_units`, `nondimensionalize`, `dimensionalize`, `compute_units`, and unpacking
helpers. Matching API and behavior do not imply matching concrete types: quantities and unit
metadata belong to the selected backend.

## Standalone use

Use qualified imports when comparing or loading both packages, since they intentionally export
the same names.

```julia
import GeoParamsUnits as UU
import GeoParamsDynamicUnits as DU

unitful_scales = UU.GEO_units()
dynamic_scales = DU.GEO_units()

UU.nondimensionalize(3 * UU.cm / UU.yr, unitful_scales)
DU.nondimensionalize(3 * DU.cm / DU.yr, dynamic_scales)
```

## DynamicQuantities interoperability

`GeoParamsUnits` remains GeoParams' default internal backend. Loading
`GeoParamsDynamicUnits` activates a package extension that accepts DynamicQuantities quantities,
`GeoUnit` objects, and characteristic scales at the GeoParams boundary.

```julia
using GeoParams
import GeoParamsDynamicUnits as DU

density = ConstantDensity(; ρ = 2900 * DU.kg / DU.m^3)
scales = DU.GEO_units()

density_nd = nondimensionalize(density, scales)
density_dim = dimensionalize(density_nd, scales)
```

The extension also supplies explicit conversions:

```julia
unitful_scales = convert(GeoParams.Units.GeoUnits, scales)
dynamic_scales = convert(DU.GeoUnits, unitful_scales)

dynamic_parameter = DU.GeoUnit(3 * DU.cm / DU.yr)
unitful_parameter = convert(GeoParams.Units.GeoUnit, dynamic_parameter)
```

Conversion covers SI base dimensions, derived and fractional dimensions, arrays, affine
temperature, and the custom geoscience units provided by the packages.

## Loading and storage behavior

Both sibling packages are installed with GeoParams so local `lib/` dependencies resolve on the
minimum supported Julia version. The dynamic backend is not loaded by `using GeoParams`; the
extension activates only after `GeoParamsDynamicUnits` is explicitly loaded.

GeoParams material structures continue to store `GeoParamsUnits.GeoUnit`. Keeping this stable
preserves the existing public types and serialized representation. The extension therefore
provides backend interoperability rather than a global runtime switch for internal storage.

See [Nondimensionalization](@ref) for characteristic-scale workflows and [`GeoUnit`](@ref) for
the stored dimensional metadata model.
