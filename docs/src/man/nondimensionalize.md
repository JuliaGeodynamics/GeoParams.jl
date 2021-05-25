# Nondimensionalization

Create a nondimensionalization object in which we specify characteristic values, which are later employed to non-dimensionalize (or dimensionalize) all model parameters. Choose between `GEO`, `SI` or `NO` units:
- `SI` units assume all input and output is in `SI` units. Very general, but for typical geodynamic simulations often not so useful (as a million years has many seconds, resulting in large numbers).
- `GEO` units uses `SI` units throughout but assumes that input/output are in a format that is more convenient for typical geodynamic use cases, such as `Myrs`,`km` or `MPa`
- `NO` units are nondimensional units. Note that for parameters to be correctly non-dimensionalized in this case, you still have to indicate units (such as that `velocity` is given in `m/s`).

A dimensional parameter can be transformed into a non-dimensional one with `Nondimensionalize`.

```@index
```

# Specify characteristic values
Characteristic values can be defined in 3 ways.

```@docs
GEO_units
SI_units
NO_units
```

# (Non)-dimensionalize parameters
Once characteristic values have been defined, you can use them to non-dimensionalize or dimensionalize any parameter.
```@docs
Nondimensionalize
```