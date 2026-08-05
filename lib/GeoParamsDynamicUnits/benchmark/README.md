# GeoParamsDynamicUnits benchmarks

From the package directory, run:

```sh
julia --project=. --startup-file=no benchmark/benchmarks.jl
```

The dependency-free harness reports warmed median runtime and allocation results as CSV.
