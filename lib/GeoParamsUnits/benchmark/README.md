# GeoParamsUnits benchmarks

Run the allocation and runtime microbenchmarks from the package directory:

```sh
julia --project=. --startup-file=no benchmark/benchmarks.jl
```

The script warms each operation and reports the median of nine samples as CSV. Array
results are normalized per complete array conversion. For comparisons, run the same
script with the same Julia version and hardware in each checkout.
