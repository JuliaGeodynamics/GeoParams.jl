using GeoParamsDynamicUnits
using DynamicQuantities

median_value(values) = sort(values)[(length(values) + 1) ÷ 2]

function measure(f, operations; samples = 9)
    f()
    times = Float64[]
    bytes = Float64[]
    for _ in 1:samples
        GC.gc()
        result = @timed f()
        push!(times, result.time * 1.0e9 / operations)
        push!(bytes, result.bytes / operations)
    end
    return median_value(times), median_value(bytes)
end

function property_access(units, iterations)
    total = 0.0
    for _ in 1:iterations
        total += ustrip(units.length)
    end
    return total
end

function scalar_nondimensionalize(values, units, iterations)
    total = 0.0
    mask = length(values) - 1
    for i in 0:(iterations - 1)
        total += nondimensionalize(values[(i & mask) + 1], units)
    end
    return total
end

function array_nondimensionalize(values, units, iterations)
    result = nondimensionalize(values, units)
    for _ in 2:iterations
        result = nondimensionalize(values, units)
    end
    return result
end

function construct_units(lengths, iterations)
    results = Vector{Any}(undef, iterations)
    mask = length(lengths) - 1
    for i in 0:(iterations - 1)
        results[i + 1] = GEO_units(; length = lengths[(i & mask) + 1])
    end
    return results
end

function report(name, f, operations)
    nanoseconds, allocated_bytes = measure(f, operations)
    return println(name, ',', nanoseconds, ',', allocated_bytes)
end

units = GEO_units()
scalar_values = [(3.0 + i / 1024) * cm / yr for i in 1:1024]
small_array = collect(1.0:16.0) .* m
large_array = collect(1.0:10_000.0) .* m
lengths = [(1000.0 + i / 1024) * km for i in 1:1024]

println("julia_version,", VERSION)
println("case,ns_per_operation,bytes_per_operation")
report("property_access", () -> property_access(units, 1_000_000), 1_000_000)
report(
    "scalar_nondimensionalize",
    () -> scalar_nondimensionalize(scalar_values, units, 100_000),
    100_000,
)
report("array_16", () -> array_nondimensionalize(small_array, units, 2_000), 2_000)
report("array_10000", () -> array_nondimensionalize(large_array, units, 20), 20)
report("GEO_units_constructor", () -> construct_units(lengths, 1_000), 1_000)
