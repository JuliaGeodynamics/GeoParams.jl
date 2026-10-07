# Precision

```@meta
DocTestSetup = :(using GeoParams)
```

Material laws evaluate in the precision of the state they are given, not in that of
their stored parameters. The shipped databases store `Float64` values, but a
`Float32` stress returns a `Float32` strain rate, so a GPU kernel working in `Float32`
stays in `Float32` across a rheology call:

```jldoctest precision
julia> law = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003);

julia> compute_εII(law, 1.0e6; T = 1600.0) isa Float64
true

julia> ε32 = compute_εII(law, 1.0f6; T = 1600.0); ε32 isa Float32
true

julia> ε32 ≈ compute_εII(law, 1.0e6; T = 1600.0)
true
```

The state argument a law is a function of (the stress, strain rate, temperature or
pressure) selects the precision; the other arguments and the stored parameters are
converted to it. Integers express no preference, so `T = 1600` works at either
precision, and a dual number follows the precision of the value it carries. Stored
parameters are never modified.

In SI units some creep laws pass through intermediates outside the `Float32` range
(prefactors near `1e-55`, stresses raised to `n` near `1e39`) even when the result is
an ordinary `Float32`. Such laws recompute in `Float64` when their `Float32` result
overflows or underflows, and return the answer as `Float32`. Nondimensionalized laws
(see [Nondimensionalization](@ref)) stay within range and never take this path.

```@docs
precision_of
convert_precision
```
