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

Devices without `Float64` support, such as Metal, cannot receive the stored `Float64`
parameters at all. Convert the laws, or the phases' `MaterialParams`, once on the host
before passing them to a kernel:

```julia
phases32 = convert_precision(Float32, phases)   # a tuple of MaterialParams
```

Loading Metal activates an extension that disables the `Float64` recomputation
described above inside Metal kernels, so a `Float32` result that leaves the range is
returned as `Inf` or `0` there, and prints a reminder of the workflow below.
Nondimensionalize the laws before converting them: nondimensional parameters and
intermediates stay within the `Float32` range, whereas SI creep-law prefactors such as
`1e-55` do not.

```julia
phases32 = convert_precision(Float32, nondimensionalize(phases, GEO_units()))
```

The seismic-velocity corrections (`melt_correction`, `melt_correction_Takei`,
`porosity_correction`, `anelastic_correction`) and `find_Xco2` run on the host only.

`GEO_units` and `SI_units` choose the characteristic amount of substance so that the
gas constant is 1. Activation energies then appear in units of `R` times the
characteristic temperature, and molar masses stay within the `Float32` range.

```@docs
precision_of
convert_precision
```
