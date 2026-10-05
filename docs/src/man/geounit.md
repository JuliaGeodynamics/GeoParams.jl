```@meta
DocTestSetup = quote
    using GeoParams
end
```

# `GeoUnit`

Most parameters have physical units (100km, 9.81 ms⁻²). Yet, within numerical solvers it is usually not a good idea to compute with the actual SI units as they (at least in geosciences) tend to result in very large/small numbers, which can cause in roundoff errors. 
Therefore, we usually convert the dimensional units into non-dimensional ones before starting computations. The results can than be converted back to more convenient units for us when we generate the output (plots, for example).

This requires a way to store the dimensions of a parameter, but also, once they are made nondimensional, the original units (otherwise we can't convert them back to dimensional values).

For that reason, we use the `GeoUnit` structure, defined as 
```julia
struct GeoUnit{T,U}
    val          :: T
    unit         :: U
    isdimensional:: Bool
end
```

You can define a GeoUnit as:
```jldoctest geounit
julia> a = GeoUnit(100km)
GeoUnit{dimensional, km}, 
100.0

julia> b = GeoUnit(100s)
GeoUnit{dimensional, s}, 
100.0
```
Dividing the two gives the expected result
```jldoctest geounit
julia> c = a/b
GeoUnit{dimensional, m s⁻¹·⁰}, 
1000.0
```
Can also be a vector:
```jldoctest geounit
julia> T = GeoUnit(100K:10.0K:1000K)
GeoUnit{dimensional, K}, 
100.0:10.0:1000.0
```
The numerical value is obtained with:
```jldoctest geounit
julia> NumValue(c)
1000.0
```
It's unit with
```jldoctest geounit
julia> Unit(c)
m s⁻¹·⁰
```
and the value with unit with
```jldoctest geounit
julia> Value(c)
1000.0 m s⁻¹·⁰
```

If you multiply or divide a `GeoUnit` with a normal number, we assume that you are no longer interested in the units and you'll get a Float64:
```jldoctest geounit
julia> c*100
100000.0
```

The value of the `GeoUnit` struct becomes clear once we non-dimensionalize it (as explained in more detail in a separate page).
For that, first specify which characteristic values you want to use for non-dimensionalisation: 
```jldoctest geounit
julia> CharDim = GEO_units(length=1000km, temperature=1000C, stress=10MPa, viscosity=1e20Pas);
```

Let's define a stress:
```jldoctest geounit
julia> σ = GeoUnit(10MPa)
GeoUnit{dimensional, MPa}, 
10.0
```
And transfer it to nondimensional units:
```jldoctest geounit
julia> σ_nd = nondimensionalize(σ,CharDim)
GeoUnit{nondimensional, MPa}, 
0.9999999999999998
```
Instead of 1e7 Pa, we now have 1.0; yet at the same time we know that the original unit of σ_nd was in MPa, which allows to transfer it back to dimensional units:
```jldoctest geounit
julia> dimensionalize(σ_nd,CharDim)
GeoUnit{dimensional, MPa}, 
10.0
```
Some times, for e.g. plotting, we need to dimensionalize some variable and strip it from its units. This can be done with `ustrip`
```jldoctest geounit
julia> ustrip(dimensionalize(σ_nd, MPa, CharDim))
10.0
```
This is a bit verbose, so we provide a convenience function `dimensionalize_and_strip` that does this in one go:
```jldoctest geounit
julia> dimensionalize_and_strip(σ_nd, MPa, CharDim)
10.0
```
or even better, a macro `@dimstrip` that does the same:
```jldoctest geounit
julia> @dimstrip σ_nd MPa CharDim
10.0
```
