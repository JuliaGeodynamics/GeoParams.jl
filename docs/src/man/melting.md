# Melting Parameterizations

# Methods
A number of melting parameterisations are implemented, which can be set with:

```@docs
GeoParams.MeltingParam.MeltingParam_Caricchi
GeoParams.MeltingParam.MeltingParam_5thOrder
GeoParams.MeltingParam.MeltingParam_4thOrder
GeoParams.MeltingParam.MeltingParam_Quadratic
```
# Computational routines
To compute the melt fraction at given `T` and `P`, use:
```@docs
GeoParams.MeltingParam.compute_meltfraction!
GeoParams.MeltingParam.compute_meltfraction
```

Also note that phase diagrams can be imported using `PerpleX_LaMEM_Diagram`.

# Plotting routines
You can use the routine `PlotMeltFraction` to create a plot, provided that the `Plots` package has been loaded
```@docs
GeoParams.Plotting.PlotMeltFraction
```