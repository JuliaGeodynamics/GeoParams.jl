# Melting Parameterizations

# Methods
A number of melting parameterisations are implemented, which can be set with:

```@docs
GeoParams.MeltingParam.MeltingParam_Caricchi
GeoParams.MeltingParam.MeltingParam_Smooth3rdOrder
GeoParams.MeltingParam.MeltingParam_5thOrder
GeoParams.MeltingParam.MeltingParam_4thOrder
GeoParams.MeltingParam.MeltingParam_Quadratic
GeoParams.MeltingParam.MeltingParam_Assimilation
GeoParams.MeltingParam.SmoothMelting
```
# Computational routines
To compute the melt fraction at given `T` and `P`, use:
```@docs
GeoParams.MeltingParam.compute_meltfraction!
GeoParams.MeltingParam.compute_meltfraction
```

You can also obtain the derivative of melt fraction versus temperature with (useful to compute latent heat effects):
```@docs
GeoParams.MeltingParam.compute_dϕdT!
GeoParams.MeltingParam.compute_dϕdT
```

Also note that phase diagrams can be imported using `PerpleX_LaMEM_Diagram`, which may also have melt content information.
The computational routines work with that as well.

# Plotting routines
You can use the routine `PlotMeltFraction` to create a plot, provided that the `GLMakie` package has been loaded
```@docs
GeoParams.PlotMeltFraction
```