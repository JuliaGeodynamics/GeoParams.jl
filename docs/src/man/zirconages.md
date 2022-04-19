# Zircon age parameterizations

# Methods
Zircons are one of the ways in which we can date the age & activity of magmatic systems. 
Here, we provide a computational routine that computes the zircon age distribution from temperature-time paths

```@docs
GeoParams.ZirconAgeData
```
# Computational routines
There is one main routine with which you can compute zircon age probability density functions from a range of temperature-ime paths:

```@docs
GeoParams.ZirconAge.compute_zircon_age_PDF
```

This, in turn, calls two other routines:
```@docs
GeoParams.ZirconAge.compute_zircons_Ttpath
GeoParams.ZirconAge.zircon_age_PDF
```

We also provide a plotting routine, provided the `Plots` package is loaded, which produces figures such as:
![subet3](./assets/img/ZirconAge_PDF.png)
