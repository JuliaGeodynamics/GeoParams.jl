# Zircon saturation parameterizations

# Methods
Zircons are one of the ways in which we can date the age & activity of magmatic systems. 
Here, we implement simple parameterizations that estimates when zircons crystallize, which can be implemented in thermomechanical simulations.

```@docs
GeoParams.ZirconSaturation.Tierney
```
# Computational routines
To compute the zircon fraction `Fzrs` at given `T` and `P`, use:
```@docs
GeoParams.ZirconSaturation.compute_zirconsaturation!
GeoParams.ZirconSaturation.compute_zirconsaturation
```


