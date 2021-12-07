# Phase Diagrams

It is possible to employ phase diagrams as lookup tables, which can be used to compute density, melt fraction or seismic properties, for example. 


### Perple_X 
A popular way to compute phase diagrams (pseudosections) for a given (fixed) chemical composition as a function of pressure and temperature is by using [Perple_X](https://www.perplex.ethz.ch). Once a stable assemblage is computed, through Gibbs energy minimisation, it allows you to output a large number of properties such as density, melt fraction or seismic velocities.

Within GeoParams, we can import this diagram (provided it is formatted in the same manner as `LaMEM` expects `Perple_X` input to be).

```@docs
GeoParams.MaterialParameters.PhaseDiagrams.PhaseDiagram_LookupTable
GeoParams.MaterialParameters.PhaseDiagrams.Read_LaMEM_Perple_X_Diagram
```

