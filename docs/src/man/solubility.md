# Solubility

## Methods
Coupled H2O–CO2 volatile solubility (dissolved-content) closures can be set with:
```@docs
GeoParams.MaterialParameters.Solubility.Liu2005_Solubility
GeoParams.MaterialParameters.Solubility.Mafic_Solubility
GeoParams.MaterialParameters.Solubility.GasMixture
```
## Computational routines
To evaluate dissolved H2O and CO2 (mass fractions) within a user routine, use this:
```@docs
GeoParams.MaterialParameters.Solubility.compute_dissolved
GeoParams.MaterialParameters.Solubility.compute_dissolved!
```
Partial derivatives (via `ForwardDiff`, useful for implicit/Newton solves) are available for pressure, temperature, and gas composition:
```@docs
GeoParams.MaterialParameters.Solubility.∂dissolved_∂P
GeoParams.MaterialParameters.Solubility.∂dissolved_∂T
GeoParams.MaterialParameters.Solubility.∂dissolved_∂Xco2
```
Given a target dissolved H2O content, the gas composition `X_co2` that produces it can be found with:
```@docs
GeoParams.MaterialParameters.Solubility.find_Xco2
```
Gas-mixture helpers for the energy equation:
```@docs
GeoParams.MaterialParameters.Solubility.compute_gas_heatcapacity
GeoParams.MaterialParameters.Solubility.effective_molar_mass
```

Note that `ε_g` (exsolved-gas volume fraction) and `X_co2` (CO2 mole fraction of the gas) are typically ODE state variables owned by the solver, not closed algebraically at a point — the closures above are the *entries* to that mass-balance, not a replacement for it.
