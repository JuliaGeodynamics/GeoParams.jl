# CreepLaws

The following viscous creep laws are implemented:
```@docs
GeoParams.LinearViscous
GeoParams.PowerlawViscous
GeoParams.DislocationCreep
GeoParams.MaterialParameters.ConstitutiveRelationships.ViscosityPartialMelt_Costa_etal_2009
GeoParams.MaterialParameters.ConstitutiveRelationships.LinearMeltViscosity
GeoParams.MaterialParameters.ConstitutiveRelationships.GiordanoMeltViscosity
GeoParams.MaterialParameters.ConstitutiveRelationships.DiffusionCreep
GeoParams.MaterialParameters.ConstitutiveRelationships.SetDiffusionCreep
GeoParams.MaterialParameters.ConstitutiveRelationships.GrainBoundarySliding
GeoParams.MaterialParameters.ConstitutiveRelationships.SetGrainBoundarySliding
```

## Computational routines for creep laws

Grain-size sensitive laws (`DiffusionCreep`, `GrainBoundarySliding`) store a grain size `d` (5 mm and 1 mm by default; change it with e.g. `SetDiffusionCreep(law; d = 100μm)`). It is used unless a grain size is passed at runtime, e.g. `compute_εII(law, τII; T, P, d = 1.0e-4)` or `compute_εII(law, τII, (; T, P, d))`, in which case the runtime value (a scalar, or an array in the in-place `!` routines) overrides the stored one.

Once a creep rheology is defined, we can use the following routines to perform computations within the solvers
```@docs
GeoParams.compute_εII
GeoParams.compute_εII!
GeoParams.compute_τII
GeoParams.compute_τII!
```

## Water fugacity and water content

Wet flow laws take the water term through the argument `f` (default `1`), as in `compute_εII(law, τII; T, P, f)`.
For laws written with the water fugacity ``f_\mathrm{H_2O}``, `f` is always passed in Pa (for example from
`compute_water_fugacity`); the database stores `A` in ``\mathrm{MPa}^{-(n+r)}`` (diffusion creep:
``\mathrm{MPa}^{-(n+r)}\,\mathrm{μm}^{-p}``) so that the conversion to SI also rescales the fugacity part.
Laws written with the water content ``C_\mathrm{OH}`` instead take a dimensionless `f = C_OH` in H/10⁶ Si; at present these
are `Dislocation.wet_olivine1_Hirth_2003` and `Diffusion.wet_olivine1_Hirth_2003` (Hirth & Kohlstedt 2003,
constant ``C_\mathrm{OH}``). Omitting `f` sets it to 1 (1 Pa or 1 H/10⁶ Si), which for a wet law is rarely intended.

## Refitted database entries

Database entries whose name ends in `_refit` hold a prefactor `A` refitted so that the law reproduces the paper's own
data where the printed equation does not; the printed law is kept under the original name, and both comments give the
evidence. At present these are `Dislocation.dry_olivine_Hirth_2003_refit` (A = 2.8e5 instead of 1.1e5
``\mathrm{MPa}^{-3.5}\,\mathrm{s}^{-1}``), `Dislocation.wet_olivine_Karato_2003_refit` (10^2.7 instead of 10^2.9
``\mathrm{MPa}^{-4.2}\,\mathrm{s}^{-1}``) and `Dislocation.wet_quartzite_Lu_2019_refit` (1.0e-15 instead of 6.0e-15
``\mathrm{MPa}^{-6.7}\,\mathrm{s}^{-1}``).
