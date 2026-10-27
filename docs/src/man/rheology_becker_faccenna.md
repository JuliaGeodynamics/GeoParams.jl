# Rheology figures of Becker & Faccenna

This page reproduces figures 5.6, 5.7, 5.10, 5.12 and 5.13 of chapter 5 ("Rheology") of
Becker, T.W. & Faccenna, C., *Tectonic Geodynamics* (draft, 2026), with the creep laws of the GeoParams
database and [CairoMakie](https://docs.makie.org/stable/). The flow laws, geotherms and other parameter
choices are the authors' (taken from the book and from the script that made the printed figures), so that the
figures can be used as a starting point for your own computations.

Each figure is made by one function of the Makie extension of GeoParams. Their keyword arguments default to the
book's choices, and the creep laws are keywords too, so any law of the database (or your own) can be plotted
instead. The functions are available once a Makie backend is loaded:

```julia
using GeoParams, CairoMakie
```

## Stress and strain-rate conventions

Laboratory flow laws relate the differential stress ``\sigma_d = \sigma_1 - \sigma_3`` to the axial
strain rate of a compression test. GeoParams converts them to the second invariants of the deviatoric stress
and strain-rate tensors (`Apparatus = AxialCompression`), using ``\tau_{II} = \sigma_d/\sqrt{3}`` and
``\dot\varepsilon_{II} = \tfrac{\sqrt{3}}{2}\dot\varepsilon_{axial}``, so that
``\dot\varepsilon_{II} = A\, \tfrac{3^{(n+1)/2}}{2}\, \tau_{II}^n \ldots`` (eq. 5.12 of the book, see also
Gerya 2019). On this page the strain rate ``\dot\varepsilon`` is ``\dot\varepsilon_{II}`` and all stresses
are differential stresses, ``\sigma_d = \sqrt{3}\,\tau_{II}``.

The printed figures use two other conventions, so the curves below differ from the printed ones by a
constant factor:
- Fig. 5.12 uses the olivine laws of Table 5.2. Its viscosities are the same as in GeoParams, but its
  stress axis shows ``2\tau_{II}``. At a given strain rate the printed stresses are therefore
  ``2/\sqrt{3} \approx 1.15`` times higher than here, and the brittle-ductile transitions are 0.3–0.8 km deeper.
- Figs. 5.6, 5.10 and 5.13 use the stress scale ``B = 2^{(1+3/n)/2} A^{-1/n}`` given below eq. (5.5),
  i.e. ``\dot\varepsilon = A\,\sigma^n/2^{(n+3)/2}``. At the same stress the printed strain rate is
  ``\sqrt{3}\, 2^{(n+1)/2}`` times lower than here (3.5 for diffusion creep, 4.9 for the Peierls law
  of Mei et al. (2010), 6.9 for ``n = 3``, 8.2 for ``n = 3.5``). At the same strain rate the printed stress of a
  power law is ``(\sqrt{3}\,2^{(n+1)/2})^{1/n}`` times higher: 1.95 (quartz, ``n = 2.8``), 1.91
  (``n = 3``), 1.83 (olivine, ``n = 3.5``). The bands of Fig. 5.6 and the creep curves of Fig. 5.13 are shifted
  accordingly; in Fig. 5.10 the strain-rate contours move by 0.5–0.9 orders of magnitude and the
  diffusion/dislocation boundary by 0.15 log units in stress.

The Byerlee limits, temperatures, water fugacities and water contents are not affected.

## Differences to the book's parameters

Apart from the conventions above, the GeoParams database differs from the parameters of the book in a few
places. Where the book uses its own choice, the default arguments of the plotting functions override the database value.

| Law | Book (Table 5.1/5.2, script) | GeoParams database | On this page |
|---|---|---|---|
| Wet olivine, dislocation, Table 5.2 (Fig. 5.12) | ``C_\mathrm{OH}^r`` with ``r = 1``, ``C_\mathrm{OH} = 1000`` | `wet_olivine1_Hirth_2003`: ``r = 1.2`` (Hirth & Kohlstedt 2003) | `f = 1000^(1/1.2)`, which gives ``f^{1.2} = 1000`` |
| Table 5.2, ``A`` in SI units | rounded (``7.4\cdot10^{-15}``, ``6\cdot10^{-15}``) | computed from Hirth & Kohlstedt (2003) | database (stress differs by 0.1 %) |
| Dry olivine, diffusion | ``V = 6\cdot10^{-6}`` m³/mol | ``V = 10\cdot10^{-6}`` m³/mol (range 2–10) | ``V = 6\cdot10^{-6}`` |
| Dry olivine, dislocation (Figs. 5.10, 5.13) | ``V = 18\cdot10^{-6}`` m³/mol | ``V = 14\cdot10^{-6}`` m³/mol (Table 5.2: 14) | ``V = 18\cdot10^{-6}`` |
| Quartz, Rutter & Brodie (2004) | ``n = 3``, ``\log_{10}A = -4.9`` | ``n = 2.97``, ``\log_{10}A = -4.93`` | database (stresses in Fig. 5.6 differ by −1 to +17 %) |
| Quartz, Tokle et al. (2019); quartz, Jaoul et al. (1984) | ``\log_{10}A = -11.9586``; ``A = 9.0\cdot10^{-6}`` | ``-11.959``; ``10^{-5.045}`` | database (< 0.1 %) |
| Quartz, Lu & Jiang (2019) (not used here) | ``n = 3`` | ``n = 4`` (Lu & Jiang 2019) | – |

## Figure 5.6: flow stresses of crustal and mantle rocks

Dislocation creep of water-saturated quartz, feldspar (anorthite) and olivine at ``p = 0.5`` GPa, with the
water fugacity of Pitzer & Sterner (1994). Each band spans strain rates from ``10^{-15}`` (low stress) to
``10^{-13}`` s⁻¹ (high stress); quartz is drawn in lighter colors above 650 °C. The piezometer estimates of
the printed figure are not shown.

```@example rb
using GeoParams, CairoMakie
PlotFlowStressTemperature()
```

Other laws, pressures or strain rates are passed as keywords. Each law is a named tuple with a label, the law
and a color; the water variable `f` defaults to water saturation:

```@example rb
quartz = [
    (label = "Tokle 2019", law = SetDislocationCreep(GeoParams.Dislocation.wet_quartzite_Tokle_2019), color = :red),
    (label = "Lusk 2021, high p", law = SetDislocationCreep(GeoParams.Dislocation.highP_wet_quartzite_Lusk_2021), color = :green),
    (label = "Lu & Jiang 2019", law = SetDislocationCreep(GeoParams.Dislocation.wet_quartzite_Lu_2019), color = :purple),
]
PlotFlowStressTemperature(; laws = quartz, P = 1.0e9)
```

## Figure 5.7: water fugacity and water content of olivine

Left: water fugacity from the equation of state of Pitzer & Sterner (1994) (colors) and from the crustal
approximation of Shinevar et al. (2015), eq. (5.7) (contours every 250 MPa), with two linear geotherms
(surface heat flow 50 and 100 mW/m², ``k = 3`` W/m/K, ``\rho = 2750`` kg/m³) between 20 and 40 km depth.
Right: water content of olivine at water saturation, eq. (5.9), with half-space cooling geotherms
(``T_m = 1350`` °C, ``\kappa = 10^{-6}`` m²/s, ``p = \rho g z`` with ``\rho = 3000`` kg/m³; these values
reproduce the printed geotherms to within about 50 K). The printed caption places the crustal geotherms
between 10 and 40 km, but the printed lines start at about 20 km, as here.

```@example rb
PlotWaterFugacity()
```

## Figure 5.10: deformation mechanism maps of dry olivine

Four dashpots in series (eq. 5.29): diffusion creep, dislocation creep and grain-boundary sliding of
Hirth & Kohlstedt (2003) and the Peierls law of Mei et al. (2010). For a prescribed differential stress, the
strain rate of each mechanism is computed with `compute_εII`; colors show the dominant mechanism and contours
the ``\log_{10}`` of the total strain rate. a) ``T = 1250`` °C and ``p = 1`` GPa; b) ``d = 3`` mm and
``p = 1`` GPa, Peierls creep only for ``T \le 1250`` °C; c) ``d = 3`` mm along a 100 Myr half-space cooling
geotherm (``T_m = 1400`` °C) to which a gradient of 0.15 K/km is added (chosen to reproduce the
printed contours), ``p = \rho g z`` with ``\rho = 3300`` kg/m³, and Peierls creep only for ``T \le 1250`` °C. The dashed
line in c) is the stress at which the Peierls law of Goetze & Evans (1979) reaches ``10^{-15}`` s⁻¹, but not
less than 300 MPa, the lower stress limit for that law in the book's script; only the 300 MPa part is drawn in
the printed figure.

```@example rb
PlotDeformationMechanismMaps()
```

The mechanisms are a list of named tuples `(name, law, color)`, optionally with `Tmax` [K] above which the
mechanism is switched off; `law` can also be a function of temperature, as for the two GBS regimes. For example,
the diffusion creep law of Faul & Jackson (2007) instead of Hirth & Kohlstedt (2003), without GBS and Peierls creep:

```@example rb
mechanisms = [
    (name = "diffusion", law = SetDiffusionCreep(GeoParams.Diffusion.dry_olivine_Faul_2007), color = RGBf(0.85, 0.85, 0.83)),
    (name = "dislocation", law = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003), color = RGBf(0.7, 0.63, 0.48)),
]
PlotDeformationMechanismMaps(; mechanisms, yield_law = nothing)
```

## Figure 5.12: strength of oceanic lithosphere

Half-space cooling geotherms (``T_m = 1400`` °C) for plate ages of 25, 50 and 150 Myr, creep at
``\dot\varepsilon = 5\cdot10^{-15}`` s⁻¹ of "wet" (solid) and "dry" (dashed) olivine in dislocation creep (Table
5.2 of the book, Hirth & Kohlstedt 2003), limited by Byerlee's law for extension with ``\lambda = 0`` (heavy
blue line). The thin blue line is Byerlee's law for compression with ``\lambda = 0``. Pressure uses
``\rho = 2750`` kg/m³ above and ``\rho = 3300`` kg/m³ below 7 km depth. Temperatures are dotted (lower
axis).

```@example rb
fig, strength = PlotOceanicStrength()
for (key, value) in sort(collect(strength))
    println(key, ": ", round(value, digits = 1), " TN/m")
end
fig
```

The depth-integrated strengths (0–100 km) printed above can be compared with the 8 (11) and 43 (61) TN/m quoted
in the book for the wet (dry) cases at 25 and 150 Myr. The values here are 5–8 % lower. Plotting ``2\tau_{II}``
as in the book raises them to 7.6 (10.4) and 42.0 (58.2) TN/m; the remaining few percent are probably due to
the integration in the book (its depth range or resolution is not stated).

## Figure 5.13: strength of continental lithosphere

Compression with ``\lambda = 0.4`` and ``\dot\varepsilon = 10^{-15}`` s⁻¹, linear geotherms for surface heat
flows of 45, 65 and 85 mW/m² (``k = 3`` W/m/K, ``T \le 1350`` °C), a 35 km thick crust (``\rho = 2750``
kg/m³) and a dry olivine mantle (``\rho = 3300`` kg/m³). a) The whole crust is wet anorthite (Rybacki et al.
2006); b) the upper third of the crust (11.7 km) is wet quartz (Tokle et al. 2019). The crustal laws use the
water fugacity at saturation.

```@example rb
PlotContinentalStrength()
```

Each panel lists its crustal layers from top to bottom, with the depth of their base; below the last layer the
`mantle` law applies. For example, a granite upper crust over a mafic granulite lower crust:

```@example rb
granite = SetDislocationCreep(GeoParams.Dislocation.granite_Tirel_2008)
granulite = SetDislocationCreep(GeoParams.Dislocation.mafic_granulite_Wilks_1990)
panels = [(title = "granite / granulite", layers = [(bottom = 20.0e3, law = granite), (bottom = 35.0e3, law = granulite)])]
PlotContinentalStrength(; panels, size = (550, 650))
```

## Plotting functions

```@docs
PlotFlowStressTemperature
PlotWaterFugacity
PlotDeformationMechanismMaps
PlotOceanicStrength
PlotContinentalStrength
```
