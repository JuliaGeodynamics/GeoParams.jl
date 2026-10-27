# Water fugacity

Flow laws of hydrous rocks and the solubility of water in nominally anhydrous minerals depend on the
fugacity of water, ``f_\mathrm{H_2O}``. GeoParams provides the following routines; all take pressure in Pa
and temperature in K and return the fugacity in Pa, or accept `Unitful` quantities.

## Pitzer & Sterner (1994)

[`compute_water_fugacity`](@ref) solves the equation of state of Pitzer & Sterner (1994) for the molar
density ``\rho`` of H2O at given ``P`` and ``T``,
```math
    \frac{P}{RT} = \rho + c_1\rho^2 - \rho^2\frac{c_3 + 2c_4\rho + 3c_5\rho^2 + 4c_6\rho^3}{(c_2 + c_3\rho + c_4\rho^2 + c_5\rho^3 + c_6\rho^4)^2} + c_7\rho^2 e^{-c_8\rho} + c_9\rho^2 e^{-c_{10}\rho},
```
where each coefficient is a function of temperature, ``c_i = c_{i1}T^{-4} + c_{i2}T^{-2} + c_{i3}T^{-1} + c_{i4} + c_{i5}T + c_{i6}T^2``.
The fugacity follows from the residual Helmholtz energy ``A_\mathrm{res}``:
```math
    \ln f_\mathrm{H_2O} = \ln(\rho R T) + \frac{A_\mathrm{res}}{RT} + \frac{P}{\rho R T} - 1,
```
```math
    \frac{A_\mathrm{res}}{RT} = c_1\rho + \frac{1}{c_2 + c_3\rho + c_4\rho^2 + c_5\rho^3 + c_6\rho^4} - \frac{1}{c_2} - \frac{c_7}{c_8}\left(e^{-c_8\rho} - 1\right) - \frac{c_9}{c_{10}}\left(e^{-c_{10}\rho} - 1\right).
```
The implementation reproduces the online calculator of A. Withers to better than 0.4% for 0.1-2.8 GPa and
100-1500 °C.

## Crustal approximation

[`compute_water_fugacity_Shinevar`](@ref) uses the closed-form approximation of Shinevar et al. (2015),
```math
    f_\mathrm{H_2O} \approx A_1 \exp\left(-\frac{A_2 + P A_3}{RT}\right),
```
which is within a factor 1.5 of the Pitzer & Sterner value for 0.2-1 GPa and 600-1000 K.

## Water content of olivine

[`compute_water_content_olivine`](@ref) gives the hydrogen solubility in olivine of Zhao et al. (2004),
```math
    C_\mathrm{OH} = A_w \exp\left(-\frac{E_w + P V_w}{RT}\right) f_\mathrm{H_2O},
```
in H/10⁶ Si (with ``f_\mathrm{H_2O}`` in MPa in this expression; the routine takes it in Pa).

## Example

```julia
julia> using GeoParams

julia> fH2O = compute_water_fugacity(1GPa, 800C);    # Pitzer & Sterner (1994), in Pa

julia> round(MPa, fH2O, digits = 1)
1642.0 MPa

julia> round(MPa, compute_water_fugacity_Shinevar(1GPa, 800C), digits = 1)   # crustal approximation
1575.3 MPa

julia> round(compute_water_content_olivine(fH2O, 1GPa, 800C))   # H/10⁶ Si, water-saturated olivine
157.0
```

## Routines

```@docs
GeoParams.WaterFugacity.compute_water_fugacity
GeoParams.WaterFugacity.compute_water_fugacity_Shinevar
GeoParams.WaterFugacity.compute_water_content_olivine
```
