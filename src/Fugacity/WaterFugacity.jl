module WaterFugacity

# Fugacity of pure H2O and the water content of olivine that follows from it.
#
# All routines take and return SI units (P in Pa, T in K, fugacity in Pa); each
# also has a method for Unitful quantities.

using Unitful
using ..Units: Pa, K

export compute_water_fugacity,
    compute_water_fugacity_Shinevar,
    compute_water_content_olivine

# Gas constant used by Pitzer & Sterner (1994), in Pa cm³ mol⁻¹ K⁻¹; densities in
# their equation of state are in mol cm⁻³, so ρ R T is a pressure in Pa.
const R_PS94 = 8.31451e6

# Gas constant [J mol⁻¹ K⁻¹] for the Arrhenius-type expressions below
const R_gas = 8.3145

# Coefficients c_ij of the H2O equation of state of Pitzer & Sterner (1994),
# J. Chem. Phys. 101, 3111-3116, Table 1 (with ρ in mol cm⁻³ and T in K). Row i
# holds the factors of T⁻⁴, T⁻², T⁻¹, 1, T and T² that make up c_i:
#   c_i = c_i1 T⁻⁴ + c_i2 T⁻² + c_i3 T⁻¹ + c_i4 + c_i5 T + c_i6 T²
# These are the values used by the online fugacity calculator of A. Withers
# (https://www.esci.umn.edu/people/researchers/withe012/fugacity.htm), which
# test/test_WaterFugacity.jl reproduces to better than 0.4%.
const PS94_COEFFICIENTS = (
    (0.0, 0.0, 0.24657688e6, 0.51359951e2, 0.0, 0.0),
    (0.0, 0.0, 0.58638965e0, -0.28646939e-2, 0.31375577e-4, 0.0),
    (0.0, 0.0, -0.6278384e1, 0.14791599e-1, 0.35779579e-3, 0.15432925e-7),
    (0.0, 0.0, 0.0, -0.42719875e0, -0.16325155e-4, 0.0),
    (0.0, 0.0, 0.56654978e4, -0.16580167e2, 0.76560762e-1, 0.0),
    (0.0, 0.0, 0.0, 0.10917883e0, 0.0, 0.0),
    (0.38878656e13, -0.13494878e9, 0.30916564e6, 0.75591105e1, 0.0, 0.0),
    (0.0, 0.0, -0.65537898e5, 0.18810675e3, 0.0, 0.0),
    (-0.14182435e14, 0.1816539e9, -0.19769068e6, -0.23530318e2, 0.0, 0.0),
    (0.0, 0.0, 0.92093375e5, 0.12246777e3, 0.0, 0.0),
)

# Temperature-dependent coefficients c_1 … c_10 at temperature T [K]
@inline function _ps94_c(T::F) where {F}
    iT = inv(T)
    iT2 = iT * iT
    return map(PS94_COEFFICIENTS) do c
        F(c[1]) * iT2 * iT2 + F(c[2]) * iT2 + F(c[3]) * iT + F(c[4]) + F(c[5]) * T + F(c[6]) * T * T
    end
end

# Denominator c_2 + c_3 ρ + c_4 ρ² + c_5 ρ³ + c_6 ρ⁴ that appears in the EOS
@inline _ps94_D(ρ, c) = c[2] + ρ * (c[3] + ρ * (c[4] + ρ * (c[5] + ρ * c[6])))

# Pressure [Pa] of H2O at molar density ρ [mol cm⁻³] and temperature T [K]
# (Pitzer & Sterner 1994, eq. 2)
@inline function _ps94_pressure(ρ, T, c)
    D = _ps94_D(ρ, c)
    ρ2 = ρ * ρ
    Z = ρ + c[1] * ρ2 -
        ρ2 * (c[3] + ρ * (2 * c[4] + ρ * (3 * c[5] + 4 * c[6] * ρ))) / (D * D) +
        c[7] * ρ2 * exp(-c[8] * ρ) + c[9] * ρ2 * exp(-c[10] * ρ)
    return oftype(T, R_PS94) * T * Z
end

# ln of the fugacity [Pa] at density ρ, from the residual Helmholtz energy
# A_res/RT (Pitzer & Sterner 1994, eq. 1):
#   ln f = ln(ρ R T) + A_res/RT + P/(ρ R T) - 1
@inline function _ps94_lnf(ρ, P, T, c)
    A_res = c[1] * ρ + inv(_ps94_D(ρ, c)) - inv(c[2]) -
        c[7] / c[8] * (exp(-c[8] * ρ) - 1) - c[9] / c[10] * (exp(-c[10] * ρ) - 1)
    ρRT = ρ * oftype(T, R_PS94) * T
    return log(ρRT) + A_res + P / ρRT - 1
end

# Root of g in the bracket [a, b]; plain bisection so that it also compiles inside GPU kernels
@inline function _bisect(g, a, b)
    ga = g(a)
    for _ in 1:200
        m = (a + b) / 2
        (m == a || m == b) && break
        gm = g(m)
        if sign(gm) == sign(ga)
            a, ga = m, gm
        else
            b = m
        end
    end
    return (a + b) / 2
end

# Molar density [mol cm⁻³] of the stable H2O phase at (P, T).
#
# Below the critical temperature the EOS can have three roots at a given pressure
# (vapour, an unstable branch and liquid). The density axis is therefore scanned on a
# logarithmic grid for the lowest- and highest-density sign changes; each is refined
# by bisection and the root with the lower fugacity (i.e. Gibbs energy) is returned.
# Above the critical temperature the two coincide.
function _ps94_density(P::F, T::F, c) where {F}
    g(ρ) = _ps94_pressure(ρ, T, c) - P
    nscan = 64
    ρmin, ρmax = F(1.0e-12), F(0.3)
    ratio = (ρmax / ρmin)^(one(F) / nscan)

    lo_a = lo_b = hi_a = hi_b = zero(F)
    found = false
    ρa, ga = zero(F), -P
    ρb = ρmin
    for _ in 0:nscan
        gb = g(ρb)
        if sign(ga) != sign(gb)
            if !found
                lo_a, lo_b = ρa, ρb
                found = true
            end
            hi_a, hi_b = ρa, ρb
        end
        ρa, ga = ρb, gb
        ρb *= ratio
    end
    found || return F(NaN)

    ρ_lo = _bisect(g, lo_a, lo_b)
    lo_a == hi_a && return ρ_lo
    ρ_hi = _bisect(g, hi_a, hi_b)
    return _ps94_lnf(ρ_lo, P, T, c) <= _ps94_lnf(ρ_hi, P, T, c) ? ρ_lo : ρ_hi
end

"""
    fH2O = compute_water_fugacity(P, T)

Fugacity of pure H2O [Pa] at pressure `P` [Pa] and temperature `T` [K], computed with the
equation of state of Pitzer & Sterner (1994). The density of H2O is found by solving
the equation of state with a bracketed root finder; the fugacity then follows from the
residual Helmholtz energy ``A_\\mathrm{res}``:
```math
    \\ln f_\\mathrm{H_2O} = \\ln(\\rho R T) + \\frac{A_\\mathrm{res}}{RT} + \\frac{P}{\\rho R T} - 1
```
The equation of state is calibrated for 0-10 GPa and 0-2000 °C. `P` and `T` can also be
given as `Unitful` quantities, in which case the fugacity is returned in Pa.
For non-positive pressure (e.g. near the free surface of a numerical model) the fugacity
is zero, its physical limit; a wet creep law then gives zero strain rate.

Reference: Pitzer, K.S. & Sterner, S.M. (1994), Equations of state valid continuously
from zero to extreme pressures for H2O and CO2, J. Chem. Phys. 101, 3111-3116.

# Example
```jldoctest
julia> round(compute_water_fugacity(1.0e9, 1073.15) / 1.0e6, digits = 1)   # MPa
1642.0

julia> round(MPa, compute_water_fugacity(1GPa, 800C), digits = 1)
1642.0 MPa
```
"""
function compute_water_fugacity(P::Real, T::Real)
    P, T = promote(float(P), float(T))
    P <= 0 && return zero(P)        # physical limit f -> P as P -> 0; a NaN input still gives NaN
    c = _ps94_c(T)
    ρ = _ps94_density(P, T, c)
    return exp(_ps94_lnf(ρ, P, T, c))
end

compute_water_fugacity(P::Quantity, T::Quantity) =
    compute_water_fugacity(ustrip(Pa, P), ustrip(K, T)) * Pa

"""
    fH2O = compute_water_fugacity_Shinevar(P, T)

Approximate fugacity of H2O [Pa] at pressure `P` [Pa] and temperature `T` [K] for crustal
conditions, following Shinevar, Behn & Hirth (2015):
```math
    f_\\mathrm{H_2O} = A_1 \\exp\\left(-\\frac{A_2 + P A_3}{RT}\\right)
```
with ``A_1 = 5.521`` GPa, ``A_2 = 31.28`` kJ/mol and ``A_3 = -2.009 \\times 10^{-5}`` m³/mol
(parameters as quoted in eq. 5.7 of Becker & Faccenna's geodynamics textbook). For
0.2-1 GPa and 600-1000 K this is within a factor 1.5 of [`compute_water_fugacity`](@ref),
the deviation grows for colder conditions at high pressure (factor ~2.3 at 1 GPa and
500 K). `P` and `T` can also be given as `Unitful` quantities.

Reference: Shinevar, W.J., Behn, M.D. & Hirth, G. (2015), Compositional dependence of lower
crustal viscosity, Geophys. Res. Lett. 42, 8333-8340.

# Example
```jldoctest
julia> round(compute_water_fugacity_Shinevar(1.0e9, 800.0) / 1.0e6, digits = 1)   # MPa
1026.6
```
"""
function compute_water_fugacity_Shinevar(P::Real, T::Real)
    P, T = promote(float(P), float(T))
    F = typeof(P)
    A1, A2, A3 = F(5.521e9), F(31.28e3), F(-2.009e-5)
    return A1 * exp(-(A2 + P * A3) / (F(R_gas) * T))
end

compute_water_fugacity_Shinevar(P::Quantity, T::Quantity) =
    compute_water_fugacity_Shinevar(ustrip(Pa, P), ustrip(K, T)) * Pa

"""
    C_OH = compute_water_content_olivine(fH2O, P, T)

Solubility of water (hydrogen) in olivine, ``C_\\mathrm{OH}`` in H/10⁶ Si, as function of water
fugacity `fH2O` [Pa], pressure `P` [Pa] and temperature `T` [K] after Zhao, Ginsberg &
Kohlstedt (2004):
```math
    C_\\mathrm{OH} = A_w \\exp\\left(-\\frac{E_w + P V_w}{RT}\\right) f_\\mathrm{H_2O}
```
with ``A_w = 26`` H/10⁶ Si/MPa, ``E_w = 40`` kJ/mol and ``V_w = 10 \\times 10^{-6}`` m³/mol
(as quoted in eq. 5.9 of Becker & Faccenna's geodynamics textbook); in this expression
``f_\\mathrm{H_2O}`` is in MPa, which is taken care of internally. The arguments can also be
given as `Unitful` quantities; the result is always a plain number in H/10⁶ Si.

For water-saturated conditions, combine it with [`compute_water_fugacity`](@ref).

Reference: Zhao, Y.-H., Ginsberg, S.B. & Kohlstedt, D.L. (2004), Solubility of hydrogen in
olivine: dependence on temperature and iron content, Contrib. Mineral. Petrol. 147, 155-161.

# Example
```jldoctest
julia> P, T = 2.0e9, 1573.15;

julia> round(compute_water_content_olivine(compute_water_fugacity(P, T), P, T))   # H/10⁶ Si
3018.0
```
"""
function compute_water_content_olivine(fH2O::Real, P::Real, T::Real)
    fH2O, P, T = promote(float(fH2O), float(P), float(T))
    F = typeof(P)
    Aw, Ew, Vw = F(26), F(40.0e3), F(10.0e-6)
    return Aw * exp(-(Ew + P * Vw) / (F(R_gas) * T)) * fH2O * F(1.0e-6)
end

compute_water_content_olivine(fH2O::Quantity, P::Quantity, T::Quantity) =
    compute_water_content_olivine(ustrip(Pa, fH2O), ustrip(Pa, P), ustrip(K, T))

end
