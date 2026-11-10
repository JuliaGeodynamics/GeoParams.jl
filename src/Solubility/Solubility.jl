module Solubility

# Coupled H2O–CO2 volatile solubility (dissolved-content) closures for
# Degruyter & Huber (2014)-type magma-chamber models. Each closure returns
# *both* dissolved H2O and dissolved CO2 mass fractions at a point, because the
# two share the partial pressures Pw = P(1-X_co2), Pc = P·X_co2.

using Parameters, Unitful, LaTeXStrings, MuladdMacro
using ForwardDiff, StaticArrays
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
using GeoParams: fastpow, pow_check, @pow
import ..Units: isdimensional
using ..MaterialParameters: No_MaterialParam, MaterialParamsInfo
using ..Density: ThreePhase_Density, compute_density
import Base.show, GeoParams.param_info

include("../Computations.jl")

"""
    AbstractSolubility{T} <: AbstractMaterialParam

Supertype of the coupled H2O-CO2 solubility closures. A subtype gives the
dissolved H2O and CO2 mass fractions of a melt as a function of pressure,
temperature and gas composition; see [`compute_dissolved`](@ref).

Implemented closures are [`Liu2005_Solubility`](@ref) for silicic melts,
[`Mafic_Solubility`](@ref) for mafic melts and [`SunYao2026_Solubility`](@ref)
for melts of given oxide composition.
"""
abstract type AbstractSolubility{T} <: AbstractMaterialParam end

export compute_dissolved, # (m_h2o, m_co2) mass fractions
    compute_dissolved!, # in-place, fills two arrays
    compute_dissolved_ratio, # phase-ratio mix of both outputs
    ∂dissolved_∂P, # ForwardDiff partials, both outputs
    ∂dissolved_∂T,
    ∂dissolved_∂Xco2,
    find_Xco2, # Newton-invert compute_dissolved for X_co2 given target m_h2o
    compute_saturation_pressure, # Newton-invert compute_dissolved for (P, X_co2) given (m_h2o, m_co2)
    param_info,
    AbstractSolubility,
    Liu2005_Solubility, # silicic (rhyolite)
    Mafic_Solubility, # mafic (basalt), CO2 block still Liu rhyolite
    SunYao2026_Solubility, # any melt composition, from oxide wt%
    GasMixture, # H2O-CO2 gas-mixture properties
    compute_gas_heatcapacity, # mass-weighted specific heat c_g(X_co2)
    effective_molar_mass, # m_g(X_co2)
    compute_gas_fraction            # equilibrium exsolved H2O gas volume fraction

# Empty routine when no solubility is defined
compute_dissolved(::No_MaterialParam{_T}; kwargs...) where {_T} = (zero(_T), zero(_T))

# Liu et al. (2005) — silicic ------------------------------------------------
"""
    Liu2005_Solubility(; coeffs=(b1..b6, c1..c4), Pref=1MPa, Tref=1K)

Coupled H2O–CO2 solubility for silicic (rhyolite) melt after Liu et al. (2005),
as used by Degruyter & Huber (2014). It applies to rhyolite only and is useful to
reproduce Degruyter & Huber-type models; for other melt compositions, for CO2 and for
pressures above ~500 MPa use [`SunYao2026_Solubility`](@ref). [`compute_dissolved`](@ref) returns the
dissolved H2O and CO2 mass fractions of the melt as a function of pressure,
temperature, and the CO2 mole fraction of the gas `X_co2`:
```math
    m_{H_2O} = (b_1 P_w^{1/2} + b_2 P_w + b_3 P_w^{3/2})\\frac{T_{ref}}{T}
             + b_4 P_w^{3/2} + P_c(b_5 P_w^{1/2} + b_6 P_w)
```
```math
    m_{CO_2} = P_c(c_1 + c_2 P_w)\\frac{T_{ref}}{T} + P_c(c_3 P_w^{1/2} + c_4 P_w^{3/2})
```
with the dimensionless partial pressures ``P_w = P(1-X_{co2})/P_{ref}`` and
``P_c = P X_{co2}/P_{ref}``. With `Pref=1MPa`, `Tref=1K` these equal the reference
partial pressures in MPa and reproduce the Liu (2005) numbers; because they are
ratios of like-dimensioned quantities, the closure is dimensionally homogeneous
and nondimensionalizes cleanly. The fitted `coeffs` are dimensionless; `isdimensional`
tracks the reference `GeoUnit`s.

# References
- Liu, Y., Zhang, Y., Behrens, H. (2005), Solubility of H2O in rhyolitic melts at low pressures and a new empirical model for mixed H2O-CO2 solubility, JVGR 143, 219-235, https://doi.org/10.1016/j.jvolgeores.2004.09.019
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
"""
struct Liu2005_Solubility{T, U1, U2} <: AbstractSolubility{T}
    coeffs::NTuple{10, T}     # b1..b6 (H2O), c1..c4 (CO2); dimensionless
    Pref::GeoUnit{T, U1}      # pressure scale of the fit (1 MPa)
    Tref::GeoUnit{T, U2}      # temperature scale of the fit (1 K)

    function Liu2005_Solubility(;
            coeffs = (
                354.94, 9.623, -1.5223, 1.2439e-3, -1.084e-4, -1.362e-5,
                5668.0, -55.99, 0.4133, 2.041e-3,
            ),
            Pref = 1.0e6Pa, Tref = 1.0K
        )
        PrefU = convert(GeoUnit, Pref)
        TrefU = convert(GeoUnit, Tref)
        T = typeof(PrefU).types[1]
        U1 = typeof(PrefU).types[2]
        U2 = typeof(TrefU).types[2]
        return new{T, U1, U2}(T.(coeffs), PrefU, TrefU)
    end
    Liu2005_Solubility(coeffs, Pref, Tref) = Liu2005_Solubility(; coeffs, Pref, Tref)
end
isdimensional(s::Liu2005_Solubility) = isdimensional(s.Pref)

function param_info(s::Liu2005_Solubility)
    return MaterialParamsInfo(; Equation = L"m_{H_2O},\ m_{CO_2} = f_{Liu2005}(P, T, X_{co2})")
end

"""
    compute_dissolved(s::AbstractSolubility, P, T, X_co2) -> (m_h2o, m_co2)

Dissolved H2O and CO2 mass fractions of the melt (`P` in Pa, `T` in K, `X_co2`
the CO2 mole fraction of the gas). Also callable as `compute_dissolved(s; P, T,
X_co2)` and `compute_dissolved(s, args::NamedTuple)`.
"""
@inline function compute_dissolved(s::Liu2005_Solubility, P, T, X_co2)
    Tc = precision_of(P)
    b1, b2, b3, b4, b5, b6, c1, c2, c3, c4 = convert_precision(Tc, s.coeffs)
    @unpack_like P Tc Pref, Tref = s
    Pw = P * (1 - X_co2) / Pref     # ∝ H2O partial pressure in MPa
    Pc = P * X_co2 / Pref           # ∝ CO2 partial pressure in MPa
    Tr = Tref / T                   # ∝ 1/T [K]

    Pw_sqrt = sqrt(Pw)     # Pw^0.5
    Pw_15 = Pw * Pw_sqrt     # Pw^1.5

    meq = @muladd (b1 * Pw_sqrt + b2 * Pw + b3 * Pw_15) * Tr + b4 * Pw_15 + Pc * (b5 * Pw_sqrt + b6 * Pw)
    Cco2 = @muladd Pc * ((c1 + c2 * Pw) * Tr + c3 * Pw_sqrt + c4 * Pw_15)

    return (meq / 100, Cco2 / 1000000)   # wt% -> fraction ; ppm -> fraction
end

function show(io::IO, s::Liu2005_Solubility)
    return print(io, "Liu et al. (2005) silicic H2O-CO2 solubility")
end

# Mafic ----------------------------------------------------------------------
"""
    Mafic_Solubility(; coeffs=(b1..b10, c1..c4), T0=273.15K, Tref=1K, Pref=1MPa)

Coupled H2O–CO2 solubility for mafic (basalt) melt. Dissolved H2O follows the
mafic polynomial of Scholz et al. (2023; MagmaChamberBoxModel, `exsolve_mafic.m`) in the dimensionless
groups ``T_C = (T-T_0)/T_{ref}`` (numerically °C) and ``P_m = P/P_{ref}``
(numerically MPa); dissolved CO2 reuses the Liu (2005) rhyolite CO2 block.
Nondimensionalizes like [`Liu2005_Solubility`](@ref). The H2O polynomial has
no floor at zero and goes negative outside its calibration (high `X_co2`,
low `P`); the dissolved H2O output is floored at zero rather than returned
as a negative mass fraction. At low pressure it does not tend to zero: with pure
H2O gas it gives 2.74 wt% at `P → 0` and 850 °C (1.52 wt% at 1100 °C), so the
law is meaningful only at magma-chamber pressures, not near fragmentation.

# References
- Scholz, K., Townsend, M., Huber, C., et al. (2023), Investigating the Impact of an Exsolved H2O-CO2 Phase on Magma Chamber Growth, Geochemistry, Geophysics, Geosystems 24, https://doi.org/10.1029/2023GC011151
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
- Liu, Y., Zhang, Y., Behrens, H. (2005), Solubility of H2O in rhyolitic melts at low pressures and a new empirical model for mixed H2O-CO2 solubility, JVGR 143, 219-235, https://doi.org/10.1016/j.jvolgeores.2004.09.019
"""
struct Mafic_Solubility{T, U1, U2} <: AbstractSolubility{T}
    coeffs::NTuple{14, T}     # b1..b10 (H2O), c1..c4 (CO2); dimensionless
    T0::GeoUnit{T, U1}        # Kelvin -> Celsius reference offset
    Tref::GeoUnit{T, U1}      # temperature scale of the fit (1 K), difference and absolute
    Pref::GeoUnit{T, U2}      # pressure scale of the fit (1 MPa)

    function Mafic_Solubility(;
            coeffs = (
                2.99622526644026, 0.00322422830627781, -9.1389095360385,
                0.0336065247530767, 0.00747236662935722, -0.0000150329805347769,
                -0.01233608521548, -4.14842647942619e-6, -0.655454303068124,
                -7.35270395041104e-6, 5668.0, -55.99, 0.4133, 2.041e-3,
            ),
            T0 = 273.15K, Tref = 1.0K, Pref = 1.0e6Pa
        )
        T0U = convert(GeoUnit, T0)
        TrefU = convert(GeoUnit, Tref)
        PrefU = convert(GeoUnit, Pref)
        T = typeof(PrefU).types[1]
        U1 = typeof(T0U).types[2]
        U2 = typeof(PrefU).types[2]
        return new{T, U1, U2}(T.(coeffs), T0U, TrefU, PrefU)
    end
    Mafic_Solubility(coeffs, T0, Tref, Pref) = Mafic_Solubility(; coeffs, T0, Tref, Pref)
end
isdimensional(s::Mafic_Solubility) = isdimensional(s.Pref)

function param_info(s::Mafic_Solubility)
    return MaterialParamsInfo(; Equation = L"m_{H_2O},\ m_{CO_2} = f_{mafic}(P, T, X_{co2})")
end

@inline function compute_dissolved(s::Mafic_Solubility, P, T, X_co2)
    Tc = precision_of(P)
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, c1, c2, c3, c4 = convert_precision(Tc, s.coeffs)
    @unpack_like P Tc T0, Tref, Pref = s
    Tn = (T - T0) / Tref            # ∝ T in °C
    Pm = P / Pref                   # ∝ P in MPa
    meq = @muladd b1 +
        Tn * (b2 + b8 * Tn + b5 * X_co2 + b6 * Pm) +
        X_co2 * (b3 + b9 * X_co2 + b7 * Pm) +
        Pm * (b4 + b10 * Pm)

    # CO2 block: Liu (2005) rhyolite
    Pw = Pm * (1 - X_co2)
    Pc = Pm * X_co2
    Tr = Tref / T
    Pw_sqrt = sqrt(Pw)        # Pw^0.5
    Pw_15 = Pw * Pw_sqrt    # Pw^1.5
    Cco2 = @muladd Pc * ((c1 + c2 * Pw) * Tr + c3 * Pw_sqrt + c4 * Pw_15)

    # The mafic H2O polynomial has no floor at zero and goes negative at high
    # X_co2 / low P, reachable well inside a magma chamber's P-T range. A
    # negative dissolved mass fraction cannot be physical, so floor rather
    # than extrapolate.
    meq = max(meq, zero(meq))

    return (meq / 100, Cco2 / 1000000)
end

function show(io::IO, s::Mafic_Solubility)
    return print(io, "Mafic H2O-CO2 solubility (CO2 block Liu 2005)")
end

# Sun & Yao (2026) — all melt compositions ----------------------------------
"""
    SunYao2026_Solubility(; oxd_wt = (77.0, 0.1, 12.8, 0.8, 0.1, 0.6, 4.0, 4.6, 0.02), a = (a0..a9), b = (b0..b20))

Coupled H2O–CO2 solubility of Sun & Yao (2026), calibrated on 2090 experiments
on melts from carbonatite to pure silica at 660–1924 °C and 1 bar–6 GPa.
`oxd_wt` is the volatile-free melt composition in wt% of
`(SiO2, TiO2, Al2O3, FeOt, MgO, CaO, Na2O, K2O, P2O5)`, the oxide order of
[`GiordanoMeltViscosity`](@ref GeoParams.MaterialParameters.ConstitutiveRelationships.GiordanoMeltViscosity) with P2O5 in place of H2O; the default is a
high-silica rhyolite. It is converted once to the oxide mole fractions ``X``, which are
renormalized over these 9 oxides: MnO and Fe2O3 are not in the set, so fold Fe2O3 into FeOt.

Dissolved H2O (wt%) is the sum of molecular H2O and hydroxyl (eq. 7),
```math
\\begin{aligned}
    \\ln S_{H_2O,m} &= a_0 \\ln P_{H_2O} + a_1 \\ln T + a_2 \\frac{P^2}{T} + (a_3\\sqrt{P} + a_4)\\frac{X_{Na_2O}}{T}, \\\\
    \\ln S_{OH^-} &= a_5 \\ln P_{H_2O} + a_6 \\ln T + (a_7 + a_8 X_{Al_2O_3} + a_9 X_{Na_2O})\\frac{P}{T},
\\end{aligned}
```
and dissolved CO2 (ppm) is the sum of molecular CO2 and carbonate (eq. 8),
```math
\\begin{aligned}
    \\ln S_{CO_2,m} &= b_0 \\ln P_{CO_2} + b_1 \\ln P + b_2 X_{Al_2O_3}\\sqrt{P} + b_3 X_{MgO}^2 + b_4 X_{K_2O}^2, \\\\
    \\ln S_{CO_3^{2-}} &= b_0 \\ln P_{CO_2} + b_5 \\ln T + b_6 \\ln P + \\Omega + F + Q\\, S_{H_2O},
\\end{aligned}
```
with ``\\Omega = (b_7 X_{SiO_2} + b_8 X_{CaO})P/T + (b_9 X_{CaO} + b_{10} X_{Na_2O})\\sqrt{P}``,
``F = b_{11} X_{SiO_2} + b_{12} X_{Al_2O_3} + b_{13} X_{FeO_t} + b_{14} X_{MgO} + b_{15} X_{K_2O} + b_{16}\\ln[1 - (X_{SiO_2} + X_{TiO_2} + X_{Al_2O_3} + X_{P_2O_5})]``
and ``Q = b_{17} X_{SiO_2} + b_{18} X_{Al_2O_3} + b_{19} X_{Na_2O} + b_{20} X_{K_2O}``.
``T`` is in K, and the total pressure ``P`` and the partial pressures
``P_{H_2O} = P(1-X_{co2})``, ``P_{CO_2} = P X_{co2}`` are in bar; the fit uses
partial pressures in place of fugacities. The defaults of `a` and `b` are Table 1.
[`compute_dissolved`](@ref) returns the totals as mass fractions; the species
split is not returned. Pressure and temperature enter through the reference
`GeoUnit`s `Pref = 1 bar` and `Tref = 1 K`, so the law nondimensionalizes like
[`Liu2005_Solubility`](@ref).

Pure gases give exact zeros: no H2O at `X_co2 = 1` and no CO2 at `X_co2 = 0`, with zero
derivatives in `P` and `T`. The slopes in `X_co2` are infinite at these endpoints
(``a_0, a_5, b_0 < 1``), so [`∂dissolved_∂Xco2`](@ref) is not finite there. Above ~2 GPa
the term ``Q\\, S_{H_2O}`` makes the dissolved CO2 of silicic melts exceed a mass fraction
of 1; this is the paper's formula outside its validity.

# References
- Sun, C., Yao, L. (2026), A unified H2O–CO2 solubility–speciation model for magmatic liquids: Constraints on magma storage architecture in continental rifts, EPSL 689, 120115, https://doi.org/10.1016/j.epsl.2026.120115
"""
struct SunYao2026_Solubility{T, U1, U2} <: AbstractSolubility{T}
    X::NTuple{9, T}           # volatile-free oxide mole fractions, in the order of oxd_wt
    a::NTuple{10, T}          # a0..a9 (H2O, eq. 7)
    b::NTuple{21, T}          # b0..b20 (CO2, eq. 8)
    Pref::GeoUnit{T, U1}      # pressure scale of the fit (1 bar)
    Tref::GeoUnit{T, U2}      # temperature scale of the fit (1 K)

    function SunYao2026_Solubility(X, a, b, Pref, Tref)
        PrefU = convert(GeoUnit, Pref)
        TrefU = convert(GeoUnit, Tref)
        T = typeof(PrefU).types[1]
        U1 = typeof(PrefU).types[2]
        U2 = typeof(TrefU).types[2]
        return new{T, U1, U2}(NTuple{9, T}(X), NTuple{10, T}(a), NTuple{21, T}(b), PrefU, TrefU)
    end
end
function SunYao2026_Solubility(;
        oxd_wt = (77.0, 0.1, 12.8, 0.8, 0.1, 0.6, 4.0, 4.6, 0.02),
        a = (
            0.807805, -0.700336, 2.949559e-6, -74.546597, 7146.275452,
            0.449387, -0.360194, -0.308259, 2.984703, -1.69371,
        ),
        b = (
            0.711812, 0.285783, -0.04444, -1019.065935, -451.908354, 1.282791, 0.651135,
            0.040306, 0.583632, -0.116595, -0.036834, -9.761378, -9.184586, -18.300993,
            -11.489536, 9.856334, 3.126584, 0.58461, -2.027465, 0.552689, -1.570268,
        )
    )
    length(oxd_wt) == 9 || throw(ArgumentError("SunYao2026_Solubility: oxd_wt must hold 9 oxide wt% (SiO2 TiO2 Al2O3 FeOt MgO CaO Na2O K2O P2O5), got $(length(oxd_wt))"))
    all(>=(0), oxd_wt) || throw(ArgumentError("SunYao2026_Solubility: oxide wt% must be non-negative, got $oxd_wt"))
    # g/mol of SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O P2O5
    n = Tuple(oxd_wt) ./ (60.0843, 79.8658, 101.961276, 71.8444, 40.3044, 56.0774, 61.97894, 94.196, 141.9445)
    sum(n) > 0 || throw(ArgumentError("SunYao2026_Solubility: oxide wt% must not all be zero"))
    return SunYao2026_Solubility(n ./ sum(n), a, b, 1.0e5Pa, 1.0K)
end
isdimensional(s::SunYao2026_Solubility) = isdimensional(s.Pref)

function param_info(s::SunYao2026_Solubility)
    return MaterialParamsInfo(; Equation = L"m_{H_2O},\ m_{CO_2} = f_{SunYao2026}(P, T, X_{co2}, X_{oxides})")
end

@inline function compute_dissolved(s::SunYao2026_Solubility, P, T, X_co2)
    Tc = precision_of(P)
    a0, a1, a2, a3, a4, a5, a6, a7, a8, a9 = convert_precision(Tc, s.a)
    b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20 =
        convert_precision(Tc, s.b)
    SiO2, TiO2, Al2O3, FeO, MgO, CaO, Na2O, K2O, P2O5 = convert_precision(Tc, s.X)
    @unpack_like P Tc Pref, Tref = s
    Pb = P / Pref                   # ∝ P in bar
    Tk = T / Tref                   # ∝ T in K
    Pw, Pc = Pb * (1 - X_co2), Pb * X_co2
    lnPw, lnPc = log(Pw), log(Pc)
    lnP, lnT, sqrtP, P_T = log(Pb), log(Tk), sqrt(Pb), Pb / Tk

    S_h2o = exp(a0 * lnPw + a1 * lnT + a2 * Pb * P_T + (a3 * sqrtP + a4) * Na2O / Tk) +
        exp(a5 * lnPw + a6 * lnT + (a7 + a8 * Al2O3 + a9 * Na2O) * P_T)
    # a zero partial pressure gives log(0) = -Inf, whose ForwardDiff partials are NaN;
    # the guard keeps the value and its P and T derivatives exactly zero
    S_h2o = ifelse(iszero(Pw), zero(S_h2o), S_h2o)

    Ω = (b7 * SiO2 + b8 * CaO) * P_T + (b9 * CaO + b10 * Na2O) * sqrtP
    F = b11 * SiO2 + b12 * Al2O3 + b13 * FeO + b14 * MgO + b15 * K2O + b16 * log(1 - (SiO2 + TiO2 + Al2O3 + P2O5))
    Q = b17 * SiO2 + b18 * Al2O3 + b19 * Na2O + b20 * K2O
    S_co2 = exp(b0 * lnPc + b1 * lnP + b2 * Al2O3 * sqrtP + b3 * MgO^2 + b4 * K2O^2) +
        exp(b0 * lnPc + b5 * lnT + b6 * lnP + Ω + F + Q * S_h2o)
    S_co2 = ifelse(iszero(Pc), zero(S_co2), S_co2)

    return (S_h2o / 100, S_co2 / 1000000)   # wt% -> fraction ; ppm -> fraction
end

function show(io::IO, s::SunYao2026_Solubility)
    return print(io, "Sun & Yao (2026) H2O-CO2 solubility")
end

# Gas-mixture properties -----------------------------------------------------
"""
    GasMixture(; Cp_h2o=3880J/kg/K, Cp_co2=1200J/kg/K, M_h2o=18.02e-3kg/mol, M_co2=44.01e-3kg/mol)

H2O–CO2 gas-mixture properties keyed on the CO2 mole fraction of the gas
`X_co2` (Degruyter & Huber 2014). The effective molar mass is
```math
    m_g = M_{H_2O}(1-X_{co2}) + M_{CO_2} X_{co2}
```
and the mass-weighted specific heat is
```math
    c_g = \\frac{M_{H_2O} c_{H_2O}(1-X_{co2}) + M_{CO_2} c_{CO_2} X_{co2}}{m_g}
```
with the reference convention ``c_g = 0`` at ``X_{co2}=0`` — a discontinuity,
since the formula's own limit as ``X_{co2} \\to 0`` is ``c_{H_2O}``.

All fields are `GeoUnit`s, so the struct nondimensionalizes, but both accessors
read values rather than `Quantity`s: neither return carries units.
"""
@with_kw_noshow struct GasMixture{T, U1, U2} <: AbstractMaterialParam
    Cp_h2o::GeoUnit{T, U1} = 3880.0J / kg / K   # specific heat of H2O gas
    Cp_co2::GeoUnit{T, U1} = 1200.0J / kg / K   # specific heat of CO2 gas
    M_h2o::GeoUnit{T, U2} = 18.02e-3kg / mol    # molar mass H2O
    M_co2::GeoUnit{T, U2} = 44.01e-3kg / mol    # molar mass CO2
end
GasMixture(args...) = promote_construct(GasMixture, args...)
isdimensional(s::GasMixture) = isdimensional(s.Cp_h2o)

function param_info(s::GasMixture)
    return MaterialParamsInfo(; Equation = L"c_g = (M_{H_2O} c_{H_2O}(1-X) + M_{CO_2} c_{CO_2} X)/m_g")
end

"""
    effective_molar_mass(s::GasMixture, X_co2)

Effective molar mass of the H2O–CO2 gas mixture at CO2 mole fraction `X_co2`.
"""
@inline function effective_molar_mass(s::GasMixture, X_co2)
    Tc = precision_of(X_co2)
    @unpack_val Tc M_h2o, M_co2 = s   # X_co2 is dimensionless; a nondimensionalized struct yields nondim output
    return M_h2o * (1 - X_co2) + M_co2 * X_co2
end

"""
    compute_gas_heatcapacity(s::GasMixture, X_co2)

Mass-weighted specific heat of the H2O–CO2 gas mixture; zero at `X_co2 == 0`
(reference convention).
"""
@inline function compute_gas_heatcapacity(s::GasMixture, X_co2)
    Tc = precision_of(X_co2)
    @unpack_val Tc Cp_h2o, Cp_co2, M_h2o, M_co2 = s
    iszero(X_co2) && return zero(Cp_h2o * X_co2)
    m_g = M_h2o * (1 - X_co2) + M_co2 * X_co2
    return (M_h2o * Cp_h2o * (1 - X_co2) + M_co2 * Cp_co2 * X_co2) * inv(m_g)
end

# Exsolved gas fraction ------------------------------------------------------
"""
    compute_gas_fraction(s::AbstractSolubility, ρ::ThreePhase_Density, args)

Equilibrium volume fraction of exsolved H2O gas in a melt + crystal + gas magma, the
`ϕ_gas` input of [`ThreePhase_Density`](@ref), [`Suspension`](@ref GeoParams.MaterialParameters.ConstitutiveRelationships.Suspension) and
[`CriticalPorosityFragmentation`](@ref GeoParams.Fragmentation.CriticalPorosityFragmentation). `args` holds the pressure `P`, the temperature
`T`, the melt fraction `ϕ` of the condensed (melt + crystal) phases as returned by a
melting law, and the total water mass fraction `m_h2o_total` of the melt + gas, i.e.
on a crystal-free basis; crystals hold no water. `args` is also passed to the `ρmelt`
and `ρgas` laws of `ρ`. Also callable as `compute_gas_fraction(s, ρ; args...)`.

The gas is pure H2O, so the dissolved water mass fraction of the melt is
``m_{eq} = `` `compute_dissolved(s, P, T, 0)[1]`. The water mass balance of
Degruyter & Huber (2014, eq. 18) gives the gas mass fraction of the melt + gas,
```math
    x_g = \\max\\left(0, \\frac{m_{h2o,total} - m_{eq}}{1 - m_{eq}}\\right),
```
and with the crystal volume fraction ``\\phi_x = (1 - \\phi_{gas})(1 - \\phi)`` of
their melting curve (eq. A.4) the gas volume fraction of the whole magma is
```math
    \\phi_{gas} = \\frac{\\phi\\, x_g \\rho_m}{\\phi\\, x_g \\rho_m + (1 - x_g)\\rho_g}.
```
Pass `ϕ_gas` together with `ϕ_x = (1 - ϕ_gas) * (1 - ϕ)` to `ThreePhase_Density` and
`Suspension`.

Undersaturated magma (``m_{h2o,total} \\le m_{eq}``) gives zero without evaluating
`ρgas`. For the H2O gas density use [`HalbachChatterjee_Density`](@ref GeoParams.MaterialParameters.Density.HalbachChatterjee_Density),
which is within 0.8 % of IAPWS-95 for 0.1–400 MPa and 600–1200 °C. The result is `NaN` for
`m_h2o_total` outside ``[0, 1)`` or `ϕ` outside ``[0, 1]``. It is dimensionless and the
same for dimensional and nondimensionalized laws and inputs.

# Example
```julia
ρ = ThreePhase_Density(; ρgas = HalbachChatterjee_Density())
args = (; P = 1.0e8, T = 1100.0, ϕ = 0.6, m_h2o_total = 0.05)
ϕ_gas = compute_gas_fraction(SunYao2026_Solubility(), ρ, args)
ρ_magma = compute_density(ρ, (; args..., ϕ_gas, ϕ_x = (1 - ϕ_gas) * (1 - args.ϕ)))
```

# References
- Degruyter, W., Huber, C. (2014), A model for the eruption frequency of upper crustal silicic magma chambers, EPSL 403, 117-130, https://doi.org/10.1016/j.epsl.2014.06.047
"""
@inline function compute_gas_fraction(s::AbstractSolubility, ρ::ThreePhase_Density; ϕ, m_h2o_total, P, T, kwargs...)
    m_eq = compute_dissolved(s, P, T, zero(precision_of(P)))[1]
    0 <= m_h2o_total < 1 && 0 <= ϕ <= 1 || return oftype(m_eq, NaN)
    x_g = max((m_h2o_total - m_eq) / (1 - m_eq), zero(m_eq))
    # A zero gas fraction must not depend on ρgas, which may be NaN where it is undefined.
    iszero(x_g) && return x_g
    args = (; P, T, kwargs...)
    scaled_gas_volume = ϕ * x_g * compute_density(ρ.ρmelt, args)
    return scaled_gas_volume / (scaled_gas_volume + (1 - x_g) * compute_density(ρ.ρgas, args))
end
@inline compute_gas_fraction(s::AbstractSolubility, ρ::ThreePhase_Density, args::NamedTuple) = compute_gas_fraction(s, ρ; args...)

# Shared keyword / NamedTuple entry points -----------------------------------
@inline function compute_dissolved(s::AbstractSolubility; P = 0.0e0, T = 0.0e0, X_co2 = zero(precision_of(P)), kwargs...)
    return compute_dissolved(s, P, T, X_co2)
end
@inline compute_dissolved(s::AbstractSolubility, args::NamedTuple) = compute_dissolved(s; args...)

# Phase-struct dispatch; a phase with no solubility model contributes nothing.
@inline function compute_dissolved(s::AbstractMaterialParamsStruct, args)
    isempty(s.Solubility) && return (0.0e0, 0.0e0)
    return compute_dissolved(s.Solubility[1], args)
end

# ForwardDiff partials — replace the reference's hand-coded analytic partials.
# Each returns (∂m_h2o/∂·, ∂m_co2/∂·).
_dissolved_svec(s::AbstractSolubility, P, T, X_co2) = SVector(compute_dissolved(s, P, T, X_co2)...)

"""
    ∂dissolved_∂P(s, P, T, X_co2) -> (∂m_h2o/∂P, ∂m_co2/∂P)

ForwardDiff partial derivatives of [`compute_dissolved`](@ref) with respect to
pressure. Companions [`∂dissolved_∂T`](@ref), [`∂dissolved_∂Xco2`](@ref).
"""
∂dissolved_∂P(s::AbstractSolubility, P, T, X_co2) = Tuple(ForwardDiff.derivative(p -> _dissolved_svec(s, p, T, X_co2), P))

"""
    ∂dissolved_∂T(s, P, T, X_co2) -> (∂m_h2o/∂T, ∂m_co2/∂T)
"""
∂dissolved_∂T(s::AbstractSolubility, P, T, X_co2) = Tuple(ForwardDiff.derivative(t -> _dissolved_svec(s, P, t, X_co2), T))

"""
    ∂dissolved_∂Xco2(s, P, T, X_co2) -> (∂m_h2o/∂X_co2, ∂m_co2/∂X_co2)
"""
∂dissolved_∂Xco2(s::AbstractSolubility, P, T, X_co2) = Tuple(ForwardDiff.derivative(x -> _dissolved_svec(s, P, T, x), X_co2))

"""
    find_Xco2(s::AbstractSolubility, P, T, m_h2o_target; X0=0.5, tol=1e-8, max_iter=50) -> X_co2

Invert [`compute_dissolved`](@ref) for the gas composition: solve for
`X_co2 ∈ [0,1]` such that `compute_dissolved(s, P, T, X_co2)[1] ==
m_h2o_target`, at fixed pressure `P` and temperature `T`. Uses a safeguarded
Newton iteration (via [`∂dissolved_∂Xco2`](@ref)) bracketed by bisection, so
an out-of-bounds Newton step always falls back to a bisection halving
instead of leaving `[0,1]`.

Throws an `ErrorException` if `m_h2o_target` is infeasible at this `P,T`
(outside the achievable range between `X_co2=0` and `X_co2=1`) or if the
iteration fails to converge within `max_iter` steps — this never returns a
clamped or out-of-tolerance answer.
"""
function find_Xco2(s::AbstractSolubility, P, T, m_h2o_target; X0 = 0.5, tol = 1.0e-8, max_iter = 50)
    m_h2o_target ≥ 0 || throw(DomainError(m_h2o_target, "find_Xco2: a dissolved H2O mass fraction cannot be negative"))
    f(x) = compute_dissolved(s, P, T, x)[1] - m_h2o_target
    lo, hi = 0.0, 1.0
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0 && abs(flo) > tol && abs(fhi) > tol
        error("find_Xco2: m_h2o_target=$m_h2o_target is infeasible at P=$P, T=$T (residuals $flo at X_co2=0, $fhi at X_co2=1 have the same sign)")
    end

    x = clamp(X0, lo, hi)
    iter, fx = 0, f(x)
    while abs(fx) > tol && iter < max_iter
        iter += 1
        # keep [lo,hi] bracketing a sign change, regardless of whether
        # compute_dissolved is increasing or decreasing in X_co2
        if fx * flo > 0
            lo, flo = x, fx
        else
            hi = x
        end
        fp = ∂dissolved_∂Xco2(s, P, T, x)[1]
        x_newton = iszero(fp) ? NaN : x - fx / fp
        x = (isnan(x_newton) || x_newton < lo || x_newton > hi) ? (lo + hi) / 2 : x_newton
        fx = f(x)
    end
    abs(fx) > tol && error("find_Xco2: iterations did not converge (|residual|=$(abs(fx)) after $max_iter iterations)")
    return x
end

"""
    compute_saturation_pressure(s::AbstractSolubility, m_h2o, m_co2, T; tol=100eps, max_iter=50) -> (P, X_co2)

Fluid-saturation pressure `P` of a melt with dissolved H2O and CO2 mass fractions
`m_h2o` and `m_co2` at temperature `T`, and the CO2 mole fraction `X_co2` of the
coexisting H2O–CO2 gas: the inverse of [`compute_dissolved`](@ref) in `P` and `X_co2`.
`P` is in Pa and `T` in K, or both nondimensional for a nondimensionalized law. Applied to
the volatile contents of melt inclusions, this is the saturation barometer of
Sun & Yao (2026, sect. 3–4).

`m_co2 = 0` gives a pure H2O gas (`X_co2 = 0`) and `m_h2o = 0` a pure CO2 gas
(`X_co2 = 1`), with `P` matching the nonzero content alone. Otherwise a Newton
iteration in ``(\\ln P, X_{co2})`` drives both misfits ``\\ln(m/m_{target})`` below `tol`,
using [`∂dissolved_∂P`](@ref) and [`∂dissolved_∂Xco2`](@ref). Each step changes `P` by at
most a factor ``e`` and moves `X_co2` at most halfway towards 0 or 1.

The solve is iterative and meant for analysis, such as melt-inclusion barometry or
setting initial conditions, not for evaluation at every grid point of a solver. It throws
a `DomainError` for negative or all-zero contents and an `ErrorException` if it does not
converge within `max_iter` iterations; it never returns an unconverged answer. Outside its
calibration, a polynomial law such as [`Mafic_Solubility`](@ref) can match the same
contents at more than one pressure.

For rhyolite with [`SunYao2026_Solubility`](@ref) the solve converges up to ~1 GPa. Above
~1.5 GPa with small `X_co2` it can fail although a solution exists, because the model's
``\\exp(Q\\, S_{H_2O})`` term diverges.

# Example
```julia
P, X_co2 = compute_saturation_pressure(SunYao2026_Solubility(), 0.04, 0.001, 1173.15)
# (2.76e8, 0.54): 4 wt% H2O and 1000 ppm CO2 saturate at 276 MPa in rhyolite at 900 °C
```

# References
- Sun, C., Yao, L. (2026), A unified H2O–CO2 solubility–speciation model for magmatic liquids: Constraints on magma storage architecture in continental rifts, EPSL 689, 120115, https://doi.org/10.1016/j.epsl.2026.120115
"""
function compute_saturation_pressure(
        s::AbstractSolubility, m_h2o, m_co2, T;
        tol = 100 * eps(precision_of((m_h2o, m_co2, T))), max_iter = 50
    )
    m_h2o ≥ 0 && m_co2 ≥ 0 && m_h2o + m_co2 > 0 ||
        throw(DomainError((m_h2o, m_co2), "compute_saturation_pressure: dissolved mass fractions must be non-negative and not both zero"))
    Tc = precision_of((m_h2o, m_co2, T))
    @unpack_val Tc Pref = s
    lnP = log(100 * Pref)   # 10–100 MPa for the laws here, in the units of their P
    if iszero(m_co2) || iszero(m_h2o)
        # pure gas: one equation in ln P at fixed X_co2
        i, X, m = iszero(m_co2) ? (1, zero(Tc), m_h2o) : (2, one(Tc), m_co2)
        for _ in 1:max_iter
            P = exp(lnP)
            mP = compute_dissolved(s, P, T, X)[i]
            F = log(mP / m)
            abs(F) ≤ tol && return (P, X)
            # ∂/∂P can be NaN at X_co2 = 1, where the H2O partial pressure is zero (Liu2005's sqrt(Pw)), so the slope is taken just inside
            dmP = ∂dissolved_∂P(s, P, T, X == 1 ? prevfloat(X) : X)[i]
            lnP -= clamp(F * mP / (P * dmP), -1, 1)
        end
    else
        X = Tc(0.5)
        for _ in 1:max_iter
            P = exp(lnP)
            h, c = compute_dissolved(s, P, T, X)
            F = SVector(log(h / m_h2o), log(c / m_co2))
            maximum(abs, F) ≤ tol && return (P, X)
            hP, cP = ∂dissolved_∂P(s, P, T, X)
            hX, cX = ∂dissolved_∂Xco2(s, P, T, X)
            ΔlnP, ΔX = SMatrix{2, 2}(P * hP / h, P * cP / c, hX / h, cX / c) \ F
            lnP -= clamp(ΔlnP, -1, 1)
            X = clamp(X - ΔX, X / 2, (1 + X) / 2)
        end
    end
    return error("compute_saturation_pressure: no convergence within $max_iter iterations for m_h2o=$m_h2o, m_co2=$m_co2, T=$T")
end

# Two output arrays, so this cannot route through the single-array `compute_param!`.
"""
    compute_dissolved!(m_h2o, m_co2, MatParam, Phases, args)

In-place dissolved H2O and CO2 over a domain. `args` is a NamedTuple of `P, T,
X_co2` index-matched to the output arrays. Also accepts a single solubility or
phase struct in place of `(MatParam, Phases)`.
"""
function compute_dissolved!(
        m_h2o::AbstractArray, m_co2::AbstractArray,
        MatParam::NTuple{N, AbstractMaterialParamsStruct}, Phases::AbstractArray, args
    ) where {N}
    for I in eachindex(m_h2o, m_co2, Phases)
        argsi = (; zip(keys(args), getindex.(values(args), I))...)
        m_h2o[I], m_co2[I] = compute_dissolved(MatParam, Phases[I], argsi)
    end
    return nothing
end

function compute_dissolved!(
        m_h2o::AbstractArray, m_co2::AbstractArray,
        s::Union{AbstractSolubility, AbstractMaterialParamsStruct}, args
    )
    for I in eachindex(m_h2o, m_co2)
        argsi = (; zip(keys(args), getindex.(values(args), I))...)
        m_h2o[I], m_co2[I] = compute_dissolved(s, argsi)
    end
    return nothing
end

# Integer-phase selection / single struct over a domain (mirrors compute_meltfraction).
compute_dissolved(MatParam, arg, args::Vararg{Any, N}) where {N} = compute_param(compute_dissolved, MatParam, arg, args...)

# Phase-ratio mix of both outputs (mirrors compute_meltfraction_ratio). The
# shared compute_param_times_frac sums a scalar, so the two-element output gets
# its own unrolled dot-product, evaluating each phase once.
"""
    compute_dissolved_ratio(PhaseRatios, MatParam::NTuple{N,AbstractMaterialParamsStruct}, args) -> (m_h2o, m_co2)

Dissolved H2O and CO2 mass fractions at a point whose composition is a mixture of
`N` phases. `PhaseRatios[i]` is the volume fraction of phase `i`, summing to one,
and the result is the fraction-weighted average of [`compute_dissolved`](@ref)
over all phases. Each phase is evaluated exactly once.
"""
compute_dissolved_ratio(args::Vararg{Any, N}) where {N} = compute_dissolved_times_frac(compute_dissolved, args...)

@generated function compute_dissolved_times_frac(
        fn::F, PhaseRatios::Union{NTuple{N}, SVector{N}}, MatParam::NTuple{N, AbstractMaterialParamsStruct}, argsi
    ) where {F <: Function, N}
    return quote
        mh = zero($(eltype(PhaseRatios)))
        mc = zero($(eltype(PhaseRatios)))
        Base.@nexprs $N i -> begin
            @inline
            hᵢ, cᵢ = fn(MatParam[i], argsi)
            r = @inbounds PhaseRatios[i]
            mh += hᵢ * r
            mc += cᵢ * r
        end
        return (mh, mc)
    end
end

end
