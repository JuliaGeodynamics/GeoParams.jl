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
import Base.show, GeoParams.param_info

include("../Computations.jl")

abstract type AbstractSolubility{T} <: AbstractMaterialParam end

export compute_dissolved, # (m_h2o, m_co2) mass fractions
    compute_dissolved!, # in-place, fills two arrays
    compute_dissolved_ratio, # phase-ratio mix of both outputs
    ∂dissolved_∂P, # ForwardDiff partials, both outputs
    ∂dissolved_∂T,
    ∂dissolved_∂Xco2,
    find_Xco2, # Newton-invert compute_dissolved for X_co2 given target m_h2o
    param_info,
    AbstractSolubility,
    Liu2005_Solubility, # silicic (rhyolite)
    Mafic_Solubility, # mafic (basalt), CO2 block still Liu rhyolite
    GasMixture, # H2O-CO2 gas-mixture properties
    compute_gas_heatcapacity, # mass-weighted specific heat c_g(X_co2)
    effective_molar_mass            # m_g(X_co2)

# Empty routine when no solubility is defined
compute_dissolved(::No_MaterialParam{_T}; kwargs...) where {_T} = (zero(_T), zero(_T))

# Liu et al. (2005) — silicic ------------------------------------------------
"""
    Liu2005_Solubility(; coeffs=(b1..b6, c1..c4), Pref=1MPa, Tref=1K)

Coupled H2O–CO2 solubility for silicic (rhyolite) melt after Liu et al. (2005),
as used by Degruyter & Huber (2014). [`compute_dissolved`](@ref) returns the
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
    b1, b2, b3, b4, b5, b6, c1, c2, c3, c4 = s.coeffs
    if P isa Quantity
        @unpack_units Pref, Tref = s
    else
        @unpack_val Pref, Tref = s
    end
    Pw = P * (1 - X_co2) / Pref     # ∝ H2O partial pressure in MPa
    Pc = P * X_co2 / Pref           # ∝ CO2 partial pressure in MPa
    Tr = Tref / T                   # ∝ 1/T [K]

    Pw_sqrt = sqrt(Pw)     # Pw^0.5
    Pw_15 = Pw * Pw_sqrt     # Pw^1.5

    meq = @muladd (b1 * Pw_sqrt + b2 * Pw + b3 * Pw_15) * Tr + b4 * Pw_15 + Pc * (b5 * Pw_sqrt + b6 * Pw)
    Cco2 = @muladd Pc * ((c1 + c2 * Pw) * Tr + c3 * Pw_sqrt + c4 * Pw_15)

    return (1.0e-2 * meq, 1.0e-6 * Cco2)   # wt% -> fraction ; ppm -> fraction
end

function show(io::IO, s::Liu2005_Solubility)
    return print(io, "Liu et al. (2005) silicic H2O-CO2 solubility")
end

# Mafic ----------------------------------------------------------------------
"""
    Mafic_Solubility(; coeffs=(b1..b10, c1..c4), T0=273.15K, Tref=1K, Pref=1MPa)

Coupled H2O–CO2 solubility for mafic (basalt) melt. Dissolved H2O follows the
mafic polynomial of the Scholz/Degruyter–Huber reference in the dimensionless
groups ``T_C = (T-T_0)/T_{ref}`` (numerically °C) and ``P_m = P/P_{ref}``
(numerically MPa); dissolved CO2 reuses the Liu (2005) rhyolite CO2 block.
Nondimensionalizes like [`Liu2005_Solubility`](@ref).

# References
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
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, c1, c2, c3, c4 = s.coeffs
    if P isa Quantity
        @unpack_units T0, Tref, Pref = s
    else
        @unpack_val T0, Tref, Pref = s
    end
    Tc = (T - T0) / Tref            # ∝ T in °C
    Pm = P / Pref                   # ∝ P in MPa
    meq = @muladd b1 +
        Tc * (b2 + b8 * Tc + b5 * X_co2 + b6 * Pm) +
        X_co2 * (b3 + b9 * X_co2 + b7 * Pm) +
        Pm * (b4 + b10 * Pm)

    # CO2 block: Liu (2005) rhyolite
    Pw = Pm * (1 - X_co2)
    Pc = Pm * X_co2
    Tr = Tref / T
    Pw_sqrt = sqrt(Pw)        # Pw^0.5
    Pw_15 = Pw * Pw_sqrt    # Pw^1.5
    Cco2 = @muladd Pc * ((c1 + c2 * Pw) * Tr + c3 * Pw_sqrt + c4 * Pw_15)

    return (1.0e-2 * meq, 1.0e-6 * Cco2)
end

function show(io::IO, s::Mafic_Solubility)
    return print(io, "Mafic H2O-CO2 solubility (CO2 block Liu 2005)")
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
with the reference convention ``c_g = 0`` at ``X_{co2}=0``. All fields are
`GeoUnit`s, so the struct nondimensionalizes.
"""
@with_kw_noshow struct GasMixture{T, U1, U2} <: AbstractMaterialParam
    Cp_h2o::GeoUnit{T, U1} = 3880.0J / kg / K   # specific heat of H2O gas
    Cp_co2::GeoUnit{T, U1} = 1200.0J / kg / K   # specific heat of CO2 gas
    M_h2o::GeoUnit{T, U2} = 18.02e-3kg / mol    # molar mass H2O
    M_co2::GeoUnit{T, U2} = 44.01e-3kg / mol    # molar mass CO2
end
GasMixture(args...) = GasMixture(convert.(GeoUnit, args)...)
isdimensional(s::GasMixture) = isdimensional(s.Cp_h2o)

function param_info(s::GasMixture)
    return MaterialParamsInfo(; Equation = L"c_g = (M_{H_2O} c_{H_2O}(1-X) + M_{CO_2} c_{CO_2} X)/m_g")
end

"""
    effective_molar_mass(s::GasMixture, X_co2)

Effective molar mass of the H2O–CO2 gas mixture at CO2 mole fraction `X_co2`.
"""
@inline function effective_molar_mass(s::GasMixture, X_co2)
    @unpack_val M_h2o, M_co2 = s   # X_co2 is dimensionless; a nondimensionalized struct yields nondim output
    return M_h2o * (1 - X_co2) + M_co2 * X_co2
end

"""
    compute_gas_heatcapacity(s::GasMixture, X_co2)

Mass-weighted specific heat of the H2O–CO2 gas mixture; zero at `X_co2 == 0`
(reference convention).
"""
@inline function compute_gas_heatcapacity(s::GasMixture, X_co2)
    @unpack_val Cp_h2o, Cp_co2, M_h2o, M_co2 = s
    iszero(X_co2) && return zero(Cp_h2o * X_co2)
    m_g = M_h2o * (1 - X_co2) + M_co2 * X_co2
    return (M_h2o * Cp_h2o * (1 - X_co2) + M_co2 * Cp_co2 * X_co2) * inv(m_g)
end

# Shared keyword / NamedTuple entry points -----------------------------------
@inline function compute_dissolved(s::AbstractSolubility; P = 0.0e0, T = 0.0e0, X_co2 = 0.0e0, kwargs...)
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
compute_dissolved(args::Vararg{Any, N}) where {N} = compute_param(compute_dissolved, args...)

# Phase-ratio mix of both outputs (mirrors compute_meltfraction_ratio). The
# shared compute_param_times_frac sums a scalar, so the two-element output gets
# its own unrolled dot-product, evaluating each phase once.
compute_dissolved_ratio(args::Vararg{Any, N}) where {N} = compute_dissolved_times_frac(compute_dissolved, args...)

@generated function compute_dissolved_times_frac(
        fn::F, PhaseRatios::Union{NTuple{N, T}, SVector{N, T}}, MatParam::NTuple{N, AbstractMaterialParamsStruct}, argsi
    ) where {F <: Function, N, T}
    return quote
        mh = zero($T)
        mc = zero($T)
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
