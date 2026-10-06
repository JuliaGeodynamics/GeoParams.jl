module SeismicVelocity

# This implements different methods to compute seismic velocities
#
# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful, MuladdMacro
using ..Units
using ..PhaseDiagrams
using ..MaterialParameters: MaterialParamsInfo
using Statistics
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
using GeoParams: LinearInterpolator, interpolate

using Roots
import Base.show, GeoParams.param_info

abstract type AbstractSeismicVelocity{T} <: AbstractMaterialParam end

export compute_wave_velocity, # calculation routines
    compute_wave_velocity!, # calculation routines
    ConstantSeismicVelocity, # constant
    melt_correction,
    porosity_correction,
    anelastic_correction,
    param_info,
    correct_wavevelocities_phasediagrams,
    melt_correction_Takei

include("../Utils.jl")
include("../Computations.jl")

# Constant Velocity -------------------------------------------------------
"""
    ConstantSeismicVelocity(Vp=8.1 km/s, Vs=4.5km/s)

Set a constant seismic P and S-wave velocity:
```math
    V_p = cst
```
```math
    V_s = cst
```
where ``V_p, V_s`` are the P-wave and S-wave velocities [``km/s``].
"""
@with_kw_noshow struct ConstantSeismicVelocity{T, U} <: AbstractSeismicVelocity{T}
    Vp::GeoUnit{T, U} = 8.1e3m / s               # P-wave velocity
    Vs::GeoUnit{T, U} = 4.5e3m / s               # S-wave velocity
end
ConstantSeismicVelocity(args...) = promote_construct(ConstantSeismicVelocity, args...)

function param_info(s::ConstantSeismicVelocity) # info about the struct
    return MaterialParamsInfo(; Equation = L"v_p = cst \\ v_s = cst")
end

# Calculation routines
function compute_wave_velocity(s::ConstantSeismicVelocity; wave, kwargs...)
    _T = precision_of(values(kwargs))
    @unpack_val _T Vp, Vs = s
    wave === :Vp && return Vp
    wave === :Vs && return Vs
    wave === :VpVs && return Vp / Vs
    throw(ArgumentError("`wave` must be :Vp, :Vs or :VpVs, got $(repr(wave))"))
end

# Print info
function show(io::IO, g::ConstantSeismicVelocity)
    return print(
        io, "Constant seismic velocity: Vp=$(UnitValue(g.Vp)), Vs=$(UnitValue(g.Vs))"
    )
end
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Phase diagrams

#function param_info(s::PhaseDiagram_LookupTable) # info about the struct
#    return MaterialParamsInfo(Equation = L"Vp = f_{PhaseDiagram}(T,P))" )
#end

"""
    compute_wave_velocity(s::PhaseDiagram_LookupTable; P, T, wave, kwargs...)

Interpolates the seismic wave velocity selected by `wave` (e.g. `:Vp`, `:Vs`, `:VpVs`)
as a function of pressure `P` and temperature `T` from the lookup table `s`.
"""
function compute_wave_velocity(s::PhaseDiagram_LookupTable; P, T, wave, kwargs...)
    fn = getfield(s, wave)
    return fn(T, P)
end
#-------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------#
# Computational routines needed for computations with the MaterialParams structure

function compute_wave_velocity(s::AbstractMaterialParamsStruct, args)
    if isempty(s.SeismicVelocity) #in case there is a phase with no melting parametrization
        return zero(typeof(args).types[1])
    else
        return compute_wave_velocity(s.SeismicVelocity[1], args)
    end
end

#-------------------------------------------------------------------------------------------------------------

#Multiple dispatch to rest of routines found in Computations.jl
for myType in (:ConstantSeismicVelocity, :PhaseDiagram_LookupTable)
    @eval begin
        function compute_wave_velocity(p::$(myType), args)
            return compute_wave_velocity(p::$(myType); args...)
        end
    end
end

"""
    compute_wave_velocity!(V, s, args)

In-place version of [`compute_wave_velocity`](@ref) that fills the array `V` with the seismic wave
velocity over the whole domain.
"""
compute_wave_velocity!(args...) = compute_param!(compute_wave_velocity, args...)
compute_wave_velocity(MatParam, arg, args...) = compute_param(compute_wave_velocity, MatParam, arg, args...)

"""
        Vp_cor,Vs_cor = melt_correction(  Kb_L, Kb_S, Ks_S, ρL, ρS, Vp0, Vs0, ϕ, α)

Corrects P- and S-wave velocities if the rock is partially molten.

Input:
====
- `Kb_L`: adiabatic bulk modulus of melt
- `Kb_S`: adiabatic bulk modulus of the solid phase
- `Ks_S`: shear modulus of the solid phase
- `ρL`  : density of the melt
- `ρS`  : density of the solid phase
- `Vp0` : initial P-wave velocity of the solid phase
- `Vs0` : initial S-wave velocity of the solid phase
- `ϕ`   : melt volume fraction
- `α`   : grain-boundary contiguity of the solid framework, from 0 (grains fully wetted by melt)
          to 1 (no melt on grain boundaries). At textural equilibrium it decreases with melt fraction,
          roughly as `1 - A√ϕ` with `A` ≈ 1 - 2.3, i.e. 0.7 - 0.9 for a few percent melt.

Output:
====
- `Vp_cor,Vs_cor` : corrected P-wave and S-wave velocities for melt fraction

The velocity reduction is first order in `ϕ`. It loses accuracy when the product of `ϕ` and the
slopes `Λ_K`, `Λ_G` of the framework moduli approaches one, which happens for small contiguities.

The routine uses the reduction formulation of Clark & Lesher, (2017) and is based on the equilibrium geometry model for the solid skeleton of Takei et al., 1998.

References:
====

- Takei (1998) Constitutive mechanical relations of solid-liquid composites in terms of grain-boundary contiguity, Journal of Geophysical Research: Solid Earth, Vol(103)(B8), 18183--18203

- Clark & Lesher (2017) Elastic properties of silicate melts: Implications for low velocity zones at the lithosphere-asthenosphere boundary. Science Advances, Vol 3 (12), e1701312

"""
function melt_correction(
        Kb_L::_T, Kb_S::_T, Ks_S::_T, ρL::_T, ρS::_T, Vp0::_T, Vs0::_T, ϕ::_T, α::_T
    ) where {_T <: Number}
    iszero(ϕ) && return Vp0, Vs0

    # Takei 1998: Approximation Formulae for Bulk and Shear Moduli of Isotropic Solid Skeleton
    ν = 0.25                         # poisson ratio

    aij = (
        0.318, 6.78, 57.56, 0.182,
        0.164, 4.29, 26.658, 0.464,
        1.549, 4.814, 8.777, -0.29,
    )
    bij = (
        -0.3238, 0.2341,
        -0.1819, 0.5103,
    )

    # Lines below are equivalent to:
    # a = zeros(3)
    # for i in 1:3
    #     a[i] =
    #         aij[i, 1] * exp(aij[i, 2] * (ν - 0.25) + aij[i, 3] * (ν - 0.25)^3) + aij[i, 4]
    # end
    a = ntuple(Val(3)) do i
        idx = 4 * i - 3 # linear offset index
        aij[idx] * exp(aij[idx + 1] * (ν - 0.25) + aij[idx + 2] * (ν - 0.25)^3) + aij[idx + 3]
    end

    # Lines below are equivalent to:
    # b = zeros(2)
    # for i in 1:2
    #     b[i] = bij[i, 1] * ν + bij[i, 2]
    # end
    b = ntuple(Val(2)) do i
        idx = 2 * i - 1 # linear offset index
        bij[idx] * ν + bij[idx + 1]
    end

    nk = a[1] * α + a[2] * (1.0 - α) + a[3] * α * (1.0 - α) * (0.5 - α)
    nμ = b[1] * α + b[2] * (1.0 - α)

    # computation of the bulk modulus ratio of the skeletal framework over the solid phase
    ksk_k = fastpow(α, nk)
    # computation of the shear modulus ratio of the skeletal framework over the solid phase
    μsk_μ = fastpow(α, nμ)

    # apply correction for the melt fraction to adiabatic bulk and shear modulii
    ksk = ksk_k * Kb_S
    μsk = μsk_μ * Ks_S

    kb = (1.0 - ϕ) * ksk
    μ = (1.0 - ϕ) * μsk

    # slopes of the framework moduli with melt fraction, K_b/K_S = 1 - ΛK ϕ and μ/G = 1 - ΛG ϕ
    # (Clark & Lesher, 2017, eqs. 3-4)
    ΛK = (1 - kb / Kb_S) / ϕ
    ΛG = (1 - μ / Ks_S) / ϕ

    # Seismic wave velocity melt correction Clark et al., 2017
    β = Kb_S / Kb_L
    γ = Ks_S / Kb_S

    # Formulation of the fraction reduction of P-wave and S-wave
    ΔVp =
        (
        (
            (((β - 1.0) * ΛK) / ((β - 1.0) + ΛK) + 4.0 / 3.0 * γ * ΛG) /
                (1.0 + 4.0 / 3.0 * γ)
        ) - (1.0 - ρL / ρS)
    ) * (ϕ * 0.5)
    ΔVs = (ΛG - (1.0 - ρL / ρS)) * (ϕ * 0.5)

    # get the correction values
    Vp_cor = Vp0 - ΔVp * Vp0
    Vs_cor = Vs0 - Vs0 * ΔVs

    return Vp_cor, Vs_cor
end

"""
        Vs_cor = porosity_correction(Kb_S, Ks_S, ρf, ρS, Vs0, depth, α)

Corrects S-wave velocity at shallow depth for fluid-filled porosity, with the porosity taken
from an empirical porosity-depth profile.

Input:
====
- `Kb_S`: adiabatic bulk modulus of the solid phase
- `Ks_S`: shear modulus of the solid phase
- `ρf`  : density of the pore fluid
- `ρS`  : density of the solid phase
- `Vs0` : initial S-wave velocity of the solid phase
- `depth`: in kilometers
- `α`   : aspect ratio of the pores, modelled as oblate spheroids (0 < α < 1)

Output:
====
- `Vs_cor` : S-wave velocity corrected for fluid-filled porosity, clamped to be non-negative

The pores are treated like melt inclusions in [`melt_correction_Takei`](@ref): the shear-modulus
reduction of the self-consistent oblate-spheroid model (Dean, 1983; Phani, 1996) enters the
velocity reduction of Clark & Lesher (2017). For `α = 1` the inclusion model is undefined and `Vs0`
is returned unchanged.

References:
====

- Chen et al. (2020) Empirical porosity-depth model for continental crust, Hydrogeology Journal
- Clark & Lesher (2017) Elastic properties of silicate melts: Implications for low velocity zones at the lithosphere-asthenosphere boundary. Science Advances, Vol 3 (12), e1701312
- Dean (1983) Elastic moduli of porous sintered materials as modeled by a variable-aspect-ratio self-consistent oblate-spheroidal-inclusion theory, Journal of the American Ceramic Society, 66(12), 847--854
- Phani (1996) Porosity-dependence of ultrasonic velocity in sintered materials - a model based on the self-consistent spheroidal inclusion theory, Journal of Materials Science, 31, 262--266

"""
function porosity_correction(
        Kb_S::_T, Ks_S::_T, ρf::_T, ρS::_T, Vs0::_T, depth::_T, α::_T
    ) where {_T <: Number}
    ϕ = porosity_Chen2020(depth)
    Λ = inclusion_Λ(Kb_S, Ks_S, ϕ, α)
    isnothing(Λ) && return Vs0
    ΛG = Λ[2]
    ΔVs = (ΛG - (1 - ρf / ρS)) * ϕ / 2 * Vs0
    return max(Vs0 - ΔVs, zero(Vs0))
end

# Empirical porosity-depth model for continental crust (Chen et al., 2020); depth in km
@inline function porosity_Chen2020(depth::_T) where {_T}
    m = _T(0.071)
    n = _T(5.989)
    ϕ0 = _T(0.474)
    return ϕ0 / fastpow(1 + depth * m, n)
end

"""
        Vs_anel = anelastic_correction(water::Integer, Vs0, P, T)

This routine computes a correction of S-wave velocity for anelasticity

Input:
====
- `water`: water flag, 0 = dry; 1 = dampened; 2 = water saturated
- `Vs0`  : S-wave velocitiy of the solid phase (with or without melt correction)
- `Pref` : pressure given in Pa
- `Tref` : temperature given in °K

Output:
====
- `Vs_anel` : corrected S-wave velocity for anelasticity

The routine uses the reduction formulation of Karato (1993), using the quality factor formulation from Behn et al. (2009)


References:
====

- Karato, S. I. (1993). Importance of anelasticity in the interpretation of seismic tomography. Geophysical research letters, 20(15), 1623-1626.

- Behn, M. D., Hirth, G., & Elsenbeck II, J. R. (2009). Implications of grain size evolution on the seismic structure of the oceanic upper mantle. Earth and Planetary Science Letters, 282(1-4), 178-189.


"""
function anelastic_correction(water::Integer, Vs0, Pref, Tref)
    Tc = precision_of(Vs0)
    Pref, Tref = convert_precision(Tc, Pref), convert_precision(Tc, Tref)
    R = Tc(8.31446261815324)     # gas constant

    # values based on fitting experimental constraints (Behn et al., 2009)
    α = Tc(0.27)
    B0 = Tc(1.28e8)          # m/s
    dref = Tc(1.24e-5)       # m
    COHref = Tc(50.0e-6)     # 50H/1e6Si

    Gref = Tc(1.09)
    Eref = Tc(505.0e3) # J/mol
    Vref = Tc(1.2e-5)  # m3*mol

    G = 1
    E = Tc(420.0e3)              # J/mol (activation energy)
    V = Tc(1.2e-5)               # m3*mol (activation volume)

    # using remaining values from Cobden et al., 2018
    ω = Tc(0.01)                 # Hz (frequency to match for studied seismic system)
    d = Tc(1.0e-2)               # m (grain size)

    if water == 0
        COH = Tc(50.0e-6)     # for dry mantle
        r = 0               # for dry mantle
    elseif water == 1
        COH = Tc(1000.0e-6)   # for damp mantle
        r = 1               # for damp mantle
    elseif water == 2
        COH = Tc(3000.0e-6)   # for wet mantle (saturated water)
        r = 2               # for wet mantle
    else
        throw(
            ArgumentError(
                "water mode $water is not implemented. Valid values are 0 (dry), 1 (dampened) and 2 (wet)",
            ),
        )
    end

    B = @muladd @pow B0 *
        dref^(G - Gref) *
        (COH / COHref)^r *
        exp(((Pref * V + E) - (Pref * Vref + Eref)) / (R * Tref))

    Qinv = @pow (B * d^(-G) * inv(ω) * exp(-(Pref * V + E) / (R * Tref)))^α

    Vs_anel = Vs0 * (1 - Qinv / (2 * tan(π * α / 2)))

    return Vs_anel
end

"""
    PD_corrected = correct_wavevelocities_phasediagrams(PD::PhaseDiagram_LookupTable;
                                apply_porosity_correction=true, ρf=1000.0, α_porosity=0.1,
                                apply_melt_correction=true, α_melt=nothing, melt_correction_takei=true,
                                apply_anelasticity_correction=false, water=0,
                                combine=:sequential)

This applies various corrections to the seismic velocities specified in the phase diagram lookup table `PD`, and returns a new lookup table in which the fields `Vp`, `Vs` and `VpVs` hold the corrected values.
The original `Vp`,`Vs` is stored in `Vp_uncorrected`,`Vs_uncorrected`; `PD` itself is not modified.

The following corrections can be applied (together with potential options)
- *apply_anelasticity_correction*: applies an anelasticity correction ([`anelastic_correction`](@ref)) to the S-wave velocity of the solid, with the optional parameter `water`: 0 = dry; 1 = dampened; 2 = water saturated
- *apply_porosity_correction*: applies a correction for fluid-filled pores ([`porosity_correction`](@ref)) to the S-wave velocity. Optional parameters are `ρf` (fluid density, [kg/m3]) and `α_porosity` (pore aspect ratio)
- *apply_melt_correction*: applies a correction to the P- and S-wave velocities for the melt fraction of the diagram. `melt_correction_takei=true` uses [`melt_correction_Takei`](@ref), where `α_melt` is the aspect ratio of the melt inclusions (default 0.1); `false` uses [`melt_correction`](@ref), where `α_melt` is the contiguity of the solid framework (default 0.84, textural equilibrium at about 1% melt)

`combine` selects how the porosity and melt corrections act together:
- `:sequential`: the melt correction is applied to the porosity-corrected velocities
- `:weighted`: both corrections start from the velocities of the solid, and the results are averaged with the porosity and the melt fraction as weights
"""
function correct_wavevelocities_phasediagrams(
        PD::PhaseDiagram_LookupTable;
        apply_porosity_correction = true,
        ρf = 1000.0,
        α_porosity = 0.1,
        apply_melt_correction = true,
        α_melt = nothing,
        melt_correction_takei = true,
        apply_anelasticity_correction = false,
        water = 0,
        combine = :sequential,
    )
    α_melt = something(α_melt, melt_correction_takei ? 0.1 : 0.84)
    combine in (:sequential, :weighted) ||
        throw(ArgumentError("`combine` must be :sequential or :weighted, got $(repr(combine))"))

    # reconstruct the T,P knot vectors from the (regular) interpolation grid
    grid = PD.solid_Vs
    T = range(grid.T0, grid.Tmax; length = grid.numT)
    P = range(grid.P0, grid.Pmax; length = grid.numP)

    Vs_corrected = copy(PD.solid_Vs.coefs)
    Vp_corrected = copy(PD.solid_Vp.coefs)

    Kb_S = PD.solid_bulkModulus.coefs
    Ks_S = PD.solid_shearModulus.coefs
    ρS = PD.rockRho.coefs
    ρ_av = mean(ρS)
    if apply_melt_correction
        Kb_L = PD.melt_bulkModulus.coefs
        ρL = PD.meltRho.coefs
        ϕ_melt = PD.meltFrac.coefs
    end

    for I in CartesianIndices(Vs_corrected)
        Vp0 = Vp_corrected[I]
        Vs0 = Vs_corrected[I]
        if apply_anelasticity_correction
            Vs0 = anelastic_correction(water, Vs0, P[I[2]], T[I[1]])
        end

        # porosity correction (S-wave only)
        ϕW, VsW = zero(Vs0), Vs0
        if apply_porosity_correction
            depth = P[I[2]] / (9.81 * ρ_av * 1.0e3)    # approximate depth in km (lithostatic P)
            ϕW = porosity_Chen2020(depth)
            VsW = porosity_correction(Kb_S[I], Ks_S[I], ρf, ρS[I], Vs0, depth, α_porosity)
        end

        # melt correction, applied to the porosity-corrected or to the solid velocities
        ϕM, VpM, VsM = zero(Vs0), Vp0, VsW
        if apply_melt_correction && ϕ_melt[I] > 0
            ϕM = ϕ_melt[I]
            Vs_in = combine === :sequential ? VsW : Vs0
            if melt_correction_takei
                VsM, VpM = melt_correction_Takei(Kb_L[I], Kb_S[I], Ks_S[I], ρL[I], ρS[I], Vp0, Vs_in, ϕM, α_melt)
            else
                VpM, VsM = melt_correction(Kb_L[I], Kb_S[I], Ks_S[I], ρL[I], ρS[I], Vp0, Vs_in, ϕM, α_melt)
            end
            VpM, VsM = max(VpM, zero(VpM)), max(VsM, zero(VsM))
        end

        if combine === :sequential || ϕM == 0
            Vp_corrected[I], Vs_corrected[I] = VpM, VsM
        else
            ϕ_tot = ϕW + ϕM
            Vp_corrected[I] = (ϕW * Vp0 + ϕM * VpM) / ϕ_tot
            Vs_corrected[I] = (ϕW * VsW + ϕM * VsM) / ϕ_tot
        end
    end

    # Store results ----

    # Create interpolation objects (reusing the grid spacing of the input diagram)
    gp = (grid.T0, grid.dT, grid.numT, grid.Tmax, grid.P0, grid.dP, grid.numP, grid.Pmax)
    Vs_corrected_intp = interpolate(gp..., Vs_corrected)
    Vp_corrected_intp = interpolate(gp..., Vp_corrected)
    VpVs_corrected_intp = interpolate(gp..., Vp_corrected ./ Vs_corrected)

    # Initialize fields in the order they are defined in the PhaseDiagram_LookupTable structure
    Struct_Fieldnames = fieldnames(PhaseDiagram_LookupTable)[3:end] # fieldnames from structure

    # Process all fields that are present in the phase diagram (and non-dimensionalize if requested)
    Struct_Fields = Vector{Union{Nothing, LinearInterpolator}}(
        nothing, length(Struct_Fieldnames)
    )

    # Loop through all fields & copy the existing field. For :Vs_corrected, :Vp_corrected, we create the new objects
    for (i, field) in enumerate(Struct_Fieldnames)
        data = getfield(PD, field)

        if data != nothing
            Struct_Fields[i] = data
        end

        # add corrected Vs/Vp fields
        if field == :Vs
            Struct_Fields[i] = Vs_corrected_intp
        end
        if field == :Vp
            Struct_Fields[i] = Vp_corrected_intp
        end
        if field == :VpVs
            Struct_Fields[i] = VpVs_corrected_intp
        end

        # Store
        if field == :Vs_uncorrected
            Struct_Fields[i] = PD.Vs
        end
        if field == :Vp_uncorrected
            Struct_Fields[i] = PD.Vp
        end
    end

    # Store in phase diagram structure
    PD_corrected = PhaseDiagram_LookupTable(
        PD.Type, PD.Name, Struct_Fields...
    )

    return PD_corrected
end

"""
    Vs,Vp = melt_correction_Takei(Kb_L, Kb_S, Ks_S, ρL, ρS, Vp0, Vs0, ϕ, α)

This corrects the Vp/Vs velocities if melt is present. The melt sits in oblate-spheroidal inclusions
whose effect on the bulk and shear moduli follows the self-consistent model of Dean (1983) and
Phani (1996); the resulting velocity reduction follows Clark & Lesher (2017). The same parameterization
is used in Takei (2002) to relate Vp/Vs to pore geometry.

Input arguments:
- Kb_L: bulk modulus liquid
- Kb_S: bulk modulus solid
- Ks_S: shear modulus solid
- ρL: density liquid
- ρS: density solid
- Vp0: P-wave velocity of solid
- Vs0: S-wave velocity of solid
- ϕ: melt content
- α: melt aspect ratio [0.001 - 1)

Output arguments:
- Vs: corrected S-wave velocity
- Vp: corrected P-wave velocity

The reduction is first order in `ϕ`; corrected velocities that would become negative are clamped to
zero. For `α = 1` the inclusion model is undefined and `Vs0`, `Vp0` are returned unchanged.
"""
function melt_correction_Takei(
        Kb_L::_T, Kb_S::_T, Ks_S::_T, ρL::_T, ρS::_T, Vp0::_T, Vs0::_T, ϕ::_T, α::_T
    ) where {_T <: Number}
    Λ = inclusion_Λ(Kb_S, Ks_S, ϕ, α)
    isnothing(Λ) && return Vs0, Vp0
    ΛK, ΛG = Λ

    # Seismic wave velocity melt correction Clark & Lesher, 2017
    β = Kb_S / Kb_L
    γ = Ks_S / Kb_S
    Δρ = 1 - ρL / ρS

    # Formulation of the fraction reduction of P-wave and S-wave
    ΔVp = ((((β - 1) * ΛK) / ((β - 1) + ΛK) + 4γ * ΛG / 3) / (1 + 4γ / 3) - Δρ) * ϕ / 2 * Vp0
    ΔVs = (ΛG - Δρ) * ϕ / 2 * Vs0

    return max(Vs0 - ΔVs, zero(Vs0)), max(Vp0 - ΔVp, zero(Vp0))
end

# Geometric factors (ΛK, ΛG) = (P0, Q0) of Clark & Lesher (2017) for inclusions of aspect ratio α
# at volume fraction ϕ, from the self-consistent oblate-spheroid model; `nothing` where the model
# is undefined (α = 1).
function inclusion_Λ(Kb_S::_T, Ks_S::_T, ϕ::_T, α::_T) where {_T}
    f(x) = R_func(x, α, ϕ, Kb_S, Ks_S)
    isnan(f(_T(0.5))) && return nothing

    # The bracket can contain two roots, in which case the bisection is not applicable.
    lo, hi = _T(1.0e-3), 1 - _T(1.0e-3)
    R = if sign(f(lo)) != sign(f(hi))
        find_zero(f, (lo, hi), Bisection())
    else
        fzero(f, _T(0.5))
    end
    return P0_func(α, R), Q0_func(α, R)
end

"""
- Dean (1983), Elastic Moduli of Porous Sintered Materials as Modeled by a
Variable-Aspect-Ratio Self-Consistent Oblate-Spheroidal-Inclusion Theory

- Phani (1996), Porosity-dependence of ultrasonic velocity in sintered
materials - a model based on the self-consistent spheroidal inclusion theory
"""
function θ_func(α)
    s = sqrt(1 - α^2)
    return α / s^3 * (acos(α) - α * s)
end

"""
- Dean (1983), Elastic Moduli of Porous Sintered Materials as Modeled by a
Variable-Aspect-Ratio Self-Consistent Oblate-Spheroidal-Inclusion Theory

- Phani (1996), Porosity-dependence of ultrasonic velocity in sintered
materials - a model based on the self-consistent spheroidal inclusion theory
"""
function f_func(α)
    θ = θ_func(α)
    return α^2 * (3 * θ - 2) / (1 - α^2)
end

function P0_func(α, R)
    f = f_func(α)
    θ = θ_func(α)
    F1 = 1 - 3 * (f + θ) / 2 + R * (3 * f / 2 + 5 * θ / 2 - 4 // 3)
    F2 = R * (2 * θ - 2 * f - 3 * θ^2 + 2R * (f - θ + 2 * θ^2))
    return F1 / F2
end

function Q0_func(α, R)
    f = f_func(α)
    θ = θ_func(α)

    F2 = R * (2 * θ - 2 * f - 3 * θ^2 + 2 * R * (f - θ + 2 * θ^2))
    F3 = f + 3 * θ / 2 - R * (f + θ)
    F4 = 1 - (f + 3 * θ - R * (f - θ)) / 4
    F5 = f - R * (f + θ - 4 // 3)
    F6 = -f + R * (f + θ)
    F7 = 2 - (3 * f + 9 * θ - R * (3 * f + 5 * θ)) / 4
    F8 = -1 + f / 2 + 3 * θ / 2 + R * (2 - f / 2 - 5 * θ / 2)
    F9 = f - R * (f - θ)

    return (2 / F3 + 1 / F4 + F5 / F2 + (F6 * F7 - F8 * F9) / (F2 * F4)) / 5
end

# Self-consistency condition for R = 3G*/(3K* + 4G*), with K* = K_m(1 - Φ P0) and
# G* = G_m(1 - Φ Q0) the moduli of the porous solid.
function R_func(R::_T, α::_T, Φ::_T, K_m::_T, G_m::_T) where {_T}
    P0 = P0_func(α, R)
    Q0 = Q0_func(α, R)
    p1 = R * (3 * K_m + 4 * G_m)
    p2 = -3 * G_m
    p3 = -Φ * (R * (3 * K_m * P0 + 4 * G_m * Q0) - 3 * G_m * Q0)
    return p1 + p2 + p3
end

end
