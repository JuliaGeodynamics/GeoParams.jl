module SeismicVelocity

# This implements different methods to compute seismic velocities 
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using ..PhaseDiagrams
using ..MaterialParameters: MaterialParamsInfo
using Interpolations, Statistics
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct

using Roots
import Base.show, GeoParams.param_info

abstract type AbstractSeismicVelocity{T} <: AbstractMaterialParam end

export compute_pwave_velocity,
    compute_wave_velocity, # calculation routines
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
@with_kw_noshow struct ConstantSeismicVelocity{T,U} <: AbstractSeismicVelocity{T}
    Vp::GeoUnit{T,U} = 8.1km / s               # P-wave velocity
    Vs::GeoUnit{T,U} = 4.5km / s               # S-wave velocity
end
ConstantSeismicVelocity(args...) = ConstantSeismicVelocity(convert.(GeoUnit, args)...)

function param_info(s::ConstantSeismicVelocity) # info about the struct
    return MaterialParamsInfo(; Equation=L"v_p = cst \\ v_s = cst")
end

# Calculation routines
function compute_wave_velocity(s::ConstantSeismicVelocity{_T}; wave, kwargs...) where {_T}
    wave != :VpVs && return getfield(s, wave)
    @unpack Vp, Vs = s
    return Vp / Vs
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
    compute_pwave_velocity(s::PhaseDiagram_LookupTable, P,T)
Interpolates Vp, Vs or VpVs velocity as a function of `T,P` from a lookup table  
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

compute_wave_velocity!(args...) = compute_param!(compute_wave_velocity, args...)
compute_wave_velocity(args...) = compute_param(compute_wave_velocity, args...)

#=
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
- `α`   : contiguity coefficient defining the geometry of the solid framework (contiguity)
          0.0 (layered melt distributed) < 0.1 (grain boundary melt) < 1.0 (melt in separated bubble pockets)

Output:
====
- `Vp_cor,Vs_cor` : corrected P-wave and S-wave velocities for melt fraction

The routine uses the reduction formulation of Clark & Lesher, (2017) and is based on the equilibrium geometry model for the solid skeleton of Takei et al., 1998.

References:
====

- Takei (1998) Constitutive mechanical relations of solid-liquid composites in terms of grain-boundary contiguity, Journal of Geophysical Research: Solid Earth, Vol(103)(B8), 18183--18203

- Clark & Lesher (2017) Elastic properties of silicate melts: Implications for low velocity zones at the lithosphere-asthenosphere boundary. Science Advances, Vol 3 (12), e1701312

"""
function melt_correction(
    Kb_L::_T, Kb_S::_T, Ks_S::_T, ρL::_T, ρS::_T, Vp0::_T, Vs0::_T, ϕ::_T, α::_T
) where {_T<:Number}

    # Takei 1998: Approximation Formulae for Bulk and Shear Moduli of Isotropic Solid Skeleton
    ν = 0.25                         # poisson ratio

    aij = (
        0.318, 6.780, 57.560,  0.182,
        0.164, 4.290, 26.658,  0.464,
        1.549, 4.814, 8.777, -0.290,
    )
    bij = (
        -0.3238, 0.2341, 
        -0.1819, 0.5103
    )

    # Lines below are equivalent to:
    # a = zeros(3)
    # for i in 1:3
    #     a[i] =
    #         aij[i, 1] * exp(aij[i, 2] * (ν - 0.25) + aij[i, 3] * (ν - 0.25)^3) + aij[i, 4]
    # end
    a = ntuple(Val(3)) do i
        idx = 4*i-3 # linear offset index
        aij[idx] * exp(aij[idx+1] * (ν - 0.25) + aij[idx+2] * (ν - 0.25)^3) + aij[idx+3]
    end

    # Lines below are equivalent to:
    # b = zeros(2)
    # for i in 1:2
    #     b[i] = bij[i, 1] * ν + bij[i, 2]
    # end
    b = ntuple(Val(2)) do i
        idx = 2*i-1 # linear offset index
        bij[idx] * ν + bij[idx+1]
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

    # ratio of skeleton adiabatic bulk modulus over solid phase bulk modulus
    ΛK = Kb_S / kb

    # ratio of skeleton shear modulus over solid phase shear modulus
    ΛG = Ks_S / μ

    # Seismic wave velocity melt correction Clark et al., 2017
    β = Kb_S / Kb_L
    γ = Ks_S / Kb_S

    # Formulation of the fraction reduction of P-wave and S-wave
    ΔVp =
        (
            (
                (((β - 1.0) * ΛK) / ((β - 1.0) + ΛK) + 4.0 / 3.0 * γ * ΛG) /
                (1.0 + 4.0 / 3.0 * γ)
            ) - (1.0 - ρL / ρL)
        ) * (ϕ * 0.5)
    ΔVs = (ΛG - (1.0 - ρL / ρS)) * (ϕ * 0.5)

    # get the correction values
    Vp_cor = Vp0 - ΔVp * Vp0
    Vs_cor = Vs0 - Vs0 * ΔVs

    return Vp_cor, Vs_cor
end
=#

#=
"""
        Vs_cor = porosity_correction(  Kb_L, Kb_S, Ks_S, ρL, ρS, Vp0, Vs0, ϕ, α)

Corrects S-wave velocity at shallow depth as function of empirical porosity-depth profile. 

Input:
====
- `Kb_S`: adiabatic bulk modulus of the solid phase
- `Ks_S`: shear modulus of the solid phase
- `ρL`  : density of the melt
- `ρS`  : density of the solid phase
- `Vs0` : initial S-wave velocity of the solid phase
- `depth`: in kilometers
- `α`   : contiguity coefficient defining the geometry of the solid framework (contiguity)
          0.0 (layered fluid distributed) < 0.1 (grain boundary melt) < 1.0 (fluid in separated bubble pockets)

Output:
====
- `Vs_cor` : corrected P-wave and S-wave velocities for water-filled porosity

The routine is based on the equilibrium geometry model for the solid skeleton of Takei et al., 1998.

References:
====

- Takei (1998) Constitutive mechanical relations of solid-liquid composites in terms of grain-boundary contiguity, Journal of Geophysical Research: Solid Earth, Vol(103)(B8), 18183--18203
- Chen et al. (2020) Empirical porosity-depth model for continental crust

"""
function porosity_correction(
    Kb_S::_T, Ks_S::_T, ρf::_T, ρS::_T, Vs0::_T, depth::_T, α::_T
) where {_T<:Number}
    # Empirical porosity-depth model for continental crust after Chen et al., 2020 (hydrogeology journal)
    m = 0.071
    n = 5.989
    ϕ0 = 0.474

    ϕ = ϕ0 / fastpow((1.0+depth*m), n)

    # Takei 1998: Approximation Formulae for Bulk and Shear Moduli of Isotropic Solid Skeleton
    ν = 0.25 # poisson ratio

    # Tuples are a better option than standard arrays for matrices/vectors of known size at compile time

    #=
    aij = (
        0.318, 6.780, 57.560,  0.182,
        0.164, 4.290, 26.658,  0.464,
        1.549, 4.814, 8.777, -0.290,
    )
    bij = (
        -0.3238, 0.2341, 
        -0.1819, 0.5103
    )
    =#

    # Takei et al. 2002:
    aij = (
        1.8625,  0.52594, -4.8397,  0.,
        4.5001, -6.1551,  -4.3634,  0.0,
       -5.6512, 6.9159,   29.595, -58.96
    )
    bij = (
        1.6122, 0.13527, 0.0,
        4.5869, 3.6086, 0.0,
        -7.5395,-4.8676, -4.3182
    )

    # Lines below are equivalent to:
    a = zeros(3)
    #for i in 1:4
    #     a[i] =
    #         aij[i, 1] * exp(aij[i, 2] * (ν - 0.25) + aij[i, 3] * (ν - 0.25)^3) + aij[i, 4]
    #end
    #a = ntuple(Val(3)) do i
    #    idx = 4*i-3 # linear offset index
    #    aij[idx] * exp(aij[idx+1] * (ν - 0.25) + aij[idx+2] * (ν - 0.25)^3) + aij[idx+3]
    #end

    # Takei (2002): 
    a = ntuple(Val(3)) do i
        idx = 4*i-3 # linear offset index
        aij[idx] + aij[idx+1]*ν^1 + aij[idx+2]*ν^2 + aij[idx+3]*ν^3
    end

    # Lines below are equivalent to:
    # b = zeros(2)
    # for i in 1:2
    #     b[i] = bij[i, 1] * ν + bij[i, 2]
    # end
    b = ntuple(Val(2)) do i
        idx = 2*i-1 # linear offset index
        bij[idx] * ν + bij[idx+1]
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

    # ratio of skeleton adiabatic bulk modulus over solid phase bulk modulus
    ΛK = Kb_S / kb

    # ratio of skeleton shear modulus over solid phase shear modulus
    ΛG = Ks_S / μ

    ΔVs = (ΛG - (1.0 - ρf / ρS)) * (ϕ * 0.5)

    # get the correction values
    Vs_cor = Vs0 - Vs0 * ΔVs

    return Vs_cor
end
=#

"""
        Vs_anel = anelastic_correction(water::Int64, Vs0::Float64,P::Float64,T::Float64)

This routine computes a correction of S-wave velocity for anelasticity

Input:
====
- `water`: water flag, 0 = dry; 1 = dampened; 2 = water saturated 
- `Vs0`  : S-wave velocitiy of the solid phase (with or without melt correction)
- `P`    : pressure given in Pa
- `T`    : temperature given in °K

Output:
====
- `Vs_anel` : corrected S-wave velocity for anelasticity

The routine uses the reduction formulation of Karato (1993), using the quality factor formulation from Behn et al. (2009)


References:
====

- Karato, S. I. (1993). Importance of anelasticity in the interpretation of seismic tomography. Geophysical research letters, 20(15), 1623-1626.

- Behn, M. D., Hirth, G., & Elsenbeck II, J. R. (2009). Implications of grain size evolution on the seismic structure of the oceanic upper mantle. Earth and Planetary Science Letters, 282(1-4), 178-189.


"""
function anelastic_correction(water::Int64, Vs0::Float64, Pref::Float64, Tref::Float64)
    R = 8.31446261815324     # gas constant

    # values based on fitting experimental constraints (Behn et al., 2009)
    α = 0.27
    B0 = 1.28e8               # m/s
    dref = 1.24e-5              # m
    COHref = 50.0 / 1e6             # 50H/1e6Si

    Gref = 1.09
    Eref = 505.0e3              # J/mol
    Vref = 1.2e-5               # m3*mol

    G = 1.00
    E = 420.0e3              # J/mol (activation energy)
    V = 1.2e-5               # m3*mol (activation volume)

    # using remaining values from Cobden et al., 2018
    ω = 0.01                 # Hz (frequency to match for studied seismic system)
    d = 1e-2                 # m (grain size)

    if water == 0
        COH = 50.0 / 1e6     # for dry mantle
        r = 0.0              # for dry mantle
    elseif water == 1
        COH = 1000.0 / 1e6   # for damp mantle    
        r = 1.0              # for damp mantle
    elseif water == 2
        COH = 3000.0 / 1e6   # for wet mantle (saturated water)
        r = 2.0              # for wet mantle
    else
        print(
            "water mode is not implemented. Valid values are 0 (dry),1 (dampened) and 2 (wet)",
        )
    end

    B =
        B0 *
        dref^(G - Gref) *
        (COH / COHref)^r *
        exp(((E + Pref * V) - (Eref + Pref * Vref)) / (R * Tref))

    Qinv = (B * d^(-G) * ω^(-1.0) * exp(-(E + Pref * V) / (R * Tref)))^α

    Vs_anel = Vs0 * (1.0 - (Qinv) / (2.0 * tan(π * α / 2.0)))

    return Vs_anel
end

"""
    PD_corrected = correct_wavevelocities_phasediagrams(PD::PhaseDiagram_LookupTable,  
                                apply_porosity_correction=true, ρf=1000.0, α_porosity=0.5,
                                apply_melt_correction=true, α_melt = 0.1,
                                apply_anelasticity_correction=true, water=2)

This applies various corrections to the seismic velocities specified in the phase diagram lookup table `PD`, and replaces the fields `Vp` and `Vs` in the diagram with the corrected ones.
The original `Vp`,`Vs` is stored in `Vp_uncorrected`,`Vs_uncorrected`

The following corrections can be applied (together with potential options)
- *apply_porosity_correction*: applies a correction for fluid-filled pores to vs velocity. Optional parameters are `ρf` (density fluid=[1000kg/m3]) and  `α_porosity` (contiguity coefficient defining the geometry of the solid framework, with 0.0 (layered fluid distributed) < 0.1 (grain boundary melt) < 1.0 (fluid in separated bubble pockets)
- *apply_melt_correction*: applies a correction to the P/S-wave velocity for the presence of melt for a given pore contiguity described by `α_melt`: 0.0 (layered fluid distributed) < 0.1 (grain boundary melt) < 1.0 (fluid in separated bubble pockets)
- *apply_anelasticity_correction*: applies an anelasticity correction to the S-wave velocity, with the optional parameter `water`: 0 = dry; 1 = dampened; 2 = water saturated

"""
function correct_wavevelocities_phasediagrams(
    PD::PhaseDiagram_LookupTable;
    apply_porosity_correction=true,
    ρf=1000.0,
    α_porosity=0.1,
    apply_melt_correction=true,
    α_melt=0.1,
    melt_correction_takei=true,
    apply_anelasticity_correction=true,
    water=0,
)

    # extract required data
    T, P = PD.Rho.itp.knots   # T,P vectors of diagrams

    # store original results correction   
    Vs_uncorrected = PD.Vs
    Vp_uncorrected = PD.Vp

    # Extract required data to apply corrections
    Vs_corrected = PD.solid_Vs.itp.coefs   # Vs velocity of solid rocks
    Vp_corrected = PD.solid_Vp.itp.coefs   # Vp velocity of solid rocks

    # Apply anelasticity correction
    if apply_anelasticity_correction == true
        for i in CartesianIndices(Vs_corrected)
            Vs_corrected[i] = anelastic_correction(water, Vs_corrected[i], P[i[2]], T[i[1]])
        end
    end

    # Apply porosity correction 
    if apply_porosity_correction == true
        Kb_S = PD.solid_bulkModulus.itp.coefs    #  bulk modulus solid
        Ks_S = PD.solid_shearModulus.itp.coefs   #  shear modulus solid
        ρS = PD.rockRho.itp.coefs
        ρ_av = mean(ρS)           # average solid density

        for i in CartesianIndices(Vs_corrected)
            depth = P[i[2]] / (9.81 * ρ_av * 1e3)         # approximate depth in km (assuming lithostatic P)
            Vs_corrected[i] = porosity_correction(
                Kb_S[i], Ks_S[i], ρf, ρS[i], Vs_corrected[i], depth, α_porosity
            )
        end
    end

    # Apply melt correction 
    if apply_melt_correction == true
        Kb_L = PD.melt_bulkModulus.itp.coefs    #  bulk modulus melt
        Kb_S = PD.solid_bulkModulus.itp.coefs   #  bulk modulus solid
        Ks_S = PD.solid_shearModulus.itp.coefs  #  shear modulus solid
        ρS = PD.rockRho.itp.coefs                 #  solid density
        ρL = PD.meltRho.itp.coefs                 #  melt density
        ϕ = PD.meltFrac.itp.coefs                #  melt fraction

        for i in CartesianIndices(Vs_corrected)
            if ϕ[i] > 0
                if melt_correction_takei == true
                    Vp_c, Vs_c = melt_correction_Takei(
                        Kb_L[i],
                        Kb_S[i],
                        Ks_S[i],
                        ρL[i],
                        ρS[i],
                        Vp_corrected[i],
                        Vs_corrected[i],
                        ϕ[i],
                        α_melt,
                    )
                else
                    Vp_c, Vs_c = melt_correction(
                        Kb_L[i],
                        Kb_S[i],
                        Ks_S[i],
                        ρL[i],
                        ρS[i],
                        Vp_corrected[i],
                        Vs_corrected[i],
                        ϕ[i],
                        α_melt,
                    )
                end

                if Vs_c < 0.0
                    Vs_c = 0.0
                end
                if Vp_c < 0.0
                    Vp_c = 0.0
                end

                Vp_corrected[i], Vs_corrected[i] = Vp_c, Vs_c
            end
        end
    end

    # Store results ----

    # Create interpolation objects
    Vs_corrected_intp = LinearInterpolation((T, P), Vs_corrected; extrapolation_bc=Flat())
    Vp_corrected_intp = LinearInterpolation((T, P), Vp_corrected; extrapolation_bc=Flat())
    VpVs_corrected_intp = LinearInterpolation(
        (T, P), Vp_corrected ./ Vs_corrected; extrapolation_bc=Flat()
    )

    # Initialize fields in the order they are defined in the PhaseDiagram_LookupTable structure 
    Struct_Fieldnames = fieldnames(PhaseDiagram_LookupTable)[4:end] # fieldnames from structure

    # Process all fields that are present in the phase diagram (and non-dimensionalize if requested)
    Struct_Fields = Vector{Union{Nothing,Interpolations.Extrapolation}}(
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
            Struct_Fields[i] = Vs_uncorrected
        end
        if field == :Vp_uncorrected
            Struct_Fields[i] = Vp_uncorrected
        end
    end

    # Store in phase diagram structure
    PD_corrected = PhaseDiagram_LookupTable(
        "Perple_X/MAGEMin/LaMEM", PD.HeaderText, PD.Name, Struct_Fields...
    )

    return PD_corrected
end

"""
    Vs,Vp = melt_correction_Takei(Kb_L, Kb_S, Ks_S, ρL, ρS, Vp0, Vs0, ϕ, α)

This corrects the Vp/Vs velocities if melt is present, following Takei (2002)

Input arguments:
- Kb_L: bulk modulus liquid
- Kb_S: bulk modulus solid
- Ks_S: shear modulus solid
- ρL: density liquid
- ρS: density solid
- Vp0: P-wave velocity of solid
- Vs0: S-wave velocity of solid
- ϕ: melt content
- α: melt aspect ratio [0.001 - 1]

Output arguments:
- Vs: corrected S-wave velocity
- Vp: corrected P-wave velocity
"""
function melt_correction_Takei(
    Kb_L::_T, Kb_S::_T, Ks_S::_T, ρL::_T, ρS::_T, Vp0::_T, Vs0::_T, ϕ::_T, α::_T
) where {_T<:Number}

    # compute R
    f(x) = R_func(x, α, ϕ, Kb_S, Ks_S)

    # Note: this bracketing algorithm sometimes fails, as there might be 2 roots
    eps = 1e-3
    if !isnan(f(0.5))
        R = zero(Kb_L)
        if (sign(f(eps)) != sign(f(1.0 - eps)))
            R = find_zero(f, (eps, 1.0 - eps), Bisection())
        else
            # if bisection fails, try fzero 
            R = fzero(f, 0.5)
        end

        # compute Q0 and P0
        ΛG = Q0_func(α, R)
        ΛK = P0_func(α, R)

        # Seismic wave velocity melt correction Clark et al., 2017
        β = Kb_S / Kb_L
        γ = Ks_S / Kb_S

        # Formulation of the fraction reduction of P-wave and S-wave
        ΔVp =
            (
                (
                    (((β - 1.0) * ΛK) / ((β - 1.0) + ΛK) + 4.0 / 3.0 * γ * ΛG) /
                    (1.0 + 4.0 / 3.0 * γ)
                ) - (1.0 - ρL / ρL)
            ) * Vp0

        ΔVs = (ΛG - (1.0 - ρL / ρS)) * (ϕ * 0.5) * Vs0
    else
        ΔVp = 0.0
        ΔVs = 0.0
    end

    Vs_new = Vs0 - ΔVs
    Vp_new = Vp0 - ΔVp

    # @show Vs_new, Vp_new, ΔVp, R
    if Vs_new < 0
        Vs_new = 0.0
    end

    if Vp_new < 0
        Vp_new = 0.0
    end

    return Vs_new, Vp_new
end

"""
- Dean (1983), Elastic Moduli of Porous Sintered Materials as Modeled by a 
Variable-Aspect-Ratio Self-Consistent Oblate-Spheroidal-Inclusion Theory

- Phani (1996), Porosity-dependence of ultrasonic velocity in sintered 
materials - a model based on the self-consistent spheroidal inclusion theory 
"""
function θ_func(α::_T) where {_T}
    return α / ((1 - α^2)^(3 / 2)) * (acos(α) - α * (1 - α^2)^0.5)
end

"""
- Dean (1983), Elastic Moduli of Porous Sintered Materials as Modeled by a 
Variable-Aspect-Ratio Self-Consistent Oblate-Spheroidal-Inclusion Theory

- Phani (1996), Porosity-dependence of ultrasonic velocity in sintered 
materials - a model based on the self-consistent spheroidal inclusion theory 
"""
function f_func(α::_T) where {_T}
    θ = θ_func(α)
    return α^2 * (3 * θ - 2) / (1 - α^2)
end

function P0_func(α::_T, R::_T) where {_T}
    f = f_func(α)
    θ = θ_func(α)
    F1 = 1 - 3 / 2 * (f + θ) + R * (3 / 2 * f + 5 / 2 * θ - 4 / 3)
    F2 = R * (2 * θ - 2 * f - 3 * θ^2 + 2R * (f - θ + 2θ^2))
    return F1 / F2
end

function Q0_func(α::_T, R::_T) where {_T}
    f = f_func(α)
    θ = θ_func(α)

    F2 = R * (2 * θ - 2 * f - 3 * θ^2 + 2 * R * (f - θ + 2 * θ^2))
    F3 = f + 3 / 2 * θ - R * (f + θ)
    F4 = 1 - 1 / 4 * (f + 3 * θ - R * (f - θ))
    F5 = f - R * (f + θ - 4 / 3)
    F6 = -f + R * (f + θ)
    F7 = 2 - 1 / 4 * (3 * f + 9 * θ - R * (3 * f + 5 * θ))
    F8 = -1 + 1 / 2 * f + 3 / 2 * θ + R * (2 - 1 / 2 * f - 5 / 2 * θ)
    F9 = f - R * (f - θ)

    return 1 / 5 * (2 / F3 + 1 / F4 + F5 / F2 + (F6 * F7 - F8 * F9) / (F2 * F4))
end

function P0_func_deriv(α::_T, R::_T) where {_T}
    """
    d/dR (P0(R)) : Evaluated through SymPy:
        
    from sympy import symbols, diff, Rational
    R, f, θ, α = symbols('R f θ α')
    F1 = 1 - Rational(3,2)*(f+θ) + R*(Rational(3,2)*f + Rational(5,2)*θ - Rational(4,3))
    F2 = R * ( 2*θ - 2*f - 3*θ^2 + 2*R*(f - θ + 2*θ^2) )
    F3 = f + Rational(3, 2)*θ - R*(f + θ)
    F4 = 1 - Rational(1, 4) * (f + 3*θ - R*(f - θ))
    F5 = f - R*(f + θ - Rational(4, 3))
    F6 = -f + R*(f + θ)
    F7 = 2 - Rational(1, 4) * (3*f + 9*θ - R*(3*f + 5*θ))
    F8 = -1 + Rational(1, 2)*f + Rational(3, 2)*θ + R*(2 - Rational(1, 2)*f -Rational(5, 2)*θ)
    F9 = f - R*(f - θ)
    P0 = F1/F2
    derivative = diff(P0, R)
    """
    f = f_func(α)
    θ = θ_func(α)
    p1 = -2 * f - 4 * θ^2 + 2 * θ
    p2 = R * (3 * f / 2 + 5 * θ / 2 - 4 / 3) - 3 * f / 2 - 3 * θ / 2 + 1
    p3 = R * (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ)^2
    p4 = 3 * f / 2 + 5 * θ / 2 - 4 / 3
    p5 = R * (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ)
    p6 = R * (3 * f / 2 + 5 * θ / 2 - 4 / 3) - 3 * f / 2 - 3 * θ / 2 + 1
    p7 = R^2 * (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ)
    return p1 * p2 / p3 + p4 / p5 - p6 / p7
end

function Q0_func_deriv(α::_T, R::_T) where {_T}
    """
    d/dR (Q0(R)) : Evaluated through SymPy
    """
    f = f_func(α)
    θ = θ_func(α)
    p1 = (-f / 4 + θ / 4) / (5 * (R * (f - θ) / 4 - f / 4 - 3 * θ / 4 + 1)^2)
    p2 = 2 * (f + θ) / (5 * (-R * (f + θ) + f + 3 * θ / 2)^2)
    p3 = -f / 4 + θ / 4
    p4 = (-R * (f - θ) + f) * (R * (-f / 2 - 5 * θ / 2 + 2) + f / 2 + 3 * θ / 2 - 1)
    p5 = (R * (f + θ) - f) * (R * (3 * f + 5 * θ) / 4 - 3 * f / 4 - 9 * θ / 4 + 2)
    p6 =
        5 *
        R *
        (R * (f - θ) / 4 - f / 4 - 3 * θ / 4 + 1)^2 *
        (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ)
    p7 = (-R * (f + θ - 4 / 3) + f) * (-2 * f - 4 * θ^2 + 2 * θ)
    p8 = 5 * R * (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ)^2
    p9 = (-R * (f - θ) + f) * (R * (-f / 2 - 5 * θ / 2 + 2) + f / 2 + 3 * θ / 2 - 1)
    p10 = (R * (f + θ) - f) * (R * (3 * f + 5 * θ) / 4 - 3 * f / 4 - 9 * θ / 4 + 2)
    p11 = -2 * f - 4 * θ^2 + 2 * θ
    p12 =
        5 *
        R *
        (R * (f - θ) / 4 - f / 4 - 3 * θ / 4 + 1) *
        (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ)^2
    p13 = (-f - θ + 4 / 3) / (5 * R * (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ))
    p14 = (3 * f / 4 + 5 * θ / 4) * (R * (f + θ) - f)
    p15 = (f - θ) * (R * (-f / 2 - 5 * θ / 2 + 2) + f / 2 + 3 * θ / 2 - 1)
    p16 = (f + θ) * (R * (3 * f + 5 * θ) / 4 - 3 * f / 4 - 9 * θ / 4 + 2)
    p17 = (R * (f - θ) - f) * (-f / 2 - 5 * θ / 2 + 2)
    p18 =
        5 *
        R *
        (R * (f - θ) / 4 - f / 4 - 3 * θ / 4 + 1) *
        (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ)
    p19 =
        (-R * (f + θ - 4 / 3) + f) /
        (5 * R^2 * (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ))
    p20 = (-R * (f - θ) + f) * (R * (-f / 2 - 5 * θ / 2 + 2) + f / 2 + 3 * θ / 2 - 1)
    p21 = (R * (f + θ) - f) * (R * (3 * f + 5 * θ) / 4 - 3 * f / 4 - 9 * θ / 4 + 2)
    p22 =
        5 *
        R^2 *
        (R * (f - θ) / 4 - f / 4 - 3 * θ / 4 + 1) *
        (2 * R * (f + 2 * θ^2 - θ) - 2 * f - 3 * θ^2 + 2 * θ)
    return p1 +
           p2 +
           p3 * (-p4 + p5) / p6 +
           p7 / p8 +
           (-p9 + p10) * p11 / p12 +
           p13 +
           (p14 + p15 + p16 + p17) / p18 - p19 - (-p20 + p21) / p22
end

function R_func(R::_T, α::_T, Φ::_T, K_m::_T, G_m::_T) where {_T}
    P0 = P0_func(α, R)
    Q0 = Q0_func(α, R)
    p1 = R * (3 * K_m + 4 * G_m)
    p2 = -3 * G_m
    p3 = -Φ * (R * (3 * K_m * P0 + 4 * G_m * Q0) - 3 * G_m * Q0)
    return p1 + p2 + p3
end

function R_func_deriv(R::_T, α::_T, Φ::_T, K_m::_T, G_m::_T) where {_T}
    P0 = P0_func(α, R)
    Q0 = Q0_func(α, R)
    P0_deriv = P0_func_deriv(α, R)
    Q0_deriv = Q0_func_deriv(α, R)
    p1 = 3 * K_m + 4 * G_m
    p2 = -3 * Φ * K_m * (P0 + R * P0_deriv)
    p3 = -4 * Φ * G_m * (Q0 + R * Q0_deriv)
    p4 = 3 * Φ * G_m * Q0_deriv
    return p1 + p2 + p3 + p4
end

"""
    simple Newton root finding algorithm as used to determine R
"""
function find_roots_R(R0::_T, α::_T, Φ::_T, K_m::_T, G_m::_T; tol=1e-5) where {_T}
    Rn = R0
    for n in 1:500
        f = R_func(Rn, α, Φ, K_m, G_m)
        if abs(f) < tol
            return Rn
        end
        df = R_func_deriv(Rn, α, Φ, K_m, G_m)

        if df == 0
            error("No solution found")
        end
        Rn = Rn - f / df
    end
end

end
