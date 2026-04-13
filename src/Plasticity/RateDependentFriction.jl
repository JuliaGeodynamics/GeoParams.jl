export RateDependentFriction,
    comput_Vp,
    compute_gamma_effective,
    compute_plastic_strain_rate,
    compute_yield_stress

# RateDependentFriction  -------------------------------------------------------
"""
    RateDependentFriction()   
"""
@with_kw_noshow struct RateDependentFriction{T, U1, U2, U3, U4} <: AbstractPlasticity{T}
    σc::GeoUnit{T, U1} = 10.0e6Pa     # compressive strength
    γd::GeoUnit{T, U2} = 0.3NoUnits # dynamic internal friction coefficient
    γs::GeoUnit{T, U2} = 0.3NoUnits # static internal friction coefficient
    D::GeoUnit{T, U3} = 0m
    Vc::GeoUnit{T, U4} = 1.0e-8m / s     # characteristic slip rate at which half
    # of the friction change occurs
end

RateDependentFriction(args::Vararg{Any, N}) = RateDependentFriction(convert.(GeoUnit, args)...)

@inline function compute_yield_stress(x::RateDependentFriction, η, ηvp, P, τII)
    εII_plastic = compute_plastic_strain_rate(x, τII, η, ηvp)
    Vp = comput_Vp(x, εII_plastic)
    γ_effective = compute_gamma_effective(x, Vp)
    return σc + γ_effective * P
end

@inline comput_Vp(x::RateDependentFriction, εII_plastic) = 2 * εII_plastic * x.D

@inline compute_gamma_effective(x::RateDependentFriction, Vp) = x.γd + (x.γs - x.γd) / (1 + Vp / x.Vc)

@inline compute_plastic_strain_rate(::RateDependentFriction, τII, η, ηvp) = 0.5 * τII * (inv(ηvp) - inv(η))
