# export RateDependentFriction,
#     comput_Vp,
#     compute_gamma_effective,
#     compute_plastic_strain_rate,
#     compute_yield_stress

# RateDependentFriction  -------------------------------------------------------
"""
    RateDependentFriction()   
"""
struct LowTemperaturePlasticity{T, U1, U2, U3, U4, U6, U7, U8} <: AbstractPlasticity{T}
    σres::GeoUnit{T, U1}
    σL::GeoUnit{T, U2}
    σK::GeoUnit{T, U3}
    σb::GeoUnit{T, U4}
    E::GeoUnit{T, U5}
    R::GeoUnit{T, U6}
    V::GeoUnit{T, U7}
    A::GeoUnit{T, U8}
end

RateDependentFriction(args::Vararg{Any, N}) = RateDependentFriction(convert.(GeoUnit, args)...)

@inline function compute_yield_stress(x::RateDependentFriction, η, ηvp, P, τII)
    εII_plastic = compute_plastic_strain_rate(x, τII, η, ηvp)
    Vp = comput_Vp(x, εII_plastic)
    γ_effective = compute_gamma_effective(x, Vp)
    return σc + γ_effective * P
end
