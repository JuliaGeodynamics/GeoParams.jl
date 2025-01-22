@with_kw_noshow struct RSF{T, U1, U2, U3, U4} <: AbstractRateAndState
    γ0::GeoUnit{T, U1} = 0NoUnits # dynamic internal friction coefficient
    Θ::GeoUnit{T, U1} = 0NoUnits
    a::GeoUnit{T, U1} = 0NoUnits
    b::GeoUnit{T, U1} = 0NoUnits
    V0::GeoUnit{T, U2} = 0m / s
    L::GeoUnit{T, U3} = 0m # characteristic slip distance
    μ::GeoUnit{T, U4} = 1.0e9Pa # shear modulus
    ν::GeoUnit{T, U1} = 0.5NoUnits # poisson ratio
end

@inline function compoute_friction(x::RSF, V)
    (; γ0, a, b, V0, L, Θ) = x
    return γ0 + a * log(V / V0) + b * log(Θ * V0 / L)
end

@inline compute_shear_stress(x::RSF, V, σn) = compoute_friction(x, V) * σn

@inline compute_dΘdt(x::RSF, V) = 1 - V * x.Θ / x.L

##########

@with_kw_noshow struct RSF_regularized{T, U1, U2, U3, U4} <: AbstractRateAndState
    γ0::GeoUnit{T, U1} = 0NoUnits # dynamic internal friction coefficient
    Θ::GeoUnit{T, U1} = 0NoUnits
    a::GeoUnit{T, U1} = 0NoUnits
    b::GeoUnit{T, U1} = 0NoUnits
    V0::GeoUnit{T, U2} = 0m / s
    L::GeoUnit{T, U3} = 0m # characteristic slip distance
    μ::GeoUnit{T, U4} = 1.0e9Pa # shear modulus
    ν::GeoUnit{T, U1} = 0.5NoUnits # poisson ratio
end

@inline function compute_shear_stress(x::RSF_regularized, V, σn)
    (; γ0, a, b, V0, L, Θ) = x
    return asinh(0.5 * (V / V0) * exp((b * log(Θ * V0 / L) + γ0) / a)) * a * σn
end

##########

@with_kw_noshow struct RSF_invariant_regularized{T, U1, U2, U3, U4, U5} <: AbstractRateAndState
    γ0::GeoUnit{T, U1} = 0NoUnits # dynamic internal friction coefficient
    Θ::GeoUnit{T, U1} = 0NoUnits
    a::GeoUnit{T, U1} = 0NoUnits
    b::GeoUnit{T, U1} = 0NoUnits
    V0::GeoUnit{T, U2} = 0m / s
    L::GeoUnit{T, U3} = 0m # characteristic slip distance
    σc::GeoUnit{T, U4} = 0Pa # compressive strength
    μ::GeoUnit{T, U5} = 1.0e9Pa # shear modulus
    ν::GeoUnit{T, U1} = 0.5NoUnits # poisson ratio
end

@inline function compute_yield_stress(x::RSF_invariant_regularized, Vp, P)
    (; σc, γ0, a, b, V0, L, Θ) = x
    Ω = log(Θ * V0) / L
    return σc + asinh(0.5 * (Vp / V0) * exp((b * Ω + γ0) / a)) * a * P
end

@inline compute_Ω(x::AbstractRateAndState) = log(x.Θ * V0) / x.L
@inline compute_dΩdt(x::RSF_invariant_regularized, Vp) = (x.V0 * exp(-compute_Ω(x)) - Vp) / x.L

## Adaptive time stepping

function compute_Δt(x::AbstractRateAndState, Vp, ζ, Δt)
    @assert ζ < 1
    return ζ * min(
        compute_Δt_weakening(x, Vp, μ, Δx),
        compute_Δt_healing(x, Vp),
        Δt,
    )
end

# weakening time step
function compute_Δt_weakening(x::AbstractRateAndState, Vp, μ, Δx)
    X_state = compute_X_state(x, P, μ, Δx)
    return X_state * x.L / Vp
end

function compute_X_state(x, P, μ, Δx)
    K = compute_K(x, μ, Δx)
    ε = compute_ε(x, K, P)
    return X_state = if ε < 0
        min(0.2, 1 - (x.b - x.a) * P / (K * L))
    else
        min(0.2, (x.a * P) / (K * L - (x.a - x.b) * P))
    end
end

@inline function compute_ε(x, K, P)
    (; a, b, L) = x
    return 0.25 * (K * L / (a * P) - (b - a) / a)^2 - K * L / (a * P)
end

@inline compute_K(x, μ, Δx) = 2 * μ * (π * Δx * (1 - x.ν))

# healing time step
@inline compute_Δt_healing(x::AbstractRateAndState, Vp) = X_state * min(x.L / Vp, x.L / V0 * compute_Ω(x))
