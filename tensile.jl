using ForwardDiff

struct TensilePlasticity{T}
    C::T
    ϕ::T
    sinϕ::T
    cosϕ::T
    Ψ::T
    sinΨ::T
    cosΨ::T
    PC1 ::T
    PC2 ::T
    τC1 ::T
    τC2 ::T
    β ::T
    ηve::T
    σT::T
    δ::T
    function TensilePlasticity(
        C::T, ϕ::T, Ψ::T, PC1::T, PC2::T, τC1::T, τC2::T, β::T, ηve::T, σT::T, δ::T
    ) where T
        sinϕ, cosϕ = sincos(ϕ)
        sinΨ, cosΨ = sincos(Ψ)
        new{T}(C, ϕ, sinϕ, cosϕ, Ψ, sinΨ, cosΨ, PC1, PC2, τC1, τC2, β, ηve, σT, δ)
    end
end

let # let block so that we do not pollute the namespace
    fn_Q = :QI, :QII, :QIII, :QIV, :QV
    fn_check_Q = :check_QI, :check_QII, :check_QIII, :check_QIV
    fn_pl = :plastic_strain_rateI, :plastic_strain_rateII, :plastic_strain_rateIII, :plastic_strain_rateIV, :plastic_strain_rateV
    fn = tuple(fn_Q..., fn_check_Q..., fn_pl...)
    for fnᵢ in fn
        @eval begin
            @inline $fnᵢ(args::NamedTuple) = $fnᵢ(; args...)
        end
    end
end

# Plastic potential
@inline QI(; Ptrial, σT, δ, kwargs...) = -Ptrial - @muladd(σT - δ * σT)
@inline QII(; Ptrial, PC1, τIItrial, τC1, β, ηve, dt, kwargs...) = √(((τIItrial - τC1) / ηve)^2 + ((-Ptrial + PC1) * β / dt)^2)
@inline QIII(; Ptrial, τIItrial, kwargs...) = τIItrial - Ptrial
@inline QIV(; Ptrial, PC2, τIItrial, τC2, β, ηve, dt, kwargs...) = √(((τIItrial - τC2) / ηve)^2 + ((-Ptrial + PC2) * β / dt)^2)
@inline QV(; Ptrial, τIItrial, sinΨ, kwargs...) = @muladd τIItrial - Ptrial * sinΨ

# Plastic potential branching checks
@inline check_QI(; τIItrial, τC1, kwargs...) = τIItrial ≤ τC1
@inline check_QII(; Ptrial, PC1, τIItrial, τC1, β, ηve, dt, kwargs...) = @muladd τC1 < τIItrial ≤ ηve * β / dt * (-Ptrial + PC1) + τC1
@inline check_QIII(; Ptrial, PC1, PC2, τIItrial, τC1, τC2, β, ηve, dt, kwargs...) = @muladd ηve * β / dt * (-Ptrial + PC1) + τC1 < τIItrial ≤ ηve * β / dt * (-Ptrial + PC2) + τC2
@inline check_QIV(; Ptrial, PC2, τIItrial, τC2, β, ηve, sinϕ, dt, kwargs...) = @muladd ηve * β / dt * (-Ptrial + PC2) + τC2 < τIItrial ≤ ηve * β / (dt * sinϕ) * (-Ptrial + PC2) + τC2 

# plastic strain rate
@inline function plastic_strain_rateI(; Ptrial, σT, δ, χ, kwargs...)
    ε_dev = 0.0
    ε_vol =  χ * QI(; Ptrial, σT, δ) 
    return ε_dev, ε_vol
end

@inline @inline function plastic_strain_rateII(; Ptrial, PC, τIItrial, τC, β, ηve, dt, χ, kwargs...)
    ε_dev = ((τIItrial - τC) / ηve) * χ
    ε_vol = ((-Ptrial + PC) * β / dt)  * χ
    return ε_dev, ε_vol
end

@inline function plastic_strain_rateIII(; Ptrial, τIItrial, σT, ηve, β, dt, χ, kwargs...)
    ε_vol = (τIItrial - Ptrial - σT) / @muladd (ηve + inv(β) * dt) * χ
    ε_dev = 0.5 * ε_vol
    return ε_dev, ε_vol
end

@inline function plastic_strain_rateIV(; Ptrial, PC2, τIItrial, τC2, ηve, β, dt, χ, kwargs...)
    ε_dev = (τIItrial - τC2) / ηve * χ
    ε_vol = (-Ptrial - PC2) * β / dt * χ
    return ε_dev, ε_vol
end

@inline function plastic_strain_rateV(; Ptrial, C, τIItrial, sinϕ, cosϕ, sinΨ, ηve, β, dt, χ, kwargs...)
    a = @muladd  ((τIItrial - Ptrial * sinϕ - C * cosϕ) / (ηve + inv(β) * dt * sinϕ * sinΨ)) * χ
    ε_dev = 0.5 * a
    ε_vol = sinϕ * a
    return ε_dev, ε_vol
end

@inline function yield_function(p::TensilePlasticity, P, τII)
    (; C, sinϕ, cosϕ, δ, σT) = p
    return max(
        @muladd(τII - P * sinϕ - C * cosϕ),
        τII - P - σT,
        - P - fma(-δ, σT, σT)
    )
end

function compute_Q(p, args::NamedTuple)
    Q = if yield_function(p, args.Ptrial, args.τIItrial) > 0
        Q = if check_QI(args)
            QI(args)

        elseif check_QII(args)
            QII(args)

        elseif check_QIII(args)
            QIII(args)

        elseif check_QIV(args)
            QIV(args)
        
        else
            QV(args)

        end
    else
        0.0
    end
end

function compute_plastic_strain(p, args::NamedTuple)
    ε_dev, ε_vol = if yield_function(p, args.Ptrial, args.τIItrial) > 0
        ε_dev, ε_vol = if check_QI(args)
            plastic_strain_rateI(args)

        elseif check_QII(args)
            plastic_strain_rateII(args)

        elseif check_QIII(args)
            plastic_strain_rateIII(args)

        elseif check_QIV(args)
            plastic_strain_rateIV(args)
        
        else
            plastic_strain_rateIV(args)

        end
    else
        0.0, 0.0
    end
end

function plastic_pressure_correction(P_tr, ε_vol_pl, β, dt)
    P = P_tr + dt / β * ε_vol_pl
    return P
end

function plastic_stressII_correction(τII_tr, εII_pl, ηve)
    τII = τII_tr - 2.0 * ηve * εII_pl
    return τII
end

function plastic_stress_corrections_normal(τ_tr, τII_tr, P_tr, P, εII_pl, ηve)
    τ = -P + (τ_tr + P_tr) * (1.0 - 2.0 * ηve * εII_pl /τII_tr)
    return τ
end

function plastic_stress_corrections_shear(τ_tr, τII_tr, εII_pl, ηve)
    τ =  τ_tr*(1.0 - 2.0 * ηve * εII_pl / τII_tr)
    return τ
end

function plastic_corrections(P_tr, τ_ii_tr::NTuple{N1, T}, τ_ij_tr::NTuple{N2, T}, ε_vol_pl, β, dt, τII_tr, εII_pl, ηve) where {N1, N2, T}
    P       = plastic_pressure_correction(P_tr, ε_vol_pl, β, dt)
    τII     = plastic_stressII_correction(τII_tr, εII_pl, ηve)
    τnormal = ntuple((Val(N1))) do i 
        plastic_stress_corrections_normal(τ_ii_tr[i], τII_tr, P_tr, P, εII_pl, ηve)
    end
    τshear  = ntuple((Val(N2))) do i 
        plastic_stress_corrections_normal(τ_ij_tr[i], τII_tr, P_tr, P, εII_pl, ηve)
    end
    return P, τII, τnormal, τshear
end


x = C, ϕ, Ψ, PC1, PC2, τC1, τC2, β, ηve, σT, δ = tuple(rand(11)...)
p = TensilePlasticity(x...)
Ptrial = 1.0
τIItrial = 1.0
dt = 1.0
(; C, ϕ, Ψ, PC1, PC2, τC1, τC2, β, ηve, σT, δ, sinϕ, cosϕ, sinΨ, cosΨ) = p
args = (; Ptrial = Ptrial, τIItrial = τIItrial, dt = dt, χ = 1/2, C, ϕ, Ψ, PC1, PC2, τC1, τC2, β, ηve, σT, δ, sinϕ, cosϕ, sinΨ, cosΨ)

@be compute_plastic_strain($p, $args)

