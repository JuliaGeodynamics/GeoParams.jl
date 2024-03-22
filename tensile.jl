using ForwardDiff, MuladdMacro, GeoParams, LinearAlgebra

# struct TensilePlasticity{T}
#     C::T
#     ϕ::T
#     sinϕ::T
#     cosϕ::T
#     Ψ::T
#     sinΨ::T
#     cosΨ::T
#     x1σ ::T
#     x2σ ::T
#     x1τ ::T
#     x2τ ::T
#     β ::T
#     ηve::T
#     σ_T::T
#     δ::T
#     function TensilePlasticity(
#         C::T, ϕ::T, Ψ::T, x1σ::T, x2σ::T, x1τ::T, x2τ::T, β::T, ηve::T, σ_T::T, δ::T
#     ) where T
#         sinϕ, cosϕ = sincos(ϕ)
#         sinΨ, cosΨ = sincos(Ψ)
#         new{T}(C, ϕ, sinϕ, cosϕ, Ψ, sinΨ, cosΨ, x1σ, x2σ, x1τ, x2τ, β, ηve, σ_T, δ)
#     end
# end

struct TensilePlasticity{T}
    C::T
    ϕ::T
    sinϕ::T
    cosϕ::T
    Ψ::T
    sinΨ::T
    cosΨ::T
    β ::T
    δσ_T::T
    σ_T_to_C::T

    function TensilePlasticity(
        C::T, ϕ::T, Ψ::T, β::T, δσ_T::T, σ_T_to_C::T
    ) where T<:AbstractFloat
        sinϕ, cosϕ = sincos(ϕ)
        sinΨ, cosΨ = sincos(Ψ)
        new{T}(C, ϕ, sinϕ, cosϕ, Ψ, sinΨ, cosΨ, β, δσ_T, σ_T_to_C)
    end
end

@inline TensilePlasticity(args...) = TensilePlasticity(promote(args...)...)

function tensile_parameters(p::TensilePlasticity)
    (; σ_T_to_C, C, δσ_T, sinϕ, cosϕ) = p

    σ_T  =  σ_T_to_C * C
    x1σ  =  σ_T - δσ_T                 # σm at the intersection of cutoff and Mode-1
    x1τ  = -x1σ + σ_T                  # τII at the intersection of cutoff and Mode-1
    x2σ  =  (σ_T - C*cosϕ)/(1.0-sinϕ)  # σm at the intersection of Drucker-Prager and Mode-1
    x2τ  = -x2σ + σ_T                  # τII at the intersection of Drucker-Prager and Mode-1

    return σ_T, x1σ, x1τ, x2σ, x2τ
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
@inline QI(; Ptrial, σ_T, δ, kwargs...) = -Ptrial - @muladd(σ_T - δσ_T)
@inline QII(; Ptrial, x1σ, τII_tr, x1τ, β, ηve, dt, kwargs...) = √(((τII_tr - x1τ) / ηve)^2 + ((-Ptrial + x1σ) * β / dt)^2)
@inline QIII(; Ptrial, τII_tr, kwargs...) = τII_tr - Ptrial
@inline QIV(; Ptrial, x2σ, τII_tr, x2τ, β, ηve, dt, kwargs...) = √(((τII_tr - x2τ) / ηve)^2 + ((-Ptrial + x2σ) * β / dt)^2)
@inline QV(; Ptrial, τII_tr, sinΨ, kwargs...) = @muladd τII_tr - Ptrial * sinΨ

# Plastic potential branching checks
@inline check_QI(; τII_tr, x1τ, kwargs...) = τII_tr ≤ x1τ
@inline check_QII(; Ptrial, x1σ, τII_tr, x1τ, β, ηve, dt, kwargs...) = @muladd x1τ < τII_tr ≤ ηve * β / dt * (-Ptrial + x1σ) + x1τ
@inline check_QIII(; Ptrial, x1σ, x2σ, τII_tr, x1τ, x2τ, β, ηve, dt, kwargs...) = @muladd ηve * β / dt * (-Ptrial + x1σ) + x1τ < τII_tr ≤ ηve * β / dt * (-Ptrial + x2σ) + x2τ
@inline check_QIV(; Ptrial, x2σ, τII_tr, x2τ, β, ηve, sinϕ, dt, kwargs...) = @muladd ηve * β / dt * (-Ptrial + x2σ) + x2τ < τII_tr ≤ ηve * β / (dt * sinϕ) * (-Ptrial + x2σ) + x2τ 

# plastic strain rate
@inline function plastic_strain_rateI(; Ptrial, σ_T, δσ_T, χ, kwargs...)
    ε_dev = 0.0
    ε_vol =  χ * QI(; Ptrial, σ_T, δσ_T) 
    return ε_dev, ε_vol
end

@inline @inline function plastic_strain_rateII(; Ptrial, x1σ, τII_tr, x1τ, β, ηve, dt, χ, kwargs...)
    ε_dev = ((τII_tr - x1τ) / ηve) * χ
    ε_vol = ((-Ptrial + x1σ) * β / dt)  * χ
    return ε_dev, ε_vol
end

@inline function plastic_strain_rateIII(; Ptrial, τII_tr, σ_T, ηve, β, dt, χ, kwargs...)
    ε_vol = (τII_tr - Ptrial - σ_T) / @muladd (ηve + inv(β) * dt) * χ
    ε_dev = 0.5 * ε_vol
    return ε_dev, ε_vol
end

@inline function plastic_strain_rateIV(; Ptrial, x2σ, τII_tr, x2τ, ηve, β, dt, χ, kwargs...)
    ε_dev = (τII_tr - x2τ) / ηve * χ
    ε_vol = (-Ptrial - x2σ) * β / dt * χ
    return ε_dev, ε_vol
end

@inline function plastic_strain_rateV(; Ptrial, C, τII_tr, sinϕ, cosϕ, sinΨ, ηve, β, dt, χ, kwargs...)
    a = @muladd  ((τII_tr - Ptrial * sinϕ - C * cosϕ) / (ηve + inv(β) * dt * sinϕ * sinΨ)) * χ
    ε_dev = 0.5 * a
    ε_vol = sinϕ * a
    return ε_dev, ε_vol
end

@inline function yield_function(p::TensilePlasticity, P, τII, σ_T, δσ_T)
    (; C, sinϕ, cosϕ, δσ_T) = p
    return max(
        @muladd(τII - P * sinϕ - C * cosϕ),
        τII - P - σ_T,
        @muladd(-P -(σ_T - δσ_T))
    )
end

function compute_Q(p, args::NamedTuple)
    Q = if yield_function(p, args.Ptrial, args.τII_tr) > 0
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
    ε_dev, ε_vol = if yield_function(p, args.Ptrial, args.τII_tr, args.σ_T, args.δσ_T) > 0
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

function plastic_pressure_correction(Ptrial, ε_vol, β, dt)
    P = Ptrial + dt / β * ε_vol
    return P
end

function plastic_stressII_correction(τII_tr, ε_dev, ηve)
    τII = τII_tr - 2.0 * ηve * ε_dev
    return τII
end

function plastic_stress_corrections_normal(τ_tr, τII_tr, Ptrial, P, ε_dev, ηve)
    τ = -P + (τ_tr + Ptrial) * (1.0 - 2.0 * ηve * ε_dev /τII_tr)
    return τ
end

function plastic_stress_corrections_shear(τ_tr, τII_tr, ε_dev, ηve)
    τ =  τ_tr*(1.0 - 2.0 * ηve * ε_dev / τII_tr)
    return τ
end

function plastic_corrections(Ptrial, τ_normal_tr::NTuple{N1, T}, τ_shear_tr::NTuple{N2, T}, τII_tr, ε_dev, ε_vol, β, ηve, dt) where {N1, N2, T}
    P       = plastic_pressure_correction(Ptrial, ε_vol, β, dt)
    τII     = plastic_stressII_correction(τII_tr, ε_dev, ηve)
    τnormal = ntuple((Val(N1))) do i 
        plastic_stress_corrections_normal(τ_normal_tr[i], τII_tr, Ptrial, P, ε_dev, ηve)
    end
    τshear  = ntuple((Val(N2))) do i 
        plastic_stress_corrections_normal(τ_shear_tr[i], τII_tr, Ptrial, P, ε_dev, ηve)
    end
    return P, τII, (τnormal..., τshear...)
end


# # function foo()
    # C, ϕ, Ψ, β, δσ_T, σ_T_to_C = rand(6)
    yr       = 365.25 * 24 * 3600
    MPa      = 1e6
    dt       = 1e3 * yr
    C        = 15e6
    ϕ        = 30.0
    Ψ        = 15.0
    β        = 3*2.8e-11
    G        = 3.0/β
    η        = 1e30
    ηve      = 1/(1/η + 1/G/dt)
    δσ_T     = 1.5
    σ_T_to_C = 0.5

    p = TensilePlasticity(C, ϕ, Ψ, β, δσ_T, σ_T_to_C)
    σ_T, x1σ, x1τ, x2σ, x2τ = tensile_parameters(p)

    τ_xx, τ_yy, τ_xy = (500, -500, 0) .* MPa
    P = 1e3 * MPa
    τII = second_invariant(τ_xx, τ_yy, τ_xy)
    (; C, ϕ, Ψ, β, δσ_T, sinϕ, cosϕ, sinΨ, cosΨ) = p

    for it in 1:100
        Ptrial = P
        τII_tr = τII
        args = (; Ptrial = Ptrial, τII_tr = τII_tr, dt = dt, χ = 1/2, C, ϕ, Ψ, x1σ, x2σ, x1τ, x2τ, β, ηve, δσ_T, σ_T, sinϕ, cosϕ, sinΨ, cosΨ)

        ε_dev, ε_vol = compute_plastic_strain(p, args)
        P, τII, τij = plastic_corrections(Ptrial, (τ_xx, τ_yy), (τ_xy,), τII, ε_dev, ε_vol, p.β, args.ηve, dt)
        τ_xx, τ_yy, τ_xy = τij
        rP   = norm(P-Ptrial, 2)
        rτII = norm(τII-τII_tr, 2)

        if rP < 1e-6 && rτII < 1e-6
            println("Converged at iteration $it")
            println("    rP = $rP, rτII = $rτII")
            break
        end
    end
# # end


# C, ϕ, Ψ, x1σ, x2σ, x1τ, x2τ, β, ηve, σ_T, δ = rand(11)
# p = TensilePlasticity(C, ϕ, Ψ, x1σ, x2σ, x1τ, x2τ, β, ηve, σ_T, δ)
# Ptrial = 1.0
# τII_tr = 1.0
# dt = 1.0
# (; C, ϕ, Ψ, x1σ, x2σ, x1τ, x2τ, β, ηve, σ_T, δ, sinϕ, cosϕ, sinΨ, cosΨ) = p
# args = (; Ptrial = Ptrial, τII_tr = τII_tr, dt = dt, χ = 1/2, C, ϕ, Ψ, x1σ, x2σ, x1τ, x2τ, β, ηve, σ_T, δ, sinϕ, cosϕ, sinΨ, cosΨ)

# compute_plastic_strain(p, args)

# # @b foo()


# # t_vec    = range(t[1], t[2], length=nt)
# # τII_vec  = zero(t_vec)
# # ε_vec  = [ε for i = 1:nt]
# # τ_vec  = [τ0 for i = 1:nt]

# # time_τII_0D!(τ_vec, τII_vec, v, ε_vec, args, t_vec, verbose=verbose);

# # function time_τII_0D!(τ_vec::Vector{T}, v::Union{CompositeRheology,Tuple, Parallel}, εII_vec::Vector{T}, args, t_vec::AbstractVector{T}; verbose=false) where {T}

# #     nt  = length(τ_vec)
# #     τII = τ_vec[1]

# #     for i=2:nt  
# #         dt      = t_vec[i]-t_vec[i-1]
# #         args    = merge(args, (; τII_old=τII, dt=dt))
# #         τII,    = compute_τII(v, εII_vec[i-1], args, verbose=verbose)
        
# #         τ_vec[i] = τII
# #     end

# #     return nothing
# # end