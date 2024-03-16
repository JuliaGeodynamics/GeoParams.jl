using GeoParams, ForwardDiff, MuladdMacro
# # Plasticity
# @. F    = τii - Coh*cos(ϕ) - Pt*sin(ϕ)
# @. Ptc  = Pt
# @. ηvep = ηve
# @. ispl = 0
# @. ispl[F>=0] = 1
# @. ε̇iiᵉᶠᶠ   = sqrt( (ε̇xy + τxy0/2/ηe)^2 + 0.5*( (ε̇xxd + τxx0/(2*ηe))^2 + ((ε̇yyd + τyy0/(2*ηe))).^2 + ((ε̇zzd + τzz0/(2*ηe))).^2 ) ) 
# for it=1:50  
#     @. ηvep   = (Coh*cos(ϕ) + Ptc*sin(ϕ) + ηvp*λ̇rel) / 2.0 / ε̇iiᵉᶠᶠ # out
#     @. ε̇xxd_pl = λ̇rel*(τxx/2/τii)
#     @. ε̇yyd_pl = λ̇rel*(τyy/2/τii)
#     @. ε̇zzd_pl = λ̇rel*(τzz/2/τii) 
#     @. ε̇xy_pl  = λ̇rel*(τxy/2/τii) 
#     @. ∇v_pl   = sin(ψ)*λ̇rel
#     if params.oop == :Vermeer1990
#         @. ε̇zz_pl  = λ̇rel*(τxx/2/τii + τyy/2/τii)/2   # dqdτzz*λ̇rel
#         @. ∇v_pl   = 3/2*sin(ψ)*λ̇rel
#     end
#     @. Ptc    = Pt0  - Kb*Δt*(∇v - ∇v_pl) # out
#     @. τxx    =  2*ηve * (ε̇xxd + τxx0/(2*ηe) -  ε̇xxd_pl)         # out
#     @. τyy    =  2*ηve * (ε̇yyd + τyy0/(2*ηe) -  ε̇yyd_pl)         # out
#     @. τzz    =  2*ηve * (ε̇zzd + τzz0/(2*ηe) -  ε̇zzd_pl)         # out
#     @. τxy    =  2*ηve * (ε̇xy  + τxy0/(2*ηe) -  ε̇xy_pl)          # out
#     @. τii    = sqrt(τxy^2 + 0.5*(τyy^2 + τxx^2 + τzz^2))        # out
#     @. Fc     = τii - Coh*cos(ϕ) - Ptc*sin(ϕ) - ηvp*λ̇rel         # out
#     @. λ̇rel  += (F.>0) .* Fc / (ηvp + ηve + Kb*Δt*sin(ϕ)*sin(ψ)) # out
#     if maximum(Fc) < ϵ break end 
# end

@inline plastic_strain_rate(p, τij, λ) = ∂Q∂τ(p, τij) .* λ
@inline plastic_volume(p::DruckerPrager_regularised, λ) = p.sinΨ.val * λ

@inline function plastic_pressure_correction(p::DruckerPrager_regularised; P0 = 0e0, λ = 0e0, ∇v = 0e0, dt = 0e0, Kb = 0e0)
    P = P0 - Kb*dt*(∇v - plastic_volume(p, λ))
    return P
end

@inline function plastic_stress_corrections(τ0, ε, ε_pl, ηve, ηe)
    τ = 2 * ηve * (ε + τ0/(2*ηe) - ε_pl)
    return τ
end

@inline function compute_lambda(p::DruckerPrager_regularised, F, ηve,  Kb, dt)
    (; sinϕ, sinΨ, η_vp) = p
    λ = F / (η_vp + ηve + Kb * dt * sinϕ * sinΨ)
    return λ
end

p    = DruckerPrager_regularised()
P    = 0e0
τij  = τxx, τyy, τxy = eps(), eps(), eps()
τ0ij = τxx0, τyy0, τxy0 = 0e0, 0e0, 0e0 
εij  = εxx, εyy, εxy = eps(), eps(), eps()
τII  = second_invariant(τxx, τyy, τxy)
εII  = second_invariant(εxx, εyy, εxy)
λ    = 0e0
ηve, ηe = 1, 1
Kb, dt = 0e0, 1e0

# (1) trial yield function where λ = 0
F = compute_yieldfunction(p, (; P, τII, λ = λ)) 
# (2) plastic strain rate components
εij_pl = εxx_pl, εyy_pl, εxy_pl = plastic_strain_rate(p, τij, λ) 
# (3) correct pressure
Pnew = plastic_pressure_correction(p; P0 = P, λ = λ, ∇v = 0e0, dt = 1e0, Kb = 0e0)
# (4) correct stress
τij  = τxx, τyy, τxy = ntuple(Val(3)) do i 
    plastic_stress_corrections(τ0ij[i], εij[i], εij_pl[i], ηve, ηe)
end
# (5) compute stress second invariant
τII  = second_invariant(τxx, τyy, τxy)
# (6) update Fc
F = compute_yieldfunction(p, (; P, τII, λ = λ))
# (6) update λ
compute_lambda(p, F, ηve,  Kb, dt)


# function plastic_corrections(P_tr, τ_ii_tr::NTuple{N1, T}, τ_ij_tr::NTuple{N2, T}, ε_vol_pl, β, dt, τII_tr, εII_pl, ηve) where {N1, N2, T}
#     P       = plastic_pressure_correction(P_tr, ε_vol_pl, β, dt)
#     τII     = plastic_stressII_correction(τII_tr, εII_pl, ηve)
#     τnormal = ntuple((Val(N1))) do i 
#         plastic_stress_corrections_normal(τ_ii_tr[i], τII_tr, P_tr, P, εII_pl, ηve)
#     end
#     τshear  = ntuple((Val(N2))) do i 
#         plastic_stress_corrections_normal(τ_ij_tr[i], τII_tr, P_tr, P, εII_pl, ηve)
#     end
#     return P, τII, τnormal, τshear
# end