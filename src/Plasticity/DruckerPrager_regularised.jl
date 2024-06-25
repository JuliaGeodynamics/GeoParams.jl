export DruckerPrager_regularised

# DruckerPrager_regularised  -------------------------------------------------------
"""
    DruckerPrager_regularised(ϕ=30, Ψ=0, C=10e6Pa, η_vp=1e20Pa*s)

Sets parameters for reularised Drucker-Prager plasticity, where the yield stress ``\\sigma_{y}`` is computed by
```math  
    \\sigma_{y} = (P-P_f)\\tan(ϕ) + C + 2η_vpε̇II_pl 
```
with ``\\phi`` being the friction angle (in degrees), ``C`` cohesion, ``P`` dynamic pressure, ``P_f`` the fluid pressure (both positive under compression), ``η_vp`` the regularization viscosity and ``ε̇II_pl`` the invariant of the plastic strainrate

*Yielding* occurs when the second invariant of the deviatoric stress tensor, ``\\tau_{II}=(0.5\\tau_{ij}\\tau_{ij})^{0.5}`` touches the yield stress. 
This can be computed with the yield function ``F`` and the plastic flow potential ``Q``, which are respectively given by 
```math  
    F = \\tau_{II} - \\cos(ϕ)C - \\sin(ϕ)(P-P_f) - 2 \\eta_{vp} \\dot{\\varepsilon}ε̇^{pl}_{II}
```
```math  
    Q = \\tau_{II} - \\sin(Ψ)(P-P_f) 
```
Here, Ψ is the dilation angle, which must be zero for incompressible setups.

Plasticity is activated when ``F(\\tau_{II}^{trial})`` (the yield function computed with a trial stress) is >0. In that case, plastic strainrate ``\\dot{\\varepsilon}^{pl}_{ij}`` is computed by:
```math  
    \\dot{\\varepsilon}^{pl}_{ij} =\\dot{\\lambda} {\\partial Q \\over \\partial \\sigma_{ij}}
```
where ``\\dot{\\lambda}`` is a (scalar) that is nonzero and chosen such that the resulting stress gives ``F(\\tau_{II}^{final})=0``, and ``\\sigma_{ij}=-P + \\tau_{ij}`` denotes the total stress tensor.   
        
"""
@with_kw_noshow struct DruckerPrager_regularised{T, U, U1, U2, S1<:AbstractSoftening, S2<:AbstractSoftening} <: AbstractPlasticity{T}
    softening_ϕ::S1 = NoSoftening()
    softening_C::S2 = NoSoftening()
    ϕ::GeoUnit{T,U} = 30NoUnits         # Friction angle
    Ψ::GeoUnit{T,U} = 0NoUnits          # Dilation angle
    sinϕ::GeoUnit{T,U} = sind(ϕ)NoUnits # Friction angle
    cosϕ::GeoUnit{T,U} = cosd(ϕ)NoUnits # Friction angle
    sinΨ::GeoUnit{T,U} = sind(Ψ)NoUnits # Dilation angle
    cosΨ::GeoUnit{T,U} = cosd(Ψ)NoUnits # Dilation angle
    C::GeoUnit{T,U1} = 10e6Pa           # Cohesion
    η_vp::GeoUnit{T,U2} = 1e20Pa*s      # regularisation viscosity
end

DruckerPrager_regularised(args...) = DruckerPrager_regularised(args[1:2]..., convert.(GeoUnit, args[3:end])...)
DruckerPrager_regularised(softening_ϕ::AbstractSoftening, softening_C::AbstractSoftening, args::Vararg{GeoUnit, N}) where N = DruckerPrager_regularised(softening_ϕ, softening_C, convert.(GeoUnit, args)...)

function isvolumetric(s::DruckerPrager_regularised)
    @unpack_val Ψ = s
    return Ψ == 0 ? false : true
end

number_plastic_variables(s::DruckerPrager_regularised) = Val(3)


function param_info(s::DruckerPrager_regularised) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"F = \\tau_{II} - \\cos(ϕ)C - \\sin(ϕ)(P-P_f) - 2η_vpε̇II_pl ; Q=\\tau_{II} - \\sin(Ψ)(P-P_f)",
    )
end

# Calculation routines
function (s::DruckerPrager_regularised{_T})(;
    P=zero(_T), τII=zero(_T), Pf=zero(_T), λ= zero(_T), EII=zero(_T), kwargs...
) where {_T}
    @unpack_val sinϕ, cosϕ, ϕ, C, η_vp = s
    ϕ = s.softening_ϕ(EII, ϕ)
    C = s.softening_C(EII, C)
    cosϕ, sinϕ = iszero(EII) ? (cosϕ, sinϕ) : (cosd(ϕ), sind(ϕ))
    ε̇II_pl = λ*∂Q∂τII(s, τII)  # plastic strainrate
    F = τII - cosϕ * C - sinϕ * (P - Pf)  - 2 * η_vp * ε̇II_pl # with fluid pressure (set to zero by default)
    return F
end

function (s::DruckerPrager_regularised{_T,U,U1,NoSoftening, S})(;
    P=zero(_T), τII=zero(_T), Pf=zero(_T), λ= zero(_T), EII=zero(_T), kwargs...
) where {_T,U,U1,S}
    @unpack_val sinϕ, cosϕ, ϕ, C, η_vp = s
    C = s.softening_C(EII, C)
    ε̇II_pl = λ*∂Q∂τII(s, τII)  # plastic strainrate
    F = τII - cosϕ * C - sinϕ * (P - Pf)  - 2 * η_vp * ε̇II_pl # with fluid pressure (set to zero by default)
    return F
end

function (s::DruckerPrager_regularised{_T,U,U1,S, NoSoftening})(;
    P=zero(_T), τII=zero(_T), Pf=zero(_T), λ= zero(_T), EII=zero(_T), kwargs...
) where {_T,U,U1,S}
    @unpack_val sinϕ, cosϕ, ϕ, C, η_vp = s
    ϕ = s.softening_ϕ(EII, ϕ)
    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    ε̇II_pl = λ*∂Q∂τII(s, τII)  # plastic strainrate
    F = τII - cosϕ * C - sinϕ * (P - Pf)  - 2 * η_vp * ε̇II_pl # with fluid pressure (set to zero by default)
    return F
end

function (s::DruckerPrager_regularised{_T,U,U1,NoSoftening,NoSoftening})(;
    P=zero(_T), τII=zero(_T), Pf=zero(_T), λ= zero(_T), kwargs...
) where {_T,U,U1}
    @unpack_val sinϕ, cosϕ, ϕ, C, η_vp = s
    ε̇II_pl = λ*∂Q∂τII(s, τII)  # plastic strainrate
    F = τII - cosϕ * C - sinϕ * (P - Pf)  - 2 * η_vp * ε̇II_pl # with fluid pressure (set to zero by default)
    # F = fma(-cosϕ, C, τII) - fma(sinϕ, (P - Pf),  -2.0 * η_vp * ε̇II_pl) # with fluid pressure (set to zero by default)
    return F
end


"""
    compute_yieldfunction(s::DruckerPrager_regularised; P, τII, Pf, λ, kwargs...) 

Computes the plastic yield function `F` for a given second invariant of the deviatoric stress tensor `τII`,  `P` pressure, and `Pf` fluid pressure.
"""
function compute_yieldfunction(
    s::DruckerPrager_regularised{_T}; P=zero(_T), τII=zero(_T), Pf=zero(_T), λ=zero(_T), EII=zero(_T)
) where {_T}
    return s(; P=P, τII=τII, Pf=Pf, λ=λ, EII=EII)
end

"""
    compute_yieldfunction!(F::AbstractArray{_T,N}, s::DruckerPrager_regularised{_T}; P::AbstractArray{_T,N}, τII::AbstractArray{_T,N}, Pf=zero(P)::AbstractArray{_T,N}, λ=zero(P)::AbstractArray{_T,N}, kwargs...) 

Computes the plastic yield function `F` for Drucker-Prager plasticity in an in-place manner.
Required input arrays are pressure `P` and the second invariant of the deviatoric stress tensor `τII` at every point. 
You can optionally provide an array with fluid pressure `Pf` as well. 
"""
function compute_yieldfunction!(
    F::AbstractArray{_T,N},
    s::DruckerPrager_regularised{_T};
    P::AbstractArray{_T,N},
    τII::AbstractArray{_T,N},
    Pf=zero(P)::AbstractArray{_T,N},
    λ=zero(P)::AbstractArray{_T,N},
    EII::AbstractArray{_T,N} = zero(P),
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(P)
        F[i] = compute_yieldfunction(s; P=P[i], τII=τII[i], Pf=Pf[i], λ=λ[i], EII=EII[i])
    end

    return nothing
end

# Plastic Potential 

# Derivatives w.r.t pressure
∂Q∂P(p::DruckerPrager_regularised, args; kwargs...) = NumValue(p.sinΨ)

# Derivatives of yield function
∂F∂τII(p::DruckerPrager_regularised, τII::_T; kwargs...) where _T             = _T(1)
∂F∂P(p::DruckerPrager_regularised,     P::_T; kwargs...) where _T             = -NumValue(p.sinϕ)
∂F∂λ(p::DruckerPrager_regularised,   τII::_T; P=zero(_T), kwargs...) where _T = -2 * NumValue(p.η_vp)*∂Q∂τII(p, τII, P=P) 


"""
This adds the plasticity contributions to the local residual vector. 
"""
function add_plastic_residual!(r::_T, x::_T1, c::DruckerPrager_regularised, kwargs...) where {_T, _T1}   
    τ   = x[1]  # second invariant of stress; always the 1th entry 
    P   = x[2]  # second invariant of stress; always the 1th entry 
    λ   = x[3]
    @show λ
    
    args = merge(args, (λ=λ,τII=τ, P=P))
    F = compute_yieldfunction(v, args)
    
    if F>0
        εII_pl  = λ*∂Q∂τII(c, τ, args)
        εvol_pl = λ*∂Q∂P(c, P, args)
    else
        εII_pl = 0.0
        εvol_pl = 0.0
    end 
    @show εvol_pl εII_pl ∂Q∂P(c[3], P, args) ∂Q∂τII(c[3], τ, args)
    r[1]  -= εII_pl
    r[2]  -= εvol_pl 
    r[3]   = F

    return nothing
end

"""
    add_plastic_jacobian!(J::_T, x::_T1, c::DruckerPrager_regularised, args)

This adds the plasticity contributions to the local jacobian matrix. Note that defining this function is optional; one can also compute it using AD.
"""
function add_plastic_jacobian!(J::_T, x::_T1, v::DruckerPrager_regularised, args) where {_T, _T1}   
    @unpack_val η_vp, kf, pf, b = v
    τ   = x[1]  # second invariant of stress; always the 1th entry 
    P   = x[2]  # pressure; always the 2nd entry 
    λ   = x[3]  # plastic multiplier
    
    dQdτ = ∂Q∂τII(v, τ, args)
    dQdP = ∂Q∂P(v, P, args)
    dFdτ = ∂F∂τII(v, τ, args)
    dFdP = ∂F∂P(v, P, args)

    F = compute_yieldfunction(v, args)

    J[1, 1] += 0.0 
    J[1, 2] += 0.0
    J[1, 3] += -dQdτ
    
    J[2, 1] += 0.0
    J[2, 2] += 0.0
    J[2, 3] += -dQdP
    
    J[3, 1] +=  dFdτ
    J[3, 2] +=  dFdP
    J[3, 3] +=  -η_vp
    return nothing
end



# Derivatives w.r.t stress tensor

# Hard-coded partial derivatives of the plastic potential Q
for t in (:NTuple,:SVector)
    @eval begin
        ## 3D derivatives 
        ∂Q∂τxx(::DruckerPrager_regularised, τij::$(t){6, T}) where T = 0.5 * τij[1] / second_invariant(τij)
        ∂Q∂τyy(::DruckerPrager_regularised, τij::$(t){6, T}) where T = 0.5 * τij[2] / second_invariant(τij)
        ∂Q∂τzz(::DruckerPrager_regularised, τij::$(t){6, T}) where T = 0.5 * τij[3] / second_invariant(τij)
        ∂Q∂τyz(::DruckerPrager_regularised, τij::$(t){6, T}) where T = τij[4] / second_invariant(τij)
        ∂Q∂τxz(::DruckerPrager_regularised, τij::$(t){6, T}) where T = τij[5] / second_invariant(τij)
        ∂Q∂τxy(::DruckerPrager_regularised, τij::$(t){6, T}) where T = τij[6] / second_invariant(τij) 
        ## 2D derivatives 
        ∂Q∂τxx(::DruckerPrager_regularised, τij::$(t){3, T}) where T = 0.5 * τij[1] / second_invariant(τij)
        ∂Q∂τyy(::DruckerPrager_regularised, τij::$(t){3, T}) where T = 0.5 * τij[2] / second_invariant(τij)
        ∂Q∂τxy(::DruckerPrager_regularised, τij::$(t){3, T}) where T = τij[3] / second_invariant(τij) 
    end
end

∂Q∂τII(p::DruckerPrager_regularised, τII::_T; P=zero(_T), kwargs...) where _T = 0.5

"""
    compute_εII(p::DruckerPrager_regularised{_T,U,U1}, λdot::_T, τII::_T,  P) 

This computes plastic strain rate invariant for a given ``λdot``
"""
function compute_εII(p::DruckerPrager_regularised{_T,U,U1}, λdot::_T, τII::_T, kwargs...) where {_T, U, U1}
    args = merge(kwargs, (λ=λdot,))
    F = compute_yieldfunction(p, args)
    if F>0
        ε_pl = λdot*∂Q∂τII(p, τII)

    else
        ε_pl = 0.0
    end 

    return ε_pl
end

# Print info 
function show(io::IO, g::DruckerPrager_regularised)
    return print(
        io,
        "Regularized Drucker-Prager plasticity with: C = $(UnitValue(g.C)), ϕ = $(UnitValue(g.ϕ))ᵒ, Ψ = $(UnitValue(g.Ψ))ᵒ, η_vp=$(UnitValue(g.η_vp))",
    )
end
#-------------------------------------------------------------------------

