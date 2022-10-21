export DruckerPrager

# DruckerPrager  -------------------------------------------------------
"""
    DruckerPrager(ϕ=30, Ψ=0, C=10e6Pa)

Sets parameters for Drucker-Prager plasticity, where the yield stress ``\\sigma_{y}`` is computed by
```math  
    \\sigma_{y} = (P-P_f)\\tan(ϕ) + C
```
with ``\\phi`` being the friction angle (in degrees), ``C`` cohesion, ``P`` dynamic pressure and ``P_f`` the fluid pressure (both positive under compression).  

*Yielding* occurs when the second invariant of the deviatoric stress tensor, ``\\tau_{II}=(0.5\\tau_{ij}\\tau_{ij})^{0.5}`` touches the yield stress. 
This can be computed with the yield function ``F`` and the plastic flow potential ``Q``, which are respectively given by 
```math  
    F = \\tau_{II} - \\cos(ϕ)C - \\sin(ϕ)(P-P_f)
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
@with_kw_noshow struct DruckerPrager{T, U, U1} <: AbstractPlasticity{T}
    ϕ::GeoUnit{T,U} = 30NoUnits # Friction angle
    Ψ::GeoUnit{T,U} = 0NoUnits # Dilation angle
    sinϕ::GeoUnit{T,U} = sind(ϕ)NoUnits # Friction angle
    cosϕ::GeoUnit{T,U} = cosd(ϕ)NoUnits # Friction angle
    sinΨ::GeoUnit{T,U} = sind(Ψ)NoUnits # Dilation angle
    cosΨ::GeoUnit{T,U} = cosd(Ψ)NoUnits # Dilation angle
    C::GeoUnit{T,U1} = 10e6Pa # Cohesion
end
DruckerPrager(args...) = DruckerPrager(convert.(GeoUnit, args)...)

function isvolumetric(s::DruckerPrager)
    @unpack_val Ψ = s
    return Ψ == 0 ? false : true
end

function param_info(s::DruckerPrager) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"F = \\tau_{II} - \\cos(ϕ)C - \\sin(ϕ)(P-P_f); Q=\\tau_{II} - \\sin(Ψ)(P-P_f)",
    )
end

# Calculation routines
function (s::DruckerPrager{_T,U,U1})(;
    P::_T=zero(_T), τII::_T=zero(_T), Pf::_T=zero(_T), kwargs...
) where {_T,U,U1}
    @unpack_val sinϕ, cosϕ, ϕ, C = s
    F = τII - cosϕ * C - sinϕ * (P - Pf)   # with fluid pressure (set to zero by default)
    return F
end

"""
    compute_yieldfunction(s::DruckerPrager; P, τII, Pf, kwargs...) 

Computes the plastic yield function `F` for a given second invariant of the deviatoric stress tensor `τII`,  `P` pressure, and `Pf` fluid pressure.
"""
function compute_yieldfunction(
    s::DruckerPrager{_T}; P::_T=zero(_T), τII::_T=zero(_T), Pf::_T=zero(_T)
) where {_T}
    return s(; P=P, τII=τII, Pf=Pf)
end

"""
    compute_yieldfunction!(F::AbstractArray{_T,N}, s::DruckerPrager{_T}; P::AbstractArray{_T,N}, τII::AbstractArray{_T,N}, Pf=zero(P)::AbstractArray{_T,N}, kwargs...) 

Computes the plastic yield function `F` for Drucker-Prager plasticity in an in-place manner.
Required input arrays are pressure `P` and the second invariant of the deviatoric stress tensor `τII` at every point. 
You can optionally provide an array with fluid pressure `Pf` as well. 
"""
function compute_yieldfunction!(
    F::AbstractArray{_T,N},
    s::DruckerPrager{_T};
    P::AbstractArray{_T,N},
    τII::AbstractArray{_T,N},
    Pf=zero(P)::AbstractArray{_T,N},
    kwargs...,
) where {N,_T}
    @inbounds for i in eachindex(P)
        F[i] = compute_yieldfunction(s; P=P[i], τII=τII[i], Pf=Pf[i])
    end

    return nothing
end

# Plastic Potential 

# Derivatives w.r.t pressure
∂Q∂P(p::DruckerPrager, args) = -NumValue(p.sinΨ)

# Derivatives of yield function
∂F∂τII(p::DruckerPrager, τII::_T) where _T  = _T(1)
∂F∂P(p::DruckerPrager, P::_T) where _T      = -NumValue(p.sinϕ)
∂F∂λ(p::DruckerPrager, τII::_T) where _T    = _T(0)


# Derivatives w.r.t stress tensor

# Hard-coded partial derivatives of the plastic potential Q
for t in (:NTuple,:SVector)
    @eval begin
        ## 3D derivatives 
        ∂Q∂τxx(p::DruckerPrager, τij::$(t){6, T}) where T = 0.5 * τij[1] / second_invariant(τij)
        ∂Q∂τyy(p::DruckerPrager, τij::$(t){6, T}) where T = 0.5 * τij[2] / second_invariant(τij)
        ∂Q∂τzz(p::DruckerPrager, τij::$(t){6, T}) where T = 0.5 * τij[3] / second_invariant(τij)
        ∂Q∂τyz(p::DruckerPrager, τij::$(t){6, T}) where T = τij[4] / second_invariant(τij)
        ∂Q∂τxz(p::DruckerPrager, τij::$(t){6, T}) where T = τij[5] / second_invariant(τij)
        ∂Q∂τxy(p::DruckerPrager, τij::$(t){6, T}) where T = τij[6] / second_invariant(τij) 
        ## 2D derivatives 
        ∂Q∂τxx(p::DruckerPrager, τij::$(t){3, T}) where T = 0.5 * τij[1] / second_invariant(τij)
        ∂Q∂τyy(p::DruckerPrager, τij::$(t){3, T}) where T = 0.5 * τij[2] / second_invariant(τij)
        ∂Q∂τxy(p::DruckerPrager, τij::$(t){3, T}) where T = τij[3] / second_invariant(τij) 
    end
end

∂Q∂τII(p::DruckerPrager, τII::T) where T = 0.5

"""
    compute_εII(p::DruckerPrager{_T,U,U1}, λdot::_T, τII::_T,  P) 

This computes plastic strain rate invariant for a given ``λdot``
"""
function compute_εII(p::DruckerPrager{_T,U,U1}, λdot::_T, τII::_T, kwargs...) where {_T, U, U1}
    F = compute_yieldfunction(p, kwargs)
    @show F, kwargs
    if F>0
        ε_pl = λdot*∂Q∂τII(p, τII)

    else
        ε_pl = 0.0
    end 

    return ε_pl
end

# Print info 
function show(io::IO, g::DruckerPrager)
    return print(
        io,
        "Drucker-Prager plasticity with: C = $(UnitValue(g.C)), ϕ = $(UnitValue(g.ϕ))ᵒ, Ψ = $(UnitValue(g.Ψ))ᵒ",
    )
end
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------

# Plastic multiplier 
# NOTE: this is not used in the nonlinear iterations, so may not have to be defined for all cases

"""
    lambda(F::T, p::DruckerPrager, ηve::T, ηvp::T; K=zero(T), dt=zero(T), h=zero(T), τij=(one(T), one(T), one(T)))
    
    Compute the plastic multiplier λ for a Drucker-Prager yield surface. `F` is the trial yield surface, 
    `ηve` is the visco-elastic effective viscosity (i.e. `(1/G/dt + 1/η)⁻¹`), `ηvp` is a regularization term,
    `K` is the elastic bulk modulus, h is the harderning, and `τij`` is the stress tensor in Voigt notation.
    Equations from Duretz et al. 2019 G3
"""
@inline function lambda(F::T, p::DruckerPrager, ηve::T, ηvp::T; K=zero(T), dt=zero(T), h=zero(T), τij=(one(T), one(T), one(T))) where T
    F * inv(ηve + ηvp + K * dt * p.sinΨ * p.sinϕ + h * p.cosϕ * plastic_strain(p, τij, zero(T)))
end
#-------------------------------------------------------------------------