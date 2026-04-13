export DruckerPragerCap
export compute_tensile_cap, compute_flowpotential


# DruckerPragerCap  -------------------------------------------------------
"""
    DruckerPragerCap(ϕ=30, Ψ=0, C=10e6Pa, η_vp=1e20Pa*s, Pt=-1e5Pa)

Sets parameters for Drucker-Prager-Cap plasticity for mode-1 and mode-2 plasticity,
as described in Popov et al. (2025), Geoscientific Model Development.

# Fields
- `C::T`: The cohesion parameter.
- `ϕ::T`: The friction angle (in degrees).
- `Ψ::T`: The dilatancy angle (in degrees).
- `η_vp::T`: The Duvaut-Lions regularisation viscosity for the plasticity model.
- `pT::T`: The tensile strength (should be < 0).
"""
@with_kw_noshow struct DruckerPragerCap{T, U, U1, U2, S1 <: AbstractSoftening, S2 <: AbstractSoftening, S3 <: AbstractSoftening} <: AbstractPlasticity{T}
    softening_ϕ::S1 = NoSoftening()
    softening_C::S2 = NoSoftening()
    softening_Ψ::S3 = NoSoftening()
    ϕ::GeoUnit{T, U} = 30NoUnits         # Friction angle
    Ψ::GeoUnit{T, U} = 0NoUnits          # Dilation angle
    # computational parameters (precomputed, to speed up later calculations)
    sinϕ::GeoUnit{T, U} = sind(ϕ)NoUnits # Friction angle
    cosϕ::GeoUnit{T, U} = cosd(ϕ)NoUnits # Friction angle
    sinΨ::GeoUnit{T, U} = sind(Ψ)NoUnits # Dilation angle
    cosΨ::GeoUnit{T, U} = cosd(Ψ)NoUnits # Dilation angle
    C::GeoUnit{T, U1} = 10.0e6Pa           # Cohesion
    η_vp::GeoUnit{T, U2} = 1.0e20Pa * s      # regularisation viscosity
    pT::GeoUnit{T, U1} = -1.0e5Pa           # Tensile strength
end


DruckerPragerCap(args...) = DruckerPragerCap(args[1:3]..., convert.(GeoUnit, promote(args[4:end]...))...)
DruckerPragerCap(softening_ϕ::AbstractSoftening, softening_C::AbstractSoftening, softening_Ψ::AbstractSoftening, args...) = DruckerPragerCap(softening_ϕ, softening_C, softening_Ψ, convert.(GeoUnit, promote(args...))...)

function isvolumetric(s::DruckerPragerCap)
    @unpack_val Ψ = s
    return !iszero(Ψ)
end

function param_info(s::DruckerPragerCap)
    return MaterialParamsInfo(;
        Equation = L"F = \tau_{II} - kP - c \;\;\mathrm{or}\;\; a(\sqrt{\tau_{II}^2 + (P-p_y)^2} - R_y),\; Q = \tau_{II} - k_q P - \mathrm{const} \;\mathrm{or}\;\; b(\sqrt{\tau_{II}^2 + (P-p_q)^2} - R_f)",
    )
end


"""
    compute_tensile_cap(sinϕ, cosϕ, sinψ, C, pT)

Compute smooth tensile-cap parameters from Popov et al. (2025).
"""
@inline function compute_tensile_cap(sinϕ::T, cosϕ::T, sinψ::T, C::T, pT::T) where {T}
    k = sinϕ
    kq = sinψ
    c = C * cosϕ
    a = sqrt(one(T) + k^2)
    b = sqrt(one(T) + kq^2)
    cosa = inv(a)
    sina = k * cosa
    py = (pT + c * cosa) * inv(one(T) - sina)
    Ry = py - pT
    pd = @muladd py - Ry * sina
    τd = @muladd k * pd + c # delimiter point
    pq = @muladd pd + kq * τd
    ## plotting parameters (not needed for the model, but useful for plotting the yield surface and flow potential)
    Rf = pq - pT
    normvRf = sqrt((pd - pq)^2 + τd^2) * inv(Rf)
    pdf = @muladd (pd - pq) * inv(normvRf) + pq
    sdf = τd * inv(normvRf)

    return (; k, kq, c, a, b, py, Ry, pd, τd, pq, Rf, pdf, sdf, normvRf)
end

@inline ismode2_yield(py, pd, τd, τII, P) = τII * (py - pd) ≥ τd * (py - P)
@inline ismode2_flowpotential(pq, pd, τd, τII, P) = τII * (pq - pd) ≥ τd * (pq - P)

@inline function compute_F(k, c, py, a, Ry, pd, τd, τII, P)
    F = if ismode2_yield(py, pd, τd, τII, P)
        @muladd τII - k * P - c
    else
        a * (hypot(τII, P - py) - Ry)
    end
    return F
end

@inline function compute_Q(cp, τII, P)
    return Q = if ismode2_flowpotential(cp.pq, cp.pd, cp.τd, τII, P)
        cons = cp.sdf - cp.kq * cp.pdf
        @muladd τII - cp.kq * P - cons
    else
        cons = cp.Rf
        Rq = hypot(τII, (P - cp.pq))
        cp.b * (Rq - cons)
    end
end


"""
    _Aτ_Ap(cp, τII, P)

Flow-potential scalar coefficients from Popov et al. (2025), Eq. 43.

The full tensor gradient of Q is decomposed in Eq. 21-22 as:

    ∂Q/∂σ_ij = B_τ · τ_ij + B_p · δ_ij

where B_τ and B_p are the tensor coefficients. Eq. 43 defines the scalar
invariant-level coefficients A_τ and A_p, which relate to actual derivatives as:

- A_τ = ∂Q/∂τII / 2     (half the actual τII-derivative of Q)
- A_p = -∂Q/∂P          (negated actual P-derivative of Q)

To recover actual derivatives:  ∂Q/∂τII = 2A_τ,  ∂Q/∂P = -A_p.
"""
@inline function _Aτ_Ap(cp, τII, P)
    @unpack pq, pd, τd, b, kq = cp

    if ismode2_flowpotential(pq, pd, τd, τII, P)
        return 0.5, kq
    else
        Rq = hypot(τII, P - pq)
        if iszero(Rq)
            return zero(τII), zero(P)
        end
        inv_Rq = inv(Rq)
        return 0.5 * b * τII * inv_Rq, - b * (P - pq) * inv_Rq
    end
end


# DPCap yield function
function (s::DruckerPragerCap)(;
        P = 0.0, τII = 0.0, Pf = 0.0, EII = 0.0, perturbation_C = 1.0, kwargs...
    )
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    ϕ = s.softening_ϕ(EII, ϕ)
    Ψ = s.softening_Ψ(EII, Ψ)
    C = s.softening_C(EII, C)
    C *= perturbation_C

    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    sinΨ = iszero(EII) ? sinΨ : sind(Ψ)
    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    F = compute_F(cp.k, cp.c, cp.py, cp.a, cp.Ry, cp.pd, cp.τd, τII, P - Pf)
    return F
end

function (s::DruckerPragerCap{_T, U, U1, U2, NoSoftening, NoSoftening, NoSoftening})(;
        P = 0.0, τII = 0.0, Pf = 0.0, EII = 0.0, perturbation_C = 1.0, kwargs...,
    ) where {_T, U, U1, U2}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    C *= perturbation_C

    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    F = compute_F(cp.k, cp.c, cp.py, cp.a, cp.Ry, cp.pd, cp.τd, τII, P - Pf)
    return F
end

function (s::DruckerPragerCap{_T, U, U1, U2, NoSoftening, NoSoftening, AbstractSoftening})(;
        P = 0.0, τII = 0.0, Pf = 0.0, EII = 0.0, perturbation_C = 1.0, kwargs...,
    ) where {_T, U, U1, U2, AbstractSoftening}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    Ψ = s.softening_Ψ(EII, Ψ)
    C *= perturbation_C
    sinΨ = iszero(EII) ? sinΨ : sind(Ψ)

    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    F = compute_F(cp.k, cp.c, cp.py, cp.a, cp.Ry, cp.pd, cp.τd, τII, P - Pf)
    return F
end

function (s::DruckerPragerCap{_T, U, U1, U2, NoSoftening, AbstractSoftening, AbstractSoftening})(;
        P = 0.0, τII = 0.0, Pf = 0.0, EII = 0.0, perturbation_C = 1.0, kwargs...,
    ) where {_T, U, U1, U2, AbstractSoftening}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    C = s.softening_C(EII, C)
    Ψ = s.softening_Ψ(EII, Ψ)
    C *= perturbation_C
    sinΨ = iszero(EII) ? sinΨ : sind(Ψ)

    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    F = compute_F(cp.k, cp.c, cp.py, cp.a, cp.Ry, cp.pd, cp.τd, τII, P - Pf)
    return F
end

function (s::DruckerPragerCap{_T, U, U1, U2, NoSoftening, AbstractSoftening, NoSoftening})(;
        P = 0.0, τII = 0.0, Pf = 0.0, EII = 0.0, perturbation_C = 1.0, kwargs...,
    ) where {_T, U, U1, U2, AbstractSoftening}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    C = s.softening_C(EII, C)
    C *= perturbation_C

    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    F = compute_F(cp.k, cp.c, cp.py, cp.a, cp.Ry, cp.pd, cp.τd, τII, P - Pf)
    return F
end

function (s::DruckerPragerCap{_T, U, U1, U2, AbstractSoftening, AbstractSoftening, NoSoftening})(;
        P = 0.0, τII = 0.0, Pf = 0.0, EII = 0.0, perturbation_C = 1.0, kwargs...,
    ) where {_T, U, U1, U2, AbstractSoftening}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    ϕ = s.softening_ϕ(EII, ϕ)
    C = s.softening_C(EII, C)
    C *= perturbation_C

    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    F = compute_F(cp.k, cp.c, cp.py, cp.a, cp.Ry, cp.pd, cp.τd, τII, P - Pf)
    return F
end

function (s::DruckerPragerCap{_T, U, U1, U2, AbstractSoftening, NoSoftening, NoSoftening})(;
        P = 0.0, τII = 0.0, Pf = 0.0, EII = 0.0, perturbation_C = 1.0, kwargs...,
    ) where {_T, U, U1, U2}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    ϕ = s.softening_ϕ(EII, ϕ)
    C *= perturbation_C

    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    F = compute_F(cp.k, cp.c, cp.py, cp.a, cp.Ry, cp.pd, cp.τd, τII, P - Pf)
    return F
end

function (s::DruckerPragerCap{_T, U, U1, U2, AbstractSoftening, NoSoftening, AbstractSoftening})(;
        P = 0.0, τII = 0.0, Pf = 0.0, EII = 0.0, perturbation_C = 1.0, kwargs...,
    ) where {_T, U, U1, U2}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    ϕ = s.softening_ϕ(EII, ϕ)
    Ψ = s.softening_Ψ(EII, Ψ)
    C *= perturbation_C

    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    sinΨ = iszero(EII) ? sinΨ : sind(Ψ)
    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    F = compute_F(cp.k, cp.c, cp.py, cp.a, cp.Ry, cp.pd, cp.τd, τII, P - Pf)
    return F
end
"""
    compute_yieldfunction(s::DruckerPragerCap; P, τII, Pf, kwargs...)

Computes the plastic yield function `F` for a given second invariant of the deviatoric stress tensor `τII`,  `P` pressure, and `Pf` fluid pressure.
"""
function compute_yieldfunction(
        s::DruckerPragerCap;
        P = 0.0,
        τII = 0.0,
        Pf = 0.0,
        EII = 0.0,
        perturbation_C = 1.0,
        kwargs...,
    )
    return s(; P = P, τII = τII, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
end

"""
    compute_yieldfunction!(F, s::DruckerPragerCap; P, τII, Pf, kwargs...)

Computes the plastic yield function `F` for a given second invariant of the deviatoric stress tensor `τII`,  `P` pressure, and `Pf` fluid pressure, and stores the result in `F`.
"""
function compute_yieldfunction!(
        F::AbstractArray{_T, N},
        s::DruckerPragerCap;
        P::AbstractArray{_T, N},
        τII::AbstractArray{_T, N},
        Pf::AbstractArray{_T, N} = zero(P),
        EII::AbstractArray{_T, N} = zero(P),
        perturbation_C = 1.0,
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(P)
        F[i] = compute_yieldfunction(s; P = P[i], τII = τII[i], Pf = Pf[i], EII = EII[i], perturbation_C = perturbation_C)
    end
    return nothing
end

function compute_flowpotential(
        s::DruckerPragerCap;
        P = 0.0,
        τII = 0.0,
        Pf = 0.0,
        EII = 0.0,
        perturbation_C = 1.0,
        kwargs...,
    )
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    ϕ = s.softening_ϕ(EII, ϕ)
    Ψ = s.softening_Ψ(EII, Ψ)
    C = s.softening_C(EII, C)
    C *= perturbation_C
    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    sinΨ = iszero(EII) ? sinΨ : sind(Ψ)
    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    Q = compute_Q(cp, τII, P - Pf)
    return Q
end

@inline compute_flowpotential(s::DruckerPragerCap, args) = compute_flowpotential(s; args...)

function compute_flowpotential!(
        Q::AbstractArray{_T, N},
        s::DruckerPragerCap;
        P::AbstractArray{_T, N},
        τII::AbstractArray{_T, N},
        Pf::AbstractArray{_T, N} = zero(P),
        EII::AbstractArray{_T, N} = zero(P),
        perturbation_C = 1.0,
        kwargs...,
    ) where {N, _T}
    @inbounds for i in eachindex(P)
        Q[i] = compute_flowpotential(s; P = P[i], τII = τII[i], Pf = Pf[i], EII = EII[i], perturbation_C = perturbation_C)
    end
    return nothing
end

# ∂Q∂P: actual derivative ∂Q/∂P = -Ap
function ∂Q∂P(
        s::DruckerPragerCap,
        P::_T;
        τII = zero(_T),
        Pf = zero(_T),
        EII = zero(_T),
        perturbation_C = one(_T),
        kwargs...,
    ) where {_T}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    ϕ = s.softening_ϕ(EII, ϕ)
    Ψ = s.softening_Ψ(EII, Ψ)
    C = s.softening_C(EII, C)
    C *= perturbation_C

    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    sinΨ = iszero(EII) ? sinΨ : sind(Ψ)

    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)
    _, Ap = _Aτ_Ap(cp, τII, P - Pf)
    return -Ap  # ∂Q/∂P = -Ap
end


function ∂Q∂τII(
        s::DruckerPragerCap,
        τII::_T;
        P = zero(_T),
        Pf = zero(_T),
        EII = zero(_T),
        perturbation_C = one(_T),
        kwargs...,
    ) where {_T}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    ϕ = s.softening_ϕ(EII, ϕ)
    Ψ = s.softening_Ψ(EII, Ψ)
    C = s.softening_C(EII, C)
    C *= perturbation_C
    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    sinΨ = iszero(EII) ? sinΨ : sind(Ψ)

    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)
    Aτ, _ = _Aτ_Ap(cp, τII, P - Pf)
    return Aτ
end

function ∂F∂τII(
        s::DruckerPragerCap,
        τII::_T;
        P = zero(_T),
        Pf = zero(_T),
        EII = zero(_T),
        perturbation_C = one(_T),
        kwargs...,
    ) where {_T}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    ϕ = s.softening_ϕ(EII, ϕ)
    Ψ = s.softening_Ψ(EII, Ψ)
    C = s.softening_C(EII, C)
    C *= perturbation_C
    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    sinΨ = iszero(EII) ? sinΨ : sind(Ψ)
    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    if ismode2_yield(cp.py, cp.pd, cp.τd, τII, P - Pf)
        return one(_T)
    else
        Ry = hypot(τII, (P - Pf) - cp.py)
        return iszero(Ry) ? zero(_T) : cp.a * τII * inv(Ry)
    end
end

function ∂F∂P(
        s::DruckerPragerCap,
        P::_T;
        τII = zero(_T),
        Pf = zero(_T),
        EII = zero(_T),
        perturbation_C = one(_T),
        kwargs...,
    ) where {_T}
    @unpack_val sinϕ, cosϕ, sinΨ, ϕ, Ψ, C, pT = s
    ϕ = s.softening_ϕ(EII, ϕ)
    Ψ = s.softening_Ψ(EII, Ψ)
    C = s.softening_C(EII, C)
    C *= perturbation_C
    sinϕ, cosϕ = iszero(EII) ? (sinϕ, cosϕ) : sincosd(ϕ)
    sinΨ = iszero(EII) ? sinΨ : sind(Ψ)

    cp = compute_tensile_cap(sinϕ, cosϕ, sinΨ, C, pT)

    if ismode2_yield(cp.py, cp.pd, cp.τd, τII, (P - Pf))
        return -cp.k
    else
        Ry = hypot(τII, (P - Pf) - cp.py)
        return iszero(Ry) ? zero(_T) : cp.a * ((P - Pf) - cp.py) * inv(Ry)
    end
end

∂F∂λ(s::DruckerPragerCap, τII::_T; P = zero(_T), kwargs...) where {_T} = zero(_T)

# Component gradient functions for DruckerPragerCap
#
# These return actual derivatives ∂Q/∂τij, computed from the Popov Eq. 43
# coefficients Aτ = ∂Q/∂τII / 2 via the chain rule:
#
#   ∂Q/∂τij = ∂Q/∂τII · ∂τII/∂τij = 2Aτ · ∂τII/∂τij
#
# where ∂τII/∂τij = 0.5·τij/τII  for diagonal (xx,yy,zz)
#       ∂τII/∂τij =     τij/τII  for shear    (yz,xz,xy)
#
# Result:  diagonal → Aτ · τij / τII
#          shear    → 2Aτ · τij / τII

for t in (:NTuple, :SVector)
    @eval begin
        # --- 3D (6-component Voigt) ---
        # diagonal components: ∂Q/∂τij = Aτ * τij / τII
        function ∂Q∂τxx(p::DruckerPragerCap, τij::$(t){6, T}; P = zero(T), Pf = zero(T), EII = zero(T), perturbation_C = one(T), kwargs...) where {T}
            τII = second_invariant(τij)
            Aτ = ∂Q∂τII(p, τII; P = P, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
            return iszero(τII) ? zero(T) : Aτ * τij[1] / τII
        end
        function ∂Q∂τyy(p::DruckerPragerCap, τij::$(t){6, T}; P = zero(T), Pf = zero(T), EII = zero(T), perturbation_C = one(T), kwargs...) where {T}
            τII = second_invariant(τij)
            Aτ = ∂Q∂τII(p, τII; P = P, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
            return iszero(τII) ? zero(T) : Aτ * τij[2] / τII
        end
        function ∂Q∂τzz(p::DruckerPragerCap, τij::$(t){6, T}; P = zero(T), Pf = zero(T), EII = zero(T), perturbation_C = one(T), kwargs...) where {T}
            τII = second_invariant(τij)
            Aτ = ∂Q∂τII(p, τII; P = P, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
            return iszero(τII) ? zero(T) : Aτ * τij[3] / τII
        end
        # shear components: ∂Q/∂τij = 2Aτ * τij / τII
        function ∂Q∂τyz(p::DruckerPragerCap, τij::$(t){6, T}; P = zero(T), Pf = zero(T), EII = zero(T), perturbation_C = one(T), kwargs...) where {T}
            τII = second_invariant(τij)
            Aτ = ∂Q∂τII(p, τII; P = P, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
            return iszero(τII) ? zero(T) : 2 * Aτ * τij[4] / τII
        end
        function ∂Q∂τxz(p::DruckerPragerCap, τij::$(t){6, T}; P = zero(T), Pf = zero(T), EII = zero(T), perturbation_C = one(T), kwargs...) where {T}
            τII = second_invariant(τij)
            Aτ = ∂Q∂τII(p, τII; P = P, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
            return iszero(τII) ? zero(T) : 2 * Aτ * τij[5] / τII
        end
        function ∂Q∂τxy(p::DruckerPragerCap, τij::$(t){6, T}; P = zero(T), Pf = zero(T), EII = zero(T), perturbation_C = one(T), kwargs...) where {T}
            τII = second_invariant(τij)
            Aτ = ∂Q∂τII(p, τII; P = P, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
            return iszero(τII) ? zero(T) : 2 * Aτ * τij[6] / τII
        end

        # --- 2D (3-component Voigt) ---
        # diagonal components: ∂Q/∂τij = Aτ * τij / τII
        function ∂Q∂τxx(p::DruckerPragerCap, τij::$(t){3, T}; P = zero(T), Pf = zero(T), EII = zero(T), perturbation_C = one(T), kwargs...) where {T}
            τII = second_invariant(τij)
            Aτ = ∂Q∂τII(p, τII; P = P, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
            return iszero(τII) ? zero(T) : Aτ * τij[1] / τII
        end
        function ∂Q∂τyy(p::DruckerPragerCap, τij::$(t){3, T}; P = zero(T), Pf = zero(T), EII = zero(T), perturbation_C = one(T), kwargs...) where {T}
            τII = second_invariant(τij)
            Aτ = ∂Q∂τII(p, τII; P = P, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
            return iszero(τII) ? zero(T) : Aτ * τij[2] / τII
        end
        # shear component: ∂Q/∂τij = 2Aτ * τij / τII
        function ∂Q∂τxy(p::DruckerPragerCap, τij::$(t){3, T}; P = zero(T), Pf = zero(T), EII = zero(T), perturbation_C = one(T), kwargs...) where {T}
            τII = second_invariant(τij)
            Aτ = ∂Q∂τII(p, τII; P = P, Pf = Pf, EII = EII, perturbation_C = perturbation_C)
            return iszero(τII) ? zero(T) : 2 * Aτ * τij[3] / τII
        end
    end
end

function show(io::IO, s::DruckerPragerCap)
    return print(io, "DruckerPragerCap(ϕ=$(UnitValue(s.ϕ)), Ψ=$(UnitValue(s.Ψ)), C=$(UnitValue(s.C)), η_vp=$(UnitValue(s.η_vp)), pT=$(UnitValue(s.pT))")
end
