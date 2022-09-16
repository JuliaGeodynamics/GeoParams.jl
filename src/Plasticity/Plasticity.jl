# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)
abstract type AbstractPlasticity{T} <: AbstractConstitutiveLaw{T} end

export compute_yieldfunction,      # calculation routines
    compute_yieldfunction!,
    DruckerPrager,               # constant
    PlasticFlow,
    ∂Q∂τ

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
where ``\\dot{\\lambda}`` is a (scalar) that is nonzero and chosen such that the resuling stress gives ``F(\\tau_{II}^{final})=0``, and ``\\sigma_{ij}=-P + \\tau_{ij}`` denotes the total stress tensor.   
        
"""
@with_kw_noshow struct DruckerPrager{T, U, U1} <: AbstractPlasticity{T}
    ϕ::GeoUnit{T,U} = 30NoUnits # Friction angle
    Ψ::GeoUnit{T,U} = 0NoUnits # Dilation angle
    C::GeoUnit{T,U1} = 10e6Pa # Cohesion
end
DruckerPrager(args...) = DruckerPrager(convert.(GeoUnit, args[1:3])..., args[4:end]...)

struct PlasticFlow{F1,F2,F3,F4,F5,F6,F7}
    Q::F1 # Plastic potential
    ∂Q∂τxx::F2
    ∂Q∂τyy::F3
    ∂Q∂τzz::F4
    ∂Q∂τyz::F5
    ∂Q∂τxz::F6
    ∂Q∂τxy::F7

    function PlasticFlow(;
        Q::F1=second_invariant,
        ∂Q∂τxx::F2=(τij) -> 0.5 * τij[1] / second_invariant(τij),
        ∂Q∂τyy::F3=(τij) -> 0.5 * τij[2] / second_invariant(τij),
        ∂Q∂τzz::F4=nothing,
        ∂Q∂τyz::F5=nothing,
        ∂Q∂τxz::F6=nothing,
        ∂Q∂τxy::F7=(τij) -> τij[3] / second_invariant(τij),
    ) where {F1,F2,F3,F4,F5,F6,F7}
        return new{F1,F2,F3,F4,F5,F6,F7}(
           Q, ∂Q∂τxx, ∂Q∂τyy, ∂Q∂τzz, ∂Q∂τyz, ∂Q∂τxz, ∂Q∂τxy
        )
    end
end

# struct DruckerPrager{T,F1,F2,F3,F4,F5,F6,F7,U,U1} <: AbstractPlasticity{T}
#     ϕ::GeoUnit{T,U} # Friction angle
#     Ψ::GeoUnit{T,U} # Dilation angle
#     C::GeoUnit{T,U1} # Cohesion
#     Q::GeoUnit{F1, typeof(NoUnits)} # Plastic potential
#     ∂Q∂τxx::GeoUnit{F2, typeof(NoUnits)}
#     ∂Q∂τyy::GeoUnit{F3, typeof(NoUnits)}
#     ∂Q∂τzz::GeoUnit{F4, typeof(NoUnits)}
#     ∂Q∂τyz::GeoUnit{F5, typeof(NoUnits)}
#     ∂Q∂τxz::GeoUnit{F6, typeof(NoUnits)}
#     ∂Q∂τxy::GeoUnit{F7, typeof(NoUnits)}

#     function DruckerPrager(;
#         ϕ::A=30NoUnits, # Friction angle
#         Ψ::A=0NoUnits, # Dilation angle
#         C::B=10e6Pa, # Cohesion
#         Q::F1=second_invariant,
#         ∂Q∂τxx::F2=(τij) -> 0.5 * τij[1] / second_invariant(τij),
#         ∂Q∂τyy::F3=(τij) -> 0.5 * τij[2] / second_invariant(τij),
#         ∂Q∂τzz::F4=nothing,
#         ∂Q∂τyz::F5=nothing,
#         ∂Q∂τxz::F6=nothing,
#         ∂Q∂τxy::F7=(τij) -> τij[3] / second_invariant(τij),
#     ) where {A,B,F1,F2,F3,F4,F5,F6,F7}
#         ϕU = ϕ isa GeoUnit ? ϕ : convert(GeoUnit, ϕ)
#         ΨU = Ψ isa GeoUnit ? Ψ : convert(GeoUnit, Ψ)
#         CU = C isa GeoUnit ? C : convert(GeoUnit, C)

#         funs = (Q, ∂Q∂τxx, ∂Q∂τyy, ∂Q∂τzz, ∂Q∂τyz, ∂Q∂τxz, ∂Q∂τxy)
#         funsU = ntuple(Val(7)) do i 
#             GeoUnit(funs[i])
#         end

#         t1 = typeof(ϕU)
#         t2 = typeof(CU)
#         return new{t1.types[1],F1,F2,F3,F4,F5,F6,F7,t1.types[2],t2.types[2],}(
#             ϕU, ΨU, CU, funsU...
#         )
#     end
# end

function param_info(s::DruckerPrager) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"F = \\tau_{II} - \\cos(ϕ)C - \\sin(ϕ)(P-P_f); Q=\\tau_{II} - \\sin(Ψ)(P-P_f)",
    )
end

# Calculation routines
function (s::DruckerPrager{_T,U,U1})(;
    P::_T=zero(_T), τII::_T=zero(_T), Pf::_T=zero(_T), kwargs...
) where {_T,U,U1}
    @unpack_val ϕ, C = s
    sinϕ, cosϕ = sincosd(ϕ)

    F = τII - cosϕ * C - sinϕ * (P - Pf)   # with fluid pressure (set to zero by default)

    return F
end

"""
    compute_yieldfunction(s::DruckerPrager; P, τII_old, Pf, kwargs...) 

Computes the plastic yield function `F` for a given second invariant of the deviatoric stress tensor `τII`,  `P` the pressure, and `Pf` fluid pressure.
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

# Print info 
function show(io::IO, g::DruckerPrager)
    return print(
        io,
        "Drucker-Prager plasticity with: C = $(UnitValue(g.C)), ϕ = $(UnitValue(g.ϕ))ᵒ, Ψ = $(UnitValue(g.Ψ))ᵒ",
    )
end
#-------------------------------------------------------------------------

# Plastic Potential derivatives

∂Q∂τij(Q::PlasticFlow, τij::SVector{N, T}) where {N, T} = ForwardDiff.gradient(Q.Q, τij)
∂Q∂τij(Q::PlasticFlow, τij::Vector{T}) where {T} = ForwardDiff.gradient(Q.Q, τij)
function ∂Q∂τij(Q::PlasticFlow, τij::NTuple{N,T}) where {N,T}
    return ∂Q∂τij(Q, SVector{N}(τij...))
end

function ∂Q∂τ(
    Q::PlasticFlow{F1,Nothing,Nothing,F4,F5,F6,Nothing}, τij::SVector{3,T}
) where {T,F1,F4,F5,F6}
    return ∂Q∂τij(Q, τij)
end

function ∂Q∂τ(
    Q::PlasticFlow{F1,F2,F3,F4,F5,F6,F7}, τij::SVector{3,T}
) where {T,F1,F2,F3,F4,F5,F6,F7}
    return SVector{3,T}(Q.∂Q∂τxx(τij), Q.∂Q∂τyy(τij), Q.∂Q∂τxy(τij))
end

# Wrapper for NTuple inputs (NTuples not supported by ForwardDiff.jl, but @SVectors are)
function ∂Q∂τ(Q::PlasticFlow, τij::NTuple{N,T}) where {N,T} 
    tmp = ∂Q∂τij(Q, SVector{N,T}(τij))
    return ntuple(i->tmp[i], Val(N))
end
#-------------------------------------------------------------------------

# Computational routines needed for computations with the MaterialParams structure 
function compute_yieldfunction(s::AbstractMaterialParamsStruct, args)
    if isempty(s.Plasticity)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return s.Plasticity[1](args)
    end
end

# add methods programmatically
for myType in (:DruckerPrager,)
    @eval begin
        (s::$(myType))(args) = s(; args...)
        compute_yieldfunction(s::$(myType), args) = s(args)
        function compute_yieldfunction!(
            H::AbstractArray{_T,N}, s::$(myType){_T}, args
        ) where {_T,N}
            return compute_yieldfunction!(H, s; args...)
        end
    end
end

compute_yieldfunction(args...) = compute_param(compute_yieldfunction, args...)
compute_yieldfunction!(args...) = compute_param!(compute_yieldfunction, args...)
