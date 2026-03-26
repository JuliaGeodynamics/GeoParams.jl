# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)
abstract type AbstractPlasticity{T} <: AbstractConstitutiveLaw{T} end
abstract type AbstractPlasticPotential{Float64} <: AbstractConstitutiveLaw{Float64} end

export AbstractPlasticity,
    isvolumetric,
    compute_yieldfunction, # calculation routines
    compute_yieldfunction!,
    compute_plasticpotentialDerivative,
    ∂Q∂τ, ∂Q∂τII, ∂Q∂P,
    ∂F∂τII, ∂F∂P, ∂F∂λ,
    compute_εII

include("DruckerPrager.jl")    # DP plasticity
include("DruckerPrager_regularised.jl")    # regularized DP plasticity
include("DruckerPragerCap.jl")    # DP plasticity with tensile cap


# Thin convenience wrappers
# 3D
function ∂Q∂τ(p::AbstractPlasticity, τij::SVector{6}; kwargs...)
    return @SVector [∂Q∂τxx(p, τij), ∂Q∂τyy(p, τij), ∂Q∂τzz(p, τij), ∂Q∂τyz(p, τij), ∂Q∂τxz(p, τij), ∂Q∂τxy(p, τij)]
end

function ∂Q∂τ(p::AbstractPlasticity, τij::NTuple{6}; kwargs...)
    return ∂Q∂τxx(p, τij), ∂Q∂τyy(p, τij), ∂Q∂τzz(p, τij), ∂Q∂τyz(p, τij), ∂Q∂τxz(p, τij), ∂Q∂τxy(p, τij)
end

# 2D
function ∂Q∂τ(p::AbstractPlasticity, τij::SVector{3}; kwargs...)
    return @SVector [∂Q∂τxx(p, τij), ∂Q∂τyy(p, τij), ∂Q∂τxy(p, τij)]
end

function ∂Q∂τ(p::AbstractPlasticity, τij::NTuple{3}; kwargs...)
    return ∂Q∂τxx(p, τij), ∂Q∂τyy(p, τij), ∂Q∂τxy(p, τij)
end

# Compute partial derivatives of a generic user-defined Q using AD
∂Q∂τ(Q::F, args::SVector; kwargs...) where {F <: Function} = ForwardDiff.gradient(Q, args)
∂Q∂τ(Q::F, args::Vector; kwargs...) where {F <: Function} = ForwardDiff.gradient(Q, args)
function ∂Q∂τ(Q::F, args::NTuple{N}; kwargs...) where {N, F <: Function}
    tmp = ∂Q∂τ(Q, SVector{N}(args...))
    return ntuple(i -> tmp[i], Val(N))
end

# Wrapper for arbitrary args in the form of a NamedTuple
function ∂Q∂τ(p::AbstractPlasticity, args::NamedTuple; kwargs...)
    return ∂Q∂τ(Q, args.τij, kwargs...)
end
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------

# Plastic finite strain and strain rate

"""
    plastic_strain(εvp::T, p::AbstractPlasticity{T}, τij, λ̇::T, dt::T)
    
    Integrate the finite plastic strain. Equations from Duretz et al. 2019 G3
"""
function plastic_strain(εvp::T, p::AbstractPlasticity{T}, τij, λ̇::T, dt::T) where {T}
    return εvp += plastic_strain(p, τij, λ̇) * dt
end

@inline function plastic_strain(p::AbstractPlasticity{T}, τij::T, λ̇::T) where {T}
    εvp_ij = plastic_strain_rate(p, τij, λ̇)
    εvp = √((2.0 / 3.0) * dot(εvp_ij, εvp_ij))
    return εvp
end

@inline plastic_strain_rate(p::AbstractPlasticity{T}, τij::T, λ̇::T) where {T} = ∂Q∂τ(p, τij) .* λ̇
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
for myType in (:DruckerPrager, :DruckerPrager_regularised, :DruckerPragerCap)
    @eval begin
        (p::$(myType))(args) = p(; args...)
        ∂Q∂τ(p::$(myType), args, kwargs) = ∂Q∂τ(p, args; kwargs...)
        ∂Q∂P(p::$(myType), args, kwargs) = ∂Q∂P(p, args; kwargs...)
        ∂Q∂τII(p::$(myType), args, kwargs) = ∂Q∂τII(p, args; kwargs...)

        ∂F∂P(p::$(myType), args, kwargs) = ∂F∂P(p, args; kwargs...)
        ∂F∂λ(p::$(myType), args, kwargs) = ∂F∂λ(p, args; kwargs...)
        ∂F∂τII(p::$(myType), args, kwargs) = ∂F∂τII(p, args; kwargs...)

        compute_yieldfunction(p::$(myType), args) = p(args)
        compute_εII(p::$(myType), args) = compute_εII(p, args...)

        function compute_yieldfunction!(
                H::AbstractArray{_T, N}, p::$(myType){_T}, args
            ) where {_T, N}
            return compute_yieldfunction!(H, p; args...)
        end
    end
end

compute_yieldfunction(args...) = compute_param(compute_yieldfunction, args...)
compute_yieldfunction!(args...) = compute_param!(compute_yieldfunction, args...)
compute_plasticpotentialDerivative(args...) = compute_param(∂Q∂τ, args...)
∂Q∂τ(p::AbstractMaterialParamsStruct, args) = compute_plasticpotentialDerivative(p, args)
∂Q∂τ(args...) = compute_param(∂Q∂τ, args...)
∂Q∂τII(args...) = compute_param(∂Q∂τII, args...)

function compute_plasticpotentialDerivative(p::AbstractMaterialParamsStruct, args)
    return ∂Q∂τ(p.Plasticity[1], args)
end

∂Q∂P(args...) = compute_param(∂Q∂P, args...)

function ∂Q∂P(p::AbstractMaterialParamsStruct, args)
    return ∂Q∂P(p.Plasticity[1], args)
end

lambda(args...) = compute_param(lambda, args...)
plastic_strain_rate(args...) = compute_param(plastic_strain_rate, args...)
plastic_strain(args...) = compute_param(plastic_strain, args...)
