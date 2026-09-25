# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)
"""
    AbstractPlasticity{T} <: AbstractConstitutiveLaw{T}

Supertype of plastic constitutive laws such as [`DruckerPrager`](@ref), which define a yield
function and a plastic flow potential.
"""
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

"""
    ∂Q∂τ(p::AbstractPlasticity, τij; kwargs...)

Returns the gradient of the plastic flow potential `Q` with respect to the deviatoric stress tensor
`τij` (given as a 3- or 6-component `NTuple`/`SVector`), i.e. the direction of plastic flow
``\\partial Q/\\partial \\tau_{ij}``. Also accessible as [`compute_plasticpotentialDerivative`](@ref).

Keyword arguments (e.g. `P`, `Pf`, `EII`, `perturbation_C`) are forwarded to the component functions
and define the state at which the gradient is evaluated.
"""
function ∂Q∂τ end

"""
    ∂Q∂τII(p::AbstractPlasticity, τII; kwargs...)

Returns the derivative of the plastic flow potential `Q` with respect to the second invariant of the
deviatoric stress, ``\\partial Q/\\partial \\tau_{II}``.

For [`DruckerPragerCap`](@ref) the returned value is `Aτ = (∂Q/∂τII)/2`: the diagonal components of
[`∂Q∂τ`](@ref) are `Aτ * τij / τII` and the shear components are `2Aτ * τij / τII`
(engineering-shear convention), so `∂Q∂τII` is half the full derivative.
"""
function ∂Q∂τII end

"""
    ∂Q∂P(p::AbstractPlasticity, args; kwargs...)

Returns the derivative of the plastic flow potential `Q` with respect to pressure,
``\\partial Q/\\partial P``, which controls the dilatancy of the plastic flow.
"""
function ∂Q∂P end

"""
    ∂F∂τII(p::AbstractPlasticity, args; kwargs...)

Returns the derivative of the yield function `F` with respect to the second invariant of the
deviatoric stress, ``\\partial F/\\partial \\tau_{II}``.
"""
function ∂F∂τII end

"""
    ∂F∂P(p::AbstractPlasticity, args; kwargs...)

Returns the derivative of the yield function `F` with respect to pressure, ``\\partial F/\\partial P``.
"""
function ∂F∂P end

"""
    ∂F∂λ(p::AbstractPlasticity, args; kwargs...)

Returns the derivative of the yield function `F` with respect to the plastic multiplier `λ`, used in
the return-mapping iteration.
"""
function ∂F∂λ end

"""
    compute_plasticpotentialDerivative(p, args)

Returns the gradient of the plastic flow potential with respect to the deviatoric stress tensor for
the plasticity law or `MaterialParams` `p`; equivalent to [`∂Q∂τ`](@ref).
"""
function compute_plasticpotentialDerivative end

# Thin convenience wrappers
# 3D
function ∂Q∂τ(p::AbstractPlasticity, τij::SVector{6}; kwargs...)
    return @SVector [∂Q∂τxx(p, τij; kwargs...), ∂Q∂τyy(p, τij; kwargs...), ∂Q∂τzz(p, τij; kwargs...), ∂Q∂τyz(p, τij; kwargs...), ∂Q∂τxz(p, τij; kwargs...), ∂Q∂τxy(p, τij; kwargs...)]
end

function ∂Q∂τ(p::AbstractPlasticity, τij::NTuple{6}; kwargs...)
    return ∂Q∂τxx(p, τij; kwargs...), ∂Q∂τyy(p, τij; kwargs...), ∂Q∂τzz(p, τij; kwargs...), ∂Q∂τyz(p, τij; kwargs...), ∂Q∂τxz(p, τij; kwargs...), ∂Q∂τxy(p, τij; kwargs...)
end

# 2D
function ∂Q∂τ(p::AbstractPlasticity, τij::SVector{3}; kwargs...)
    return @SVector [∂Q∂τxx(p, τij; kwargs...), ∂Q∂τyy(p, τij; kwargs...), ∂Q∂τxy(p, τij; kwargs...)]
end

function ∂Q∂τ(p::AbstractPlasticity, τij::NTuple{3}; kwargs...)
    return ∂Q∂τxx(p, τij; kwargs...), ∂Q∂τyy(p, τij; kwargs...), ∂Q∂τxy(p, τij; kwargs...)
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
    return ∂Q∂τ(p, args.τij; kwargs...)
end
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------

# Plastic finite strain and strain rate

"""
    plastic_strain(εvp::T, p::AbstractPlasticity{T}, τij, λ̇::T, dt::T)
    
    Integrate the finite plastic strain. Equations from Duretz et al. 2019 G3
"""
function plastic_strain(εvp::T, p::AbstractPlasticity{T}, τij, λ̇::T, dt::T; kwargs...) where {T}
    return εvp += plastic_strain(p, τij, λ̇; kwargs...) * dt
end

@inline function plastic_strain(p::AbstractPlasticity{T}, τij, λ̇::T; kwargs...) where {T}
    εvp_ij = plastic_strain_rate(p, τij, λ̇; kwargs...)
    εvp = √((2.0 / 3.0) * dot(εvp_ij, εvp_ij))
    return εvp
end

@inline plastic_strain_rate(p::AbstractPlasticity{T}, τij, λ̇::T; kwargs...) where {T} = ∂Q∂τ(p, τij; kwargs...) .* λ̇
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
