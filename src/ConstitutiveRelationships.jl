module ConstitutiveRelationships

# This implements constitutitive relationships (typically between stress & strainrate)
#

using Base: Float64
using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import GeoParams: param_info, fastpow
import GeoParams: second_invariant, second_invariant_staggered
using BibTeX
using ..MaterialParameters: MaterialParamsInfo
import Base.show
using ForwardDiff
using StaticArrays

const AxialCompression, SimpleShear, Invariant = 1, 2, 3

abstract type AbstractConstitutiveLaw{T} <: AbstractMaterialParam end

export param_info,
    dεII_dτII,
    dτII_dεII,
    compute_εII!,
    compute_εII,
    compute_τII!,
    compute_τII,
    computeViscosity_τII,
    computeViscosity_εII,
    computeViscosity_τII!,
    computeViscosity_εII!,
    local_iterations_εII,
    computeViscosity,
    strain_rate_circuit,
    InverseCreepLaw,
    KelvinVoigt

include("Computations.jl")

include("CreepLaw/CreepLaw.jl")              # viscous Creeplaws
include("Elasticity/Elasticity.jl")          # elasticity
include("Plasticity/Plasticity.jl")          # plasticity
include("CreepLaw/Viscosity.jl")             # composite creeplaws

# add methods programatically 
for myType in (:LinearViscous, :DiffusionCreep, :DislocationCreep, :ConstantElasticity)
    @eval begin
        compute_εII(a::$(myType), TauII, args) = compute_εII(a, TauII; args...)
        function compute_εII!(
            ε::AbstractArray{_T,N}, s::$(myType){_T}, TauII::AbstractArray{_T,N}, args
        ) where {_T,N}
            return compute_εII!(ε, s, TauII; args...)
        end

        compute_τII(a::$(myType), EpsII, args) = compute_τII(a, EpsII; args...)
        function compute_τII!(
            τ::AbstractArray{_T,N}, s::$(myType){_T}, EpsII::AbstractArray{_T,N}, args
        ) where {_T,N}
            return compute_τII!(τ, s, EpsII; args...)
        end

        if Symbol($myType) !== :ConstantElasticity
            dτII_dεII(a::$(myType), EpsII, args) = dτII_dεII(a, EpsII; args...)
            dεII_dτII(a::$(myType), TauII, args) = dεII_dτII(a, TauII; args...)
        end
    end
end

end
