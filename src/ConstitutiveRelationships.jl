module ConstitutiveRelationships

# This implements constitutitive relationships (typically between stress & strainrate)
#

using Base: Float64
using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, AbstractConstitutiveLaw, AbstractComposite
import GeoParams: param_info, fastpow
import GeoParams: second_invariant, second_invariant_staggered
using BibTeX
using ..MaterialParameters: MaterialParamsInfo
import Base.show
using ForwardDiff
using StaticArrays

const AxialCompression, SimpleShear, Invariant = 1, 2, 3

#abstract type AbstractConstitutiveLaw{T} <: AbstractMaterialParam end
#abstract type AbstractComposite <: AbstractMaterialParam end

precision(v::AbstractConstitutiveLaw{T}) where T = T



include("Computations.jl")

include("CreepLaw/CreepLaw.jl")              # viscous Creeplaws
include("Elasticity/Elasticity.jl")          # elasticity
include("Plasticity/Plasticity.jl")          # plasticity
#include("CreepLaw/Viscosity.jl")             # composite creeplaws
include("CompositeRheologies.jl")            # composite constitutive relationships

export param_info,
    dεII_dτII,
    dτII_dεII,
    dεII_dτII_AD,
    dτII_dεII_AD,
    compute_εII!,
    compute_εII,
    compute_εII_AD,
    compute_τII!,
    compute_τII,
    compute_τII_AD,
    computeViscosity_τII,
    computeViscosity_τII!,
    computeViscosity_εII,
    computeViscosity_εII!,
    local_iterations_εII,
    local_iterations_εII_AD,
    local_iterations_τII,
    local_iterations_τII_AD,
    computeViscosity,
    strain_rate_circuit,
    stress_circuit,
    InverseCreepLaw,
    KelvinVoigt,
    Parallel,
    CompositeRheology,
    AbstractComposite,
    AbstractConstitutiveLaw
    AxialCompression, SimpleShear, Invariant 

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

        # Expand derivatives
        dτII_dεII(a::$(myType), EpsII, args) = dτII_dεII(a, EpsII; args...)
        dεII_dτII(a::$(myType), TauII, args) = dεII_dτII(a, TauII; args...)
        
    
    end
end

end
