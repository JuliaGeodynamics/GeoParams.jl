export CustomRheology, dεII_dτII, dτII_dεII, compute_εII, compute_τII

struct CustomRheology{F1,F2,T} <: AbstractConstitutiveLaw{Float64}
    strain::F1 # function to compute strain rate
    stress::F2 # function to compute deviatoric stress
    args::T # NamedTuple of parameters

    function CustomRheology(strain::F1, stress::F2, args::T) where {F1,F2,T}
        f_strain = has_kwargs(strain) ? strain : (a, τII; kwargs...) -> a.strain(a, τII)
        f_stress = has_kwargs(stress) ? stress : (a, εII; kwargs...) -> a.stress(a, εII)
        return new{typeof(f_strain),typeof(f_stress),T}(f_strain, f_stress, args)
    end
end

has_kwargs(f) = first(methods(f)).nkw == 0

@inline compute_εII(a::CustomRheology, τII, args) = a.strain(a, τII; args...)
@inline compute_τII(a::CustomRheology, εII, args) = a.stress(a, εII; args...)

function compute_εII!(εII::AbstractArray, a::CustomRheology, τII::AbstractArray, args)
    for i in eachindex(εII)
        @inbounds εII[i] = a.strain(a, τII[i]; ntuple_idx(args, i)...)
    end
end

function compute_τII!(τII::AbstractArray, a::CustomRheology, εII::AbstractArray, args)
    for i in eachindex(εII)
        @inbounds τII[i] = a.stress(a, εII[i]; ntuple_idx(args, i)...)
    end
end

@inline function dεII_dτII(a::CustomRheology, τII, args)
    return ForwardDiff.derivative(τII -> a.strain(a, τII; args...), τII)
end

@inline function dτII_dεII(a::CustomRheology, εII, args)
    return ForwardDiff.derivative(εII -> a.stress(a, εII; args...), εII)
end
