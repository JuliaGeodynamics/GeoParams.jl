export CustomRheology,
    dεII_dτII,
    dτII_dεII,
    compute_εII,
    compute_τII

struct CustomRheology{F1, F2, T} <: AbstractCreepLaw{Float64}
    strain::F1 # function to compute strain rate
    stress::F2 # function to compute deviatoric stress
    args::T # NamedTuple of parameters
end

function compute_εII(a::CustomRheology, TauII, args)
    a.strain(a, TauII; args...)
end

function compute_εII!(EpsII::AbstractArray, a::CustomRheology, TauII::AbstractArray, args)
    for i in eachindex(EpsII)
        @inbounds EpsII[i] = a.strain(a, TauII[i]; ntuple_idx(args)...)
    end
end

function compute_τII(a::CustomRheology, EpsII, args)
    a.stress(a, EpsII; args...)
end

function compute_τII!(TauII::AbstractArray, a::CustomRheology, EpsII::AbstractArray, args)
    for i in eachindex(EpsII)
        @inbounds TauII[i] = a.stress(a, EpsII[i]; ntuple_idx(args)...)
    end
end

function dεII_dτII(a::CustomRheology, TauII, args)
    ForwardDiff.derivative(TauII -> a.strain(a, TauII; args...), TauII)
end

function dτII_dεII(a::CustomRheology, EpsII, args)
    ForwardDiff.derivative(EpsII -> a.stress(a, EpsII; args...), EpsII)
end