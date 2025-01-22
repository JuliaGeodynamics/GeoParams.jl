using Test, GeoParams

## HOW TO DEFINE CUSTOM RHEOLOGY

# function to compute strain rate (compulsory)
@inline function custom_εII(a::CustomRheology, TauII; args...)
    η = custom_viscosity(a; args...)
    return (TauII / η) * 0.5
end

# function to compute deviatoric stress (compulsory)
@inline function custom_τII(a::CustomRheology, EpsII; args...)
    η = custom_viscosity(a; args...)
    return 2.0 * (η * EpsII)
end

# helper function (optional)
@inline function custom_viscosity(a::CustomRheology; P = 0.0, T = 273.0, depth = 0.0, kwargs...)
    η0, Ea, Va, T0, R, cutoff = a.args.η0,
        a.args.Ea, a.args.Va, a.args.T0, a.args.R,
        a.args.cutoff
    η = η0 * exp((Ea + P * Va) / (R * T) - Ea / (R * T0))
    correction = (depth ≤ 660.0e3) + (depth > 660.0e3) * 1.0e1
    return clamp(η * correction, cutoff...)
end

# constant parameters, these are typically wrapped into a struct (compulsory)
v_args = (; η0 = 5.0e20, Ea = 200.0e3, Va = 2.6e-6, T0 = 1.6e3, R = 8.3145, cutoff = (1.0e16, 1.0e25))

# create rheology struct
v1 = CustomRheology(custom_εII, custom_τII, v_args)

@testset "CustomRheology" begin
    args = (; depth = 1.0e3, P = 1.0e3 * 3300 * 9.81, T = 1.6e3)
    τII, εII = 1.0e3, 1.0e-15
    ε = compute_εII(v1, τII, args)
    τ = compute_τII(v1, εII, args)
    dεdτ = dεII_dτII(v1, τII, args)
    dτdε = dτII_dεII(v1, εII, args)
    @test (ε, τ, dεdτ, dτdε) == (
        9.936929394365625e-19,
        1.0063470920574456e6,
        9.936929394365625e-22,
        1.0063470920574455e21,
    )

    # CompositeRheologies
    v2 = ConstantElasticity()
    c = CompositeRheology(v1, v2)
    args = (; dt = Inf, depth = 1.0e3, P = 1.0e3 * 3300 * 9.81, T = 1.6e3)
    ε = compute_εII(c, τII, args)
    τ = compute_τII(c, εII, args)
    @test (ε, τ) == (9.936929394365625e-19, 1.0063470920574456e6)
end
