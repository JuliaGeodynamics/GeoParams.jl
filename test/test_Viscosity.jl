using Test, GeoParams

@testset "Viscosity" begin
    εII = 1
    τII = 1
    η0 = 1
    dt = 1
    P, T, xx, yy, xy = (rand() for i in  1:5);
    args = (; P = P, T = T, dt = dt);

    # Physical properties using GeoParams ----------------
    # create rheology struct
    el       = SetConstantElasticity(; G=1, ν=0.5)            
    creep    = LinearViscous(; η = η0)       # Arrhenius-like (T-dependant) viscosity
    rheology = SetMaterialParams(;
        CompositeRheology = CompositeRheology((creep, el)),
    )

    @test η0 == compute_viscosity_εII(rheology, τII, args) == compute_viscosity_τII(rheology, τII, args)
    @test 1(1/η0 + 1/el.G.val/dt) == compute_elastoviscosity_εII(rheology, εII, args) == compute_elastoviscosity_τII(rheology, τII, args)

    # Slightly more complex example ----------------------
    # function to compute strain rate 
    @inline function custom_εII(a::CustomRheology, TauII; args...)
        η = custom_viscosity(a; args...)
        return (TauII / η) * 0.5
    end

    # function to compute deviatoric stress 
    @inline function custom_τII(a::CustomRheology, EpsII; args...)
        η = custom_viscosity(a; args...)
        return 2.0 * (η * EpsII)
    end

    # helper function (optional)
    @inline function custom_viscosity(a::CustomRheology; P=0.0, T=273.0, depth=0.0, kwargs...)
        η0, Ea, Va, T0, R, cutoff = a.args.η0,
        a.args.Ea, a.args.Va, a.args.T0, a.args.R,
        a.args.cutoff
        η = η0 * exp((Ea + P * Va) / (R * T) - Ea / (R * T0))
        correction = (depth ≤ 660e3) + (depth > 660e3) * 1e1
        return clamp(η * correction, cutoff...)
    end

    # constant parameters, these are typically wrapped into a struct 
    v_args = (; η0=5e20, Ea=200e3, Va=2.6e-6, T0=1.6e3, R=8.3145, cutoff=(1e16, 1e25))

    args = (; depth=1e3, P=1e3 * 3300 * 9.81, T=1.6e3)
    τII, εII = 1e3, 1e-15

    # create rheology struct
    v1 = CustomRheology(custom_εII, custom_τII, v_args)
    el = SetConstantElasticity(; G=70e9, ν=0.5)            
    rheology = SetMaterialParams(;
        CompositeRheology = CompositeRheology((v1, el)),
    )

    η = 5.0317354602872275e20
    @test η == compute_viscosity_τII(rheology, τII, args) == compute_viscosity_εII(rheology, εII, args)
    @test 1(1/η + 1/el.G.val/dt) == compute_elastoviscosity_εII(rheology, εII, args) == compute_elastoviscosity_τII(rheology, τII, args)
end