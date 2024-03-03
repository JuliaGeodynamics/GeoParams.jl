using Test, GeoParams, StaticArrays

@testset "Viscosity" begin
    εII = 1
    τII = 2
    η0  = 3
    dt  = 4
    P, T, xx, yy, xy = (rand() for i in 1:5)
    args = (; P=P, T=T, dt=dt)

    # Physical properties using GeoParams ----------------
    # create rheology struct
    G        = 1
    el       = SetConstantElasticity(; G=G, ν=0.5)
    creep    = LinearViscous(; η=η0)
    rheology = SetMaterialParams(; CompositeRheology=CompositeRheology((creep, el)))

    @test compute_viscosity(el, args) == G * dt
    @test compute_viscosity(creep, args) == η0
    @test compute_viscosity(rheology, args) == 1/(1 / η0 + 1 / G / dt)
    
    @test η0 ==
        compute_viscosity_εII(rheology, 0.0, args) ==
        compute_viscosity_τII(rheology, 0.0, args)
    @test 1/(1 / η0 + 1 / el.G.val / dt) ==
        compute_elastoviscosity_εII(rheology, εII, args) ==
        compute_elastoviscosity_τII(rheology, τII, args)

    rheologies1 = (
        SetMaterialParams(; 
            CompositeRheology=CompositeRheology((LinearViscous(; η=1), ))
        ),
        SetMaterialParams(; 
            CompositeRheology=CompositeRheology((LinearViscous(; η=2), ))
        )
    )

    @test compute_viscosity(rheologies1, 1, args) == 1.0
    @test compute_viscosity(rheologies1, 2, args) == 2.0
    
    phase_ratio = (0.5, 0.5)
    @test compute_viscosity(rheologies1, phase_ratio, args) == 1.5
    @test compute_viscosity(rheologies1, SA[0.5, 0.5], args) == 1.5

    rheologies2 = (
        SetMaterialParams(; 
            CompositeRheology=CompositeRheology((LinearViscous(; η=1), el))
        ),
        SetMaterialParams(; 
            CompositeRheology=CompositeRheology((LinearViscous(; η=2), el))
        )
    )

    @test compute_viscosity(rheologies2, 1, args) == 0.8
    @test compute_viscosity(rheologies2, 2, args) == 1.3333333333333333
    
    @test compute_viscosity(rheologies2, phase_ratio, args)  == 1.0666666666666667
    @test compute_viscosity(rheologies2, SA[0.5, 0.5], args) == 1.0666666666666667

    @test compute_elasticviscosity(rheologies2, 1, args) == el.G.val * dt
    @test compute_elasticviscosity(rheologies2, 2, args) == el.G.val * dt
    @test compute_elasticviscosity(rheologies2, phase_ratio, args) == el.G.val * dt
    
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
    @inline function custom_viscosity(
        a::CustomRheology; P=0.0, T=273.0, depth=0.0, kwargs...
    )
        (; η0, Ea, Va, T0, R, cutoff) = a.args
        η = η0 * exp((Ea + P * Va) / (R * T) - Ea / (R * T0))
        correction = (depth ≤ 660e3) + (depth > 660e3) * 1e1
        return clamp(η * correction, cutoff...)
    end

    # constant parameters, these are typically wrapped into a struct 
    v_args = (; η0=5e20, Ea=200e3, Va=2.6e-6, T0=1.6e3, R=8.3145, cutoff=(1e16, 1e25))

    dt = 100e3 * 3600 * 24 * 365
    args = (; depth=1e3, P=1e3 * 3300 * 9.81, T=1.6e3, dt = dt)
    τII, εII = 1e3, 1e-15

    # create rheology struct
    v1 = CustomRheology(custom_εII, custom_τII, v_args)
    el = SetConstantElasticity(; G=70e9, ν=0.5)
    rheology = SetMaterialParams(; CompositeRheology=CompositeRheology((v1, el)))

    η = 5.0317354602872275e20
    @test η ==
        compute_viscosity_τII(rheology, τII, args) ==
        compute_viscosity_εII(rheology, εII, args)
    @test 1/(1 / η + 1 / el.G.val / dt) ==
        compute_elastoviscosity_εII(rheology, εII, args) ==
        compute_elastoviscosity_τII(rheology, τII, args)
end
