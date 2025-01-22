using Test, GeoParams, StaticArrays
import ForwardDiff as FD
@testset "Viscosity" begin
    εII = 1
    τII = 2
    η0 = 3
    dt = 4
    P, T, xx, yy, xy = (rand() for i in 1:5)
    args = (; P = P, T = T, dt = dt)

    # Physical properties using GeoParams ----------------
    # create rheology struct
    G = 1
    el = SetConstantElasticity(; G = G, ν = 0.5)
    creep = LinearViscous(; η = η0)
    rheology = SetMaterialParams(; CompositeRheology = CompositeRheology((creep, el)))

    @test compute_viscosity(el, args) == G * dt
    @test compute_viscosity(creep, args) == η0
    @test compute_viscosity(rheology, args) == 1 / (1 / η0 + 1 / G / dt)

    # Test differentiability
    @test FD.derivative(x -> compute_viscosity(el, (; P = P, T = x, dt = dt)), T) == 0.0
    @test FD.derivative(x -> compute_viscosity(creep, (; P = P, T = x, dt = dt)), T) == 0.0
    @test FD.derivative(x -> compute_viscosity(rheology, (; P = P, T = x, dt = dt)), T) == 0.0

    @test_throws "compute_viscosity only works for linear rheologies" compute_viscosity(DislocationCreep(), args)

    @test η0 ==
        compute_viscosity_εII(rheology, 0.0, args) ==
        compute_viscosity_τII(rheology, 0.0, args)
    @test 1 / (1 / η0 + 1 / G / dt) ==
        compute_elastoviscosity_εII(rheology, εII, args) ==
        compute_elastoviscosity_τII(rheology, τII, args)

    @test FD.derivative(x -> compute_viscosity_εII(rheology, εII, (; P = P, T = x, dt = dt)), T) == 0.0
    @test FD.derivative(x -> compute_viscosity_τII(rheology, τII, (; P = P, T = x, dt = dt)), T) == 0.0

    rheologies1 = (
        SetMaterialParams(;
            CompositeRheology = CompositeRheology((LinearViscous(; η = 1),))
        ),
        SetMaterialParams(;
            CompositeRheology = CompositeRheology((LinearViscous(; η = 2),))
        ),
    )

    @test compute_viscosity(rheologies1, 1, args) == 1.0
    @test compute_viscosity(rheologies1, 2, args) == 2.0

    phase_ratio = (0.5, 0.5)
    @test compute_viscosity(rheologies1, phase_ratio, args) == 1.5
    @test compute_viscosity(rheologies1, SA[0.5, 0.5], args) == 1.5

    @test FD.derivative(x -> compute_viscosity(rheologies1, 1, (; P = P, T = x, dt = dt)), T) == 0.0
    @test FD.derivative(x -> compute_viscosity(rheologies1, 2, (; P = P, T = x, dt = dt)), T) == 0.0
    @test FD.derivative(x -> compute_viscosity(rheologies1, phase_ratio, (; P = P, T = x, dt = dt)), T) == 0.0
    @test FD.derivative(x -> compute_viscosity(rheologies1, SA[0.5, 0.5], (; P = P, T = x, dt = dt)), T) == 0.0

    rheologies2 = (
        SetMaterialParams(;
            CompositeRheology = CompositeRheology((LinearViscous(; η = 1), el))
        ),
        SetMaterialParams(;
            CompositeRheology = CompositeRheology((LinearViscous(; η = 2), el))
        ),
    )

    @test compute_viscosity(rheologies2, 1, args) == 0.8
    @test compute_viscosity(rheologies2, 2, args) == 1.3333333333333333

    @test compute_viscosity(rheologies2, phase_ratio, args) == 1.0666666666666667
    @test compute_viscosity(rheologies2, SA[0.5, 0.5], args) == 1.0666666666666667

    @test compute_elasticviscosity(rheologies2, 1, args) == el.G * dt
    @test compute_elasticviscosity(rheologies2, 2, args) == el.G * dt
    @test compute_elasticviscosity(rheologies2, phase_ratio, args) == el.G * dt

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
            a::CustomRheology; P = 0.0, T = 273.0, depth = 0.0, kwargs...
        )
        (; η0, Ea, Va, T0, R, cutoff) = a.args
        η = η0 * exp((Ea + P * Va) / (R * T) - Ea / (R * T0))
        correction = (depth ≤ 660.0e3) + (depth > 660.0e3) * 1.0e1
        return clamp(η * correction, cutoff...)
    end

    # constant parameters, these are typically wrapped into a struct
    v_args = (; η0 = 5.0e20, Ea = 200.0e3, Va = 2.6e-6, T0 = 1.6e3, R = 8.3145, cutoff = (1.0e16, 1.0e25))

    dt = 100.0e3 * 3600 * 24 * 365
    args = (; depth = 1.0e3, P = 1.0e3 * 3300 * 9.81, T = 1.6e3, dt = dt)
    τII, εII = 1.0e3, 1.0e-15

    # create rheology struct
    v1 = CustomRheology(custom_εII, custom_τII, v_args)
    el = SetConstantElasticity(; G = 70.0e9, ν = 0.5)
    rheology = SetMaterialParams(; CompositeRheology = CompositeRheology((v1, el)))

    η = 5.0317354602872275e20
    @test η ==
        compute_viscosity_τII(rheology, τII, args) ==
        compute_viscosity_εII(rheology, εII, args)
    @test 1 / (1 / η + 1 / el.G / dt) ==
        compute_elastoviscosity_εII(rheology, εII, args) ==
        compute_elastoviscosity_τII(rheology, τII, args)

    # dη / dT
    @test FD.derivative(x -> compute_viscosity_τII(rheology, τII, (; depth = 1.0e3, P = 1.0e3 * 3300 * 9.81, T = x, dt = dt)), 1.6e3) == -4.729926879551493e18
    @test FD.derivative(x -> compute_viscosity_εII(rheology, εII, (; depth = 1.0e3, P = 1.0e3 * 3300 * 9.81, T = x, dt = dt)), 1.6e3) == -4.729926879551493e18
    @test FD.derivative(x -> compute_elastoviscosity_τII(rheology, τII, (; depth = 1.0e3, P = 1.0e3 * 3300 * 9.81, T = x, dt = dt)), 1.6e3) == -4.708437955242587e18
    @test FD.derivative(x -> compute_elastoviscosity_εII(rheology, εII, (; depth = 1.0e3, P = 1.0e3 * 3300 * 9.81, T = x, dt = dt)), 1.6e3) == -4.708437955242587e18
    # dη / dP
    @test FD.derivative(x -> compute_viscosity_τII(rheology, τII, (; depth = 1.0e3, P = x, T = 1.6e3, dt = dt)), args.P) == 9.834109234429906e10
    @test FD.derivative(x -> compute_viscosity_εII(rheology, εII, (; depth = 1.0e3, P = x, T = 1.6e3, dt = dt)), args.P) == 9.834109234429906e10
    @test FD.derivative(x -> compute_elastoviscosity_τII(rheology, τII, (; depth = 1.0e3, P = x, T = 1.6e3, dt = dt)), args.P) == 9.789431074626257e10
    @test FD.derivative(x -> compute_elastoviscosity_εII(rheology, εII, (; depth = 1.0e3, P = x, T = 1.6e3, dt = dt)), args.P) == 9.789431074626257e10

    T = LinRange(0, 1.6e3, 100)
    P = LinRange(0, 9.0e9, 100)
    dT = zeros(100, 100)
    dP = zeros(100, 100)

    for (j, T) in enumerate(T), (i, P) in enumerate(P)
        dT[i, j] =
            FD.derivative(
            x -> compute_viscosity_τII(rheology, τII, (; depth = 1.0e3, P = P, T = x, dt = dt)),
            T
        )
        dP[i, j] =
            FD.derivative(
            x -> compute_viscosity_τII(rheology, τII, (; depth = 1.0e3, P = x, T = T, dt = dt)),
            P
        )
    end
end

# using GLMakie

# fig=Figure()
# ax1 = Axis(fig[1, 1])
# ax2 = Axis(fig[2, 1])
# # h1 = heatmap!(ax1, @. log10(abs(dT)))
# # h2 = heatmap!(ax2, @. log10(abs(dP)))
# h1 = heatmap!(ax1, abs.(dT), colormap=:lipari)
# h2 = heatmap!(ax2, abs.(dP), colormap=:lipari)
# Colorbar(fig[1,2], h1)
# Colorbar(fig[2,2], h2)
# fig
