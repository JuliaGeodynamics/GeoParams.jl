using Test, Statistics
using GeoParams

@testset "CreepLaw" begin

    #Make sure structs are isbits
    x = LinearViscous()
    @test isbits(x)
    @test isvolumetric(x) == false

    x = PowerlawViscous()
    @test isbits(x)

    x = ArrheniusType()
    @test isbits(x)

    x = LinearMeltViscosity()
    @test isbits(x)

    x = ViscosityPartialMelt_Costa_etal_2009()
    @test isbits(x)

    x = GiordanoMeltViscosity()
    @test isbits(x)
    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 1000km)

    # Define a linear viscous creep law ---------------------------------
    x1 = LinearViscous(; η = 1.0e18Pa * s)
    @test x1.η.val == 1.0e18

    x1_ND = LinearViscous(; η = 1.0e18Pa * s)
    @test isDimensional(x1_ND) == true
    x1_ND = nondimensionalize(x1_ND, CharUnits_GEO)                # check that we can nondimensionalize all entries within the struct
    @test isDimensional(x1_ND) == false
    @test x1_ND.η * 1.0 == 0.1
    x1_D = dimensionalize(x1_ND, CharUnits_GEO)                    # check that we can dimensionalize it again
    @test x1_D.η.val == 1.0e18

    # Given stress
    args = (;)
    ε_D = compute_εII(x1, 1.0e6Pa, args)
    @test ε_D == 5.0e-13 / s                      # dimensional input
    ε_ND = compute_εII(x1_ND, nondimensionalize(1.0e6Pa, CharUnits_GEO), args)
    @test ε_ND ≈ 0.5                            # non-dimensional
    @test ε_ND * CharUnits_GEO.strainrate ≈ ε_D   # check that non-dimensiional computations give the same results

    τ = [0.0; 0.0]
    compute_εII!(τ, x1_ND, [1.0e0; 2.0], args)
    @test τ == [5.0; 10.0]    # vector input

    # Given strainrate
    @test compute_τII(x1, 1.0e-13 / s, args) == 1.0e18 * 2 * 1.0e-13Pa       # dimensional input
    @test compute_τII(x1_ND, 1.0e0, args) == 0.2                   # non-dimensional
    ε = [0.0; 0.0]
    compute_τII!(ε, x1_ND, [1.0e0; 2.0], args)
    @test ε == [0.2; 0.4]     # vector input
    # -------------------------------------------------------------------

    # -------------------------------------------------------------------
    # Define powerlaw viscous rheology
    η0 = 10
    n = 2.3
    ε0 = 1
    x2 = PowerlawViscous(; η0 = 10, n = 2.3, ε0 = 1)

    τII = 1.0e0
    εII = 1.0e0

    (τII / η0)^(1 / n)
    @test compute_τII(x2, τII) == η0 * εII^n
    @test compute_εII(x2, εII) ≈ (τII / η0)^(1 / n)

    @test compute_viscosity_εII(x2, τII, (;)) == 5.0e0
    @test compute_viscosity_τII(x2, εII, (;)) == 1.3606693841876538

    x2 = PowerlawViscous()
    x2 = nondimensionalize(x2, CharUnits_GEO)
    @test NumValue(x2.ε0) == 0.001 # powerlaw

    # Given strainrate
    compute_τII(x1, 1.0e-13 / s, args)
    @test compute_τII(x1, 1.0e-13 / s, args) == 1.0e18 * 2 * 1.0e-13Pa       # dimensional input
    @test compute_τII(x1_ND, 1.0e0, args) == 0.2                   # non-dimensional
    @test [compute_τII(x1_ND, EII, args) for EII in (1.0e0, 2.0)] == [0.2; 0.4]     # vector input
    # -------------------------------------------------------------------

    # -------------------------------------------------------------------
    # Define powerlaw viscous rheology
    x2 = PowerlawViscous()
    x2 = nondimensionalize(x2, CharUnits_GEO)
    @test x2.ε0.val == 0.001                                 # powerlaw
    # -------------------------------------------------------------------

    # ArrheniusType rheology --------------------------------------------
    args = (;)
    x3 = ArrheniusType()
    # For a given temperature
    T = 1.0
    # and stress
    τ1 = 1.0
    #eps     =   compute_εII(x3,τ1;T,args)
    #@test   eps ≈ 0.5

    args = (T = T,)
    eps = compute_εII(x3, τ1, args)
    @test eps ≈ 0.5

    # and strain rate
    ε1 = 0.5
    tau = compute_τII(x3, ε1, args)
    @test tau ≈ 1.0
    # using a vector input ------------------
    T2 = [1.0; 0.9; 0.8]
    # and stress
    ε21 = [0.0; 0.0; 0.0]
    τ21 = [1.0; 0.8; 0.9]
    args = (T = T2,)
    compute_εII!(ε21, x3, τ21, args)
    # and strain rate
    τ22 = [0.0; 0.0; 0.0]
    ε22 = [0.5; 1.0; 0.2]
    compute_τII!(τ22, x3, ε22, args)
    # derivatives
    sol1 = dεII_dτII(x3, τ22)
    @test sol1 ≈ 0.5
    sol2 = dτII_dεII(x3, ε21)
    @test sol2 ≈ 2.0

    # Test melt viscosity ---------------------------------
    x1 = LinearMeltViscosity()
    @test  x1.B.val == 13374.0


    x1_D = LinearMeltViscosity(A = -8.159, B = 2.405e+4K, T0 = -430.9606K)   # Rhyolite
    @test isDimensional(x1_D) == true
    x1_ND = nondimensionalize(x1_D, CharUnits_GEO)                 # check that we can nondimensionalize all entries within the struct
    @test isDimensional(x1_ND) == false
    @test  NumValue(x1_ND.B) ≈ 18.890154341593682
    x1_D1 = dimensionalize(x1_ND, CharUnits_GEO)                    # check that we can dimensionalize it again
    @test  Value(x1_D.B) ≈ Value(x1_D1.B)

    # Given stress
    args_D = (; T = 700K)
    args_ND = (; T = nondimensionalize(700K, CharUnits_GEO))

    ε_D = compute_εII(x1_D, 1.0e6Pa, args_D)
    @test ε_D ≈ 3.916168662376774e-8 / s                      # dimensional input
    ε_ND = compute_εII(x1_ND, nondimensionalize(1.0e6Pa, CharUnits_GEO), args_ND)
    @test ε_ND ≈ 39161.68662376807                              # non-dimensional
    @test ε_ND * CharUnits_GEO.strainrate ≈ ε_D   # check that non-dimensiional computations give the same results

    ε = [0.0; 0.0]
    compute_εII!(ε, x1_ND, [1.0e6; 2.0e6], args_ND)
    @test ε ≈ [3.9161686623768066e11; 7.832337324753613e11]    # vector input

    # Given strainrate
    @test compute_τII(x1_D, 1.0e-13 / s, args_D) ≈ 2.553516169022932Pa      # dimensional input
    @test compute_τII(x1_ND, 1.0e-13, args_ND) ≈ 2.553516169022911e-19                  # non-dimensional

    ε = [0.0; 0.0]
    compute_τII!(ε, x1_ND, [1.0e0; 2.0], args_ND)
    @test ε ≈ [2.553516169022911e-6; 5.107032338045822e-6]     # vector input

    # With vector as input
    T = Vector(800.0:1400) * K
    η_basalt = zeros(size(T)) * Pas
    η_rhyolite = zeros(size(T)) * Pas
    x_basalt = LinearMeltViscosity()
    x_rhyolite = LinearMeltViscosity(A = -8.159, B = 2.405e+4K, T0 = -430.9606K)   # Rhyolite

    for i in eachindex(T)
        args_D = (; T = T[i])
        η_basalt[i] = compute_viscosity_τII(x_basalt, 1.0e6Pa, (; T = T[i]))
        η_rhyolite[i] = compute_viscosity_τII(x_rhyolite, 1.0e6Pa, (; T = T[i]))
    end

    args_D = (; T = 1000K)
    η_basalt1 = compute_viscosity_τII(x_basalt, 1.0e6Pa, args_D)
    @test η_basalt1 ≈ 5.2471745814062805e9Pa * s
    η_rhyolite1 = compute_viscosity_τII(x_rhyolite, 1.0e6Pa, args_D)
    @test η_rhyolite1 ≈ 4.445205243727607e8Pa * s
    # -------------------------------------------------------------------

    # Test effective viscosity of partially molten rocks -----------------
    x1_D = ViscosityPartialMelt_Costa_etal_2009(η = LinearMeltViscosity(A = -8.159, B = 2.405e+4K, T0 = -430.9606K))   # Rhyolite
    @test isDimensional(x1_D) == true
    x1_ND = nondimensionalize(x1_D, CharUnits_GEO)                 # check that we can nondimensionalize all entries within the struct
    @test isDimensional(x1_ND) == false
    @test  NumValue(x1_ND.η.B) ≈ 18.890154341593682
    x1_D1 = dimensionalize(x1_ND, CharUnits_GEO)                    # check that we can dimensionalize it again
    @test  Value(x1_D.η.B) ≈ Value(x1_D1.η.B)

    # Given stress
    args_D = (; T = 700K, ϕ = 0.5)
    args_ND = (; T = nondimensionalize(700K, CharUnits_GEO), ϕ = 0.5)

    ε_D = compute_εII(x1_D, 1.0e6Pa, args_D)
    @test ε_D ≈ 1.5692922423522311e-10 / s                      # dimensional input
    ε_ND = compute_εII(x1_ND, nondimensionalize(1.0e6Pa, CharUnits_GEO), args_ND)
    @test ε_ND ≈ 156.92922423522444                              # non-dimensional
    @test ε_ND * CharUnits_GEO.strainrate ≈ ε_D   # check that non-dimensiional computations give the same results

    ε = [0.0; 0.0]
    τ = [1.0e6; 2.0e6]
    args_ND1 = (; T = 1000 * ones(size(τ)), ϕ = 0.5 * ones(size(τ)))
    compute_εII!(ε, x1_D, τ, args_ND1)
    @test ε ≈ [0.0011248074556411618; 0.0022496149112823235]    # vector input

    # Given strainrate
    @test compute_τII(x1_D, 1.0e-13 / s, args_D) ≈ 643.9044043803415Pa      # dimensional input
    @test compute_τII(x1_ND, 1.0e-13, args_ND) ≈ 5.64401178083053e-17                  # non-dimensional

    ε = [0.0; 0.0]
    compute_τII!(ε, x1_ND, [1.0e0; 2.0], args_ND1)
    @test ε ≈ [3.65254588484434e-25, 7.305079089024537e-25]    # vector input

    x1_D = ViscosityPartialMelt_Costa_etal_2009()
    args_D = (; T = 1000K, ϕ = 0.5)
    ε_D = compute_εII(x1_D, 1.0e6Pa, args_D)
    @test ustrip(ε_D) ≈ 3.4537433628446676e-6
    η = compute_viscosity_τII(x1_D, 1.0e6Pa, args_D)
    @test  ustrip(η) ≈ 1.447704555523709e11

    @test dεII_dτII(x1_ND, 1.0e6, args_ND) ≈ 31884.63277076202
    @test dτII_dεII(x1_ND, 1.0e-15, args_ND) ≈ 0.0005644011780830529
    # -----

    # ----
    # Create plot
    T = Vector(800.0:1400) * K
    η_basalt = zeros(size(T)) * Pas
    η_basalt1 = zeros(size(T)) * Pas
    η_rhyolite = zeros(size(T)) * Pas
    ϕ_basalt = zeros(size(T))
    ϕ_rhyolite = zeros(size(T))
    x_basalt = ViscosityPartialMelt_Costa_etal_2009(η = LinearMeltViscosity())
    x_rhyolite = ViscosityPartialMelt_Costa_etal_2009(η = LinearMeltViscosity(A = -8.159, B = 2.405e+4K, T0 = -430.9606K))   # Rhyolite
    mf_basalt = MeltingParam_Smooth3rdOrder()
    mf_rhyolite = MeltingParam_Smooth3rdOrder(a = 3043.0, b = -10552.0, c = 12204.9, d = -4709.0)

    εII = 1.0e-5 / s
    for i in eachindex(T)
        args_D = (; T = T[i])
        ϕ_basalt[i] = compute_meltfraction(mf_basalt, (; T = T[i] / K))
        η_basalt[i] = compute_viscosity_εII(x_basalt, εII, (; ϕ = ϕ_basalt[i], T = T[i]))
        η_basalt1[i] = compute_viscosity_εII(x_basalt, 1.0e-3 / s, (; ϕ = ϕ_basalt[i], T = T[i]))
        ϕ_rhyolite[i] = compute_meltfraction(mf_rhyolite, (; T = T[i] / K))
        η_rhyolite[i] = compute_viscosity_εII(x_rhyolite, εII, (; ϕ = ϕ_rhyolite[i], T = T[i]))

    end
    @test mean(η_rhyolite) ≈ 1.2614388007523279e19Pas
    @test mean(η_basalt) ≈ 5.822731206904198e24Pas
    @test mean(η_basalt1) ≈ 8.638313729394237e21Pas

    #=
    using Plots
    plt1 = plot(ustrip.(T) .- 273.15, log10.(ustrip.(η_basalt)), xlabel="T [C]", ylabel="log10(η [Pas])", label="basalt, εII=1e-5/s")
    plot!(plt1, ustrip.(T) .- 273.15, log10.(ustrip.(η_basalt1)), label="basalt, εII=1e-3/s")
    plot!(plt1, ustrip.(T) .- 273.15, log10.(ustrip.(η_rhyolite)), label="rhyolite, εII=1e-5/s")

    plt2 = plot(ustrip.(T) .- 273.15, ϕ_basalt, xlabel="T [C]", ylabel="Melt fraction ϕ []", label="basalt")
    plot!(plt2, ustrip.(T) .- 273.15, ϕ_rhyolite, label="rhyolite")

    plot(plt1,plt2, layout=(2,1))
    =#

    # ----

    # Test melt viscosity from Giordano et al. (2008) ---------------------
    x1 = GiordanoMeltViscosity()
    @test  x1.AT.val == -4.55
    @test  x1.oxd_wt == (50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.14, 1.0)

    x1_D = GiordanoMeltViscosity(oxd_wt = (50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.14, 1.0))   # Rhyolite
    @test isDimensional(x1_D) == true
    x1_ND = nondimensionalize(x1_D, CharUnits_GEO)                 # check that we can nondimensionalize all entries within the struct
    @test isDimensional(x1_ND) == false
    @test  NumValue(x1_ND.AT) == -4.55
    x1_D1 = dimensionalize(x1_ND, CharUnits_GEO)                    # check that we can dimensionalize it again
    @test  Value(x1_D.AT) ≈ Value(x1_D1.AT)

    # given Temp
    args_D = (; T = 1273K)
    args_ND = (; T = nondimensionalize(1273K, CharUnits_GEO))
    Tau_II = 1.0e6Pa
    ε_D = compute_εII(x1_D, Tau_II, args_D)
    @test ε_D ≈ 1630.1958046211573 / s                      # dimensional input
    ε_ND = compute_εII(x1_ND, nondimensionalize(Tau_II, CharUnits_GEO), args_ND)
    @test ε_ND ≈ 1.630195804621171e15
    @test ε_ND * CharUnits_GEO.strainrate ≈ ε_D                          # non-dimensional

    ε = [0.0; 0.0]
    τ = [1.0e6; 2.0e6]
    args_ND1 = (; T = 1000 * ones(size(τ)), ϕ = 0.5 * ones(size(τ)))
    compute_εII!(ε, x1_D, τ, args_ND1)
    @test ε ≈ [0.28267381865073904, 0.5653476373014781] rtol = 1.0e-5

    # Given strain rate
    @test compute_τII(x1_D, 1.0e-13 / s, args_D) ≈ 6.1342324472022e-11Pa rtol = 1.0e-5     # dimensional input
    @test compute_τII(x1_ND, 1.0e-13, args_ND) ≈ 6.134232447202149e-30 rtol = 1.0e-5                 # non-dimensional


    ε = [0.0; 0.0]
    compute_τII!(ε, x1_ND, [1.0e0; 2.0], args_ND1)
    @test ε ≈ [5.6932943851735906e-24, 1.1386588770347181e-23]    # vector input

    ## intermediate composition
    oxd_int = (62.4, 0.55, 20.01, 0.03, 3.22, 9.08, 3.52, 0.93, 2.0)
    x1_D_int = GiordanoMeltViscosity(oxd_wt = oxd_int)   # Rhyolite

    @test dεII_dτII(x1_ND, nondimensionalize(Tau_II, CharUnits_GEO), args_ND) ≈ 1.6301958046211708e16 rtol = 1.0e-5
    @test dεII_dτII(x1_D_int, Tau_II, args_D) ≈ 0.00010359976542338855 * inv(Pas) rtol = 1.0e-5
    @test dτII_dεII(x1_ND, ε, args_ND) ≈ 6.134232447202149e-17 rtol = 1.0e-5
    @test dτII_dεII(x1_D, 1.0e-15s^-1, args_D) ≈ 613.4232447202199Pas
    η_ref = dτII_dεII(x1_D_int, 1.0e-15s^-1, args_D) / 2
    @test log10(ustrip.(η_ref)) ≈ (3.68) rtol = 1.0e-2

    # Create plot
    T = Vector(700.0:1673) * K
    η_MORB = zeros(size(T)) * Pas
    η_intfels = zeros(size(T)) * Pas
    η_rhyolite = zeros(size(T)) * Pas
    η_MORB_scaled = zeros(size(T)) * Pas
    η_intfels_scaled = zeros(size(T)) * Pas
    η_rhyolite_scaled = zeros(size(T)) * Pas
    x_MORB = GiordanoMeltViscosity(oxd_wt = (50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.14, 1.0))   # Rhyolite
    x_int = GiordanoMeltViscosity(oxd_wt = (62.4, 0.55, 20.01, 0.03, 3.22, 9.08, 3.52, 0.93, 2.0))   # Rhyolite
    x_rhyolite = GiordanoMeltViscosity(oxd_wt = (76.38, 0.06, 11.59, 1.03, 0.36, 3.25, 2.44, 4.66, 3.0))   # Rhyolite
    x_MORB_scaled = GiordanoMeltViscosity(oxd_wt = (50.42, 1.53, 15.13, 9.81, 7.76, 11.35, 2.83, 0.14, 1.0), η0 = 1.0e15Pas)   # Rhyolite
    x_int_scaled = GiordanoMeltViscosity(oxd_wt = (62.4, 0.55, 20.01, 0.03, 3.22, 9.08, 3.52, 0.93, 2.0), η0 = 1.0e15Pas)   # Rhyolite
    x_rhyolite_scaled = GiordanoMeltViscosity(oxd_wt = (76.38, 0.06, 11.59, 1.03, 0.36, 3.25, 2.44, 4.66, 3.0), η0 = 1.0e15Pas)   # Rhyolite

    εII = 1.0e5 / s
    for i in eachindex(T)
        args_D = (; T = T[i])
        η_MORB[i] = compute_viscosity_εII(x_MORB, εII, (; T = T[i]))
        η_intfels[i] = compute_viscosity_εII(x_int, εII, (; T = T[i]))
        η_rhyolite[i] = compute_viscosity_εII(x_rhyolite, εII, (; T = T[i]))
        η_MORB_scaled[i] = compute_viscosity_εII(x_MORB_scaled, εII, (; T = T[i]))
        η_intfels_scaled[i] = compute_viscosity_εII(x_int_scaled, εII, (; T = T[i]))
        η_rhyolite_scaled[i] = compute_viscosity_εII(x_rhyolite_scaled, εII, (; T = T[i]))
    end
    T_plots = 10000 ./ ustrip.(T)

    # fig = Figure(size = (800, 600))
    # ax = Axis(fig[1, 1], xlabel = "T [C]", ylabel = "log10(η [Pas])")
    # lines!(ax, ustrip.(T) .- 273.15, log10.(ustrip.(η_MORB)), label = "MORB, εII=1e-5/s")
    # lines!(ax, ustrip.(T) .- 273.15, log10.(ustrip.(η_intfels)), label = "Intermediate, εII=1e-5/s")
    # lines!(ax, ustrip.(T) .- 273.15, log10.(ustrip.(η_rhyolite)), label = "Rhyolite, εII=1e-5/s")


    # Legend(fig[1, 2], ax)
    # fig

    # fig = Figure(size = (800, 600))
    # ax = Axis(fig[1, 1], xlabel = "T [C]", ylabel = "log10(η [Pas])")
    # lines!(ax, ustrip.(T) .- 273.15, log10.(ustrip.(η_MORB_scaled)), label = "MORB_scaled, εII=1e-5/s, η0=1e15")
    # lines!(ax, ustrip.(T) .- 273.15, log10.(ustrip.(η_intfels_scaled)), label = "Intermediate_scaled, εII=1e-5/s, η0=1e15")
    # lines!(ax, ustrip.(T) .- 273.15, log10.(ustrip.(η_rhyolite_scaled)), label = "Rhyolite_scaled, εII=1e-5/s, η0=1e15")
    # Legend(fig[1, 2], ax)
    # fig
    # fig2 = Figure(size = (800, 600))
    # ax = Axis(fig2[1, 1], xlabel = "10000/T(K)", ylabel = "log10(η [Pas])")
    # lines!(ax, T_plots, log10.(ustrip.(η_MORB)), label = "MORB, εII=1e-5/s")
    # lines!(ax, T_plots, log10.(ustrip.(η_intfels)), label = "Intermediate, εII=1e-5/s")
    # lines!(ax, T_plots, log10.(ustrip.(η_rhyolite)), label = "Rhyolite, εII=1e-5/s")
    # # Legend(fig2[1, 2], ax, merge = true, position=:rt)
    # Legend(fig2[1, 2], ax, position=:rt)
    # fig2

end
