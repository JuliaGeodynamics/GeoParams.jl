using Test
using GeoParams


@testset "MeltViscosity.jl" begin

    # test structure
    x1 = MeltViscosity(α=0.8)
    @test Value(x1.α) == 0.8

    x1 = MeltViscosity()
    @test Value(x1.α) == 1.0

    ϕ = 0.05
    η = 1e18

    args = (; η=η, ϕ=ϕ)
    TauII = 1e8

    ε = compute_εII(x1, TauII, args)
    @test ε ≈ 5.5325795313118545e-11

    # calculate effective viscosity
    η_melt = computeViscosity_εII(x1, ε, args)
    @test η_melt ≈ 9.037375733511468e17

    # check that value for TauII is the same that previously assumed
    TauII2 = compute_τII(x1, ε, args)
    @test TauII2 == TauII

    # test for arrays
    τII_array = ones(10) * 1e8
    ϕ_array = ones(size(τII_array))* 0.05
    η_array = ones(size(τII_array))* 1e18
    ε_array = zeros(size(τII_array))

    args_array = (; η=η_array, ϕ=ϕ_array)

    compute_εII!(ε_array, x1, τII_array, args_array)
    @test ε_array[1] ≈ ε

    # compute when args are scalars
    compute_εII!(ε_array, x1, τII_array, args)
    @test ε_array[1] ≈ ε


end