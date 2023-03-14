using Test
using GeoParams

x1 = MeltViscosity()



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

end