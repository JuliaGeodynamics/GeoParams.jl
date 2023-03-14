using Test
using GeoParams


@testset "MeltViscosity.jl" begin

    x1 = MeltViscosity()

    @test Value(x1.α) == 1.0

    ϕ = 0.05
    η = 1e18

    args = (; η=η, ϕ=ϕ)
    TauII = 1e6

    ε = compute_εII(x1 , TauII, args)

end