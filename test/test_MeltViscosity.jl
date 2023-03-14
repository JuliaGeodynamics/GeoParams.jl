using Test
using GeoParams


@testset "Density.jl" begin

    x1 = MeltViscosity()

    @test Value(x.α) == 1.0

    ϕ = 0.05
    η = 1e18

    args = (; η=η, ϕ=ϕ)
    TauII = 1e6

    ε = compute_εII(x1 , TauII, args)

end