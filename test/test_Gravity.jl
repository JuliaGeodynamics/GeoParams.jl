using Test, GeoParams

@testset "Gravity" begin

    g = 9.81

    gconst = ConstantGravity()
    @test compute_gravity(gconst) == 9.81

    gconst = ConstantGravity(; g = 1)
    @test compute_gravity(gconst) == 1

    d = DippingGravity(90, 0, g)
    gᵢ = compute_gravity(d)
    @test gᵢ == (0e0, 0e0, g)

    α = 45
    d = DippingGravity(α, 0, g)
    gᵢ = compute_gravity(d)
    @test gᵢ == sind(90 - α) .* (g, 0, g)

    d = DippingGravity(α, α, g)
    gᵢ = compute_gravity(d)
    @test gᵢ ==   sind(90 - α) .* (g *  cosd(45), g *  cosd(45), g)
end
