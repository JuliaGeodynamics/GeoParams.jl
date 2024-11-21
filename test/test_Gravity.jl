using Test, GeoParams

@testset "Gravity" begin

    g = 9.81

    gconst = ConstantGravity()
    @test gconst.g == convert(GeoUnit, 9.81m / s^2)

    gconst = ConstantGravity(; g = 1)
    @test gconst.g == convert(GeoUnit,1)

    d = DippingGravity(90, 0, g)
    gᵢ = compute_gravity2(d)
    @test gᵢ == (0e0, 0e0, g)

    d = DippingGravity(45, 0, g)
    gᵢ = compute_gravity2(d)
    @test gᵢ == sind(90 - α) .* (g, 0, g)

    d = DippingGravity(45, 45, g)
    gᵢ = compute_gravity2(d)
    @test gᵢ ==   sind(90 - α) .* (g *  cosd(45), g *  cosd(45), g)
end
