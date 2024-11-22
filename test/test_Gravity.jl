using Test, GeoParams

@testset "Gravity" begin

    gravity = 9.81

    gconst = ConstantGravity()
    @test compute_gravity(gconst) == 9.81

    gconst = ConstantGravity(; g = 1)
    @test compute_gravity(gconst) == 1

    d = DippingGravity(90, 0, gravity)
    gᵢ = compute_gravity(d)
    @test gᵢ == (0e0, 0e0, gravity)

    α = 45
    d = DippingGravity(α, 0, gravity)
    gᵢ = compute_gravity(d)
    @test gᵢ == sind(90 - α) .* (gravity, 0, gravity)

    d = DippingGravity(α, α, gravity)
    gᵢ = compute_gravity(d)
    @test gᵢ ==   sind(90 - α) .* (gravity *  cosd(45), gravity *  cosd(45), gravity)

    r = SetMaterialParams(; Gravity = DippingGravity(), )
    @test compute_gravity(r) == (0e0, 0e0, gravity)
end

