using Test, GeoParams

@testset "Gravity" begin

    gravity = 9.81

    gconst = ConstantGravity()
    @test compute_gravity(gconst) == 9.81
    @test isbits(gconst)
    @test param_info(gconst).Equation === L"g = 9.81 m s^{-2}"
    str = sprint(show, gconst)
    @test occursin("Gravitational acceleration:", str)

    gconst = ConstantGravity(; g = 1)
    @test compute_gravity(gconst) == 1
    @test isbits(gconst)
    @test sprint(show, gconst) == "Gravitational acceleration: g=1.0"

    d = DippingGravity(90, 0, gravity)
    gᵢ = compute_gravity(d)
    @test isbits(d)
    @test param_info(d).Equation === L"g = ($(s.gx, s.gy, s.gz)) m s^{-2}"
    @test sprint(show, d) == "Gravitational acceleration: g=(9.81, 0.0, 9.81)"
    @test gᵢ == (0.0e0, 0.0e0, gravity)

    α = 45
    d = DippingGravity(α, 0, gravity)
    gᵢ = compute_gravity(d)
    @test gᵢ == sind(90 - α) .* (gravity, 0, gravity)

    d = DippingGravity(α, α, gravity)
    gᵢ = compute_gravity(d)
    @test gᵢ == sind(90 - α) .* (gravity * cosd(45), gravity * cosd(45), gravity)

    r = SetMaterialParams(; Gravity = DippingGravity())
    @test compute_gravity(r) == (0.0e0, 0.0e0, gravity)
end
