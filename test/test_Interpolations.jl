using Test
using GeoParams
import GeoParams.Interpolations: LinearInterpolator, interpolate, interpolate_field, lerp, get_corners

@testset "Interpolations.jl" begin
    # 2x2 unit grid:  f(0,0)=1, f(1,0)=2, f(0,1)=3, f(1,1)=4
    data = [1.0 2.0; 3.0 4.0]
    itp = interpolate(0.0, 1.0, 2, 1.0, 0.0, 1.0, 2, 1.0, data)
    @test itp isa LinearInterpolator

    # interior bilinear value
    @test itp(0.5, 0.5) ≈ 2.5

    # corners
    @test itp(0.0, 0.0) ≈ 1.0
    @test itp(1.0, 1.0) ≈ 4.0

    # Flat extrapolation: clamp below range
    @test itp(-5.0, -5.0) ≈ 1.0

    # Flat extrapolation: clamp above range -> exercises the `x_clamped ≥ Tmax`
    # and `y_clamped ≥ Pmax` weight branches (t = 1, s = 1)
    @test itp(2.0, 2.0) ≈ 4.0
    @test itp(2.0, 0.0) ≈ 3.0   # x ≥ Tmax only
    @test itp(0.0, 2.0) ≈ 2.0   # y ≥ Pmax only

    # vectorized evaluation
    out = itp([0.0, 1.0], [0.0, 1.0])
    @test out ≈ [1.0, 4.0]

    # helper primitives
    @test lerp(0.0, 10.0, 0.25) ≈ 2.5
    @test get_corners(data, 1, 1) == (1.0, 3.0, 2.0, 4.0)  # column-major corners

    # interpolate_field: same math via the standalone entry point
    v = interpolate_field(0.0, 1.0, 2, 1.0, 0.0, 1.0, 2, 1.0, data, 0.5, 0.5)
    @test v ≈ 2.5
    # edge clamp through interpolate_field too
    @test interpolate_field(0.0, 1.0, 2, 1.0, 0.0, 1.0, 2, 1.0, data, 2.0, 2.0) ≈ 4.0
end
