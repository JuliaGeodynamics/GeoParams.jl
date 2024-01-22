using GeoParams, Test

@testset "LinearSoftening" begin
    min_v, max_v = 15e0, 30e0
    lo, hi = 0.0, 1.0
    
    @test LinearSoftening(min_v, max_v, lo, hi) === LinearSoftening((min_v, max_v), (lo, hi))
    
    soft = LinearSoftening(min_v, max_v, lo, hi)
    
    @test soft(1) == min_v
    @test soft(0) == max_v
    @test soft(0.5) == 0.5 * (min_v + max_v)
    
end

min_v, max_v = 15e0, 30e0
lo, hi = 0.0, 1.0

@test LinearSoftening(min_v, max_v, lo, hi) === LinearSoftening((min_v, max_v), (lo, hi))

soft = LinearSoftening(min_v, max_v, lo, hi)

@test soft(1) == min_v
@test soft(0) == max_v
@test soft(0.5) == 0.5 * (min_v + max_v)
