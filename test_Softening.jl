using GeoParams, Test

@testset "LinearSoftening" begin

    # Test NoSoftening
    soft = NoSoftening()
    x = rand()
    @test x === soft(x, rand())


    # Test LinearSoftening
    min_v, max_v = 15e0, 30e0
    lo, hi = 0.0, 1.0
    
    @test LinearSoftening(min_v, max_v, lo, hi) === LinearSoftening((min_v, max_v), (lo, hi))
    
    soft = LinearSoftening(min_v, max_v, lo, hi)
    
    @test soft(30, 1) == min_v
    @test soft(30, 0) == max_v
    @test soft(30, 0.5) == 0.5 * (min_v + max_v)

    # test plasticity structs
    d = DruckerPrager(; softening_C = soft)
    
end
