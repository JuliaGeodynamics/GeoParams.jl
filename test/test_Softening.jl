using GeoParams, Test

@testset "LinearSoftening" begin

    # Test NoSoftening
    soft = NoSoftening()
    x = rand()
    @test x === soft(rand(), x)

    # Test LinearSoftening
    min_v, max_v = rand()*15, (rand()+1)*15
    lo, hi = 0.0, 1.0
    
    @test LinearSoftening(min_v, max_v, lo, hi) === LinearSoftening((min_v, max_v), (lo, hi))
    
    soft_ϕ = LinearSoftening(min_v, max_v, lo, hi)
    
    @test soft_ϕ(1  , max_v) == min_v
    @test soft_ϕ(0  , max_v) == max_v
    @test soft_ϕ(0.5, max_v) ≈ 0.5 * (min_v + max_v)

    min_v, max_v = 20e0, 20e0
    soft_ϕ = LinearSoftening(min_v, max_v, lo, hi)
    
    @test soft_ϕ(1  , max_v) == 20e0
    @test soft_ϕ(0  , max_v) == 20e0
    @test soft_ϕ(0.5, max_v) == 20e0

    # test Drucker-Prager with softening
    min_v, max_v = 15e0, 30e0
    lo, hi = 0.0, 1.0
    soft_ϕ = LinearSoftening(min_v, max_v, lo, hi)
  
    τII = 20e6
    P = 1e6
    args = (P=P, τII=τII)

    p = DruckerPrager()
    @test compute_yieldfunction(p, args) ≈ 1.0839745962155614e7
    args = (P=P, τII=τII, EII=1e0)

    p1 = DruckerPrager(; softening_ϕ = soft_ϕ)
    p2 = DruckerPrager(; softening_C = LinearSoftening((0e0, 10e6), (lo, hi)))
    p3 = DruckerPrager(; softening_ϕ = soft_ϕ, softening_C = LinearSoftening((0e0, 10e6), (lo, hi)))

    @test compute_yieldfunction(p1, args) ≈ 1.0081922692006797e7
    @test compute_yieldfunction(p2, args) ≈ 1.95e7
    @test compute_yieldfunction(p3, args) ≈ 1.974118095489748e7

    # test regularized Drucker-Prager with softening
    p = DruckerPrager_regularised()
    @test compute_yieldfunction(p, args) ≈ 1.0839745962155614e7

    args = (P=P, τII=τII, EII=1e0)

    p1 = DruckerPrager(; softening_ϕ = soft_ϕ)
    p2 = DruckerPrager(; softening_C = LinearSoftening((0e0, 10e6), (lo, hi)))
    p3 = DruckerPrager(; softening_ϕ = soft_ϕ, softening_C = LinearSoftening((0e0, 10e6), (lo, hi)))

    @test compute_yieldfunction(p1, args) ≈ 1.0081922692006797e7
    @test compute_yieldfunction(p2, args) ≈ 1.95e7
    @test compute_yieldfunction(p3, args) ≈ 1.974118095489748e7
    
    # non linear softening 
    p4 = DruckerPrager(; softening_C = NonLinearSoftening())
    p5 = DruckerPrager(; softening_C = NonLinearSoftening(ξ₀ = 30, Δ = 10))

    @test compute_yieldfunction(p4, args) ≈ 1.95e7
    @test compute_yieldfunction(p5, (P=P, τII=τII, EII=0e0)) ≈ 1.9499974039493073e7
    @test compute_yieldfunction(p5, (P=P, τII=τII, EII=10e0)) ≈ 1.9499982679491926e7
end