using GeoParams, Test

@testset "Softening" begin

    # Test NoSoftening
    soft = NoSoftening()
    x = rand()
    @test x === soft(rand(), x)

    # Test LinearSoftening
    min_v, max_v = rand() * 15, (rand() + 1) * 15
    lo, hi = 0.0, 1.0

    @test LinearSoftening(min_v, max_v, lo, hi) === LinearSoftening((min_v, max_v), (lo, hi))

    soft_ϕ = LinearSoftening(min_v, max_v, lo, hi)

    @test soft_ϕ(1, max_v) == min_v
    @test soft_ϕ(0, max_v) == max_v
    @test soft_ϕ(0.5, max_v) ≈ 0.5 * (min_v + max_v)

    min_v, max_v = 20.0e0, 20.0e0
    soft_ϕ = LinearSoftening(min_v, max_v, lo, hi)

    @test soft_ϕ(1, max_v) == 20.0e0
    @test soft_ϕ(0, max_v) == 20.0e0
    @test soft_ϕ(0.5, max_v) == 20.0e0

    # test Drucker-Prager with softening
    min_v, max_v = 15.0e0, 30.0e0
    lo, hi = 0.0, 1.0
    soft_ϕ = LinearSoftening(min_v, max_v, lo, hi)
    softening_C = LinearSoftening((0.0e0, 10.0e6), (lo, hi))

    τII = 20.0e6
    P = 1.0e6
    args = (P = P, τII = τII)

    p = DruckerPrager()
    @test isbits(p)
    @test compute_yieldfunction(p, args) ≈ 1.0839745962155614e7
    args = (P = P, τII = τII, EII = 1.0e0)

    p1 = DruckerPrager(; softening_ϕ = soft_ϕ)
    p2 = DruckerPrager(; softening_C = LinearSoftening((0.0e0, 10.0e6), (lo, hi)))
    p3 = DruckerPrager(; softening_ϕ = soft_ϕ, softening_C = softening_C)
    p4 = DruckerPrager(; softening_ϕ = DecaySoftening())
    p5 = DruckerPrager(; softening_C = softening_C, softening_ϕ = DecaySoftening())

    @test compute_yieldfunction(p1, args) ≈ 1.0081922692006797e7
    @test compute_yieldfunction(p2, args) ≈ 1.95e7
    @test compute_yieldfunction(p3, args) ≈ 1.974118095489748e7
    @test compute_yieldfunction(p4, args) ≈ 9.977203951679487e6
    @test compute_yieldfunction(p5, args) ≈ 1.997376090963723e7

    # test regularized Drucker-Prager with softening
    p = DruckerPrager_regularised()
    @test isbits(p)
    @test compute_yieldfunction(p, args) ≈ 1.0839745962155614e7

    args = (P = P, τII = τII, EII = 1.0e0)

    p1 = DruckerPrager_regularised(; softening_ϕ = soft_ϕ)
    p2 = DruckerPrager_regularised(; softening_C = LinearSoftening((0.0e0, 10.0e6), (lo, hi)))
    p3 = DruckerPrager_regularised(; softening_ϕ = soft_ϕ, softening_C = LinearSoftening((0.0e0, 10.0e6), (lo, hi)))

    @test compute_yieldfunction(p1, args) ≈ 1.0081922692006797e7
    @test compute_yieldfunction(p2, args) ≈ 1.95e7
    @test compute_yieldfunction(p3, args) ≈ 1.974118095489748e7

    # non linear softening
    p4 = DruckerPrager_regularised(; softening_C = NonLinearSoftening())
    p5 = DruckerPrager_regularised(; softening_C = NonLinearSoftening(ξ₀ = 30, Δ = 10))

    @test compute_yieldfunction(p4, args) ≈ 1.95e7
    @test compute_yieldfunction(p5, (P = P, τII = τII, EII = 0.0e0)) ≈ 1.9499974039493073e7
    @test compute_yieldfunction(p5, (P = P, τII = τII, EII = 10.0e0)) ≈ 1.9499982679491926e7

    ### Test nondimensionalization
    CharDim   = GEO_units(;length = 100km, viscosity = 1e21Pa*s)

    Coh = 1.0e6Pa
    ls=LinearSoftening(Coh/2, Coh, 0.0NoUnits, 0.5NoUnits)

    ls_nd = nondimensionalize(ls, CharDim)
    @test ls_nd.min_value.val == 0.05 
    @test ls_nd.max_value.val == 0.1
    @test ls_nd.hi.val == 0.5
    @test ls_nd.lo.val == 0

    nls = NonLinearSoftening(ξ₀=30MPa, Δ=10MPa)
    nls_nd = nondimensionalize(nls, CharDim)

    @test nls_nd.ξ₀.val == 3
    @test nls_nd.Δ.val == 1

    ds    = DecaySoftening(εref = 1e-15/s)
    ds_nd = nondimensionalize(ds, CharDim)

    @test ds_nd.εref.val == 0.1

end
