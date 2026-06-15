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

    # test DPCap with dilation-angle softening
    softening_Ψ = LinearSoftening(0.0e0, 20.0e0, lo, hi)
    p_cap = DruckerPragerCap(; Ψ = 20.0, softening_Ψ = softening_Ψ)

    F_no_soft = compute_yieldfunction(p_cap, (P = P, τII = τII, EII = 0.0))
    F_soft = compute_yieldfunction(p_cap, (P = P, τII = τII, EII = 1.0))
    Q_no_soft = compute_flowpotential(p_cap, (P = P, τII = τII, EII = 0.0))
    Q_soft = compute_flowpotential(p_cap, (P = P, τII = τII, EII = 1.0))

    @test F_no_soft ≈ F_soft
    @test !isapprox(Q_no_soft, Q_soft)

    ### Test nondimensionalization
    CharDim = GEO_units(; length = 100km, viscosity = 1.0e21Pa * s)

    Coh = 1.0e6Pa
    ls = LinearSoftening(Coh / 2, Coh, 0.0NoUnits, 0.5NoUnits)

    ls_nd = nondimensionalize(ls, CharDim)
    @test ls_nd.min_value.val == 0.05
    @test ls_nd.max_value.val == 0.1
    @test ls_nd.hi.val == 0.5
    @test ls_nd.lo.val == 0

    nls = NonLinearSoftening(ξ₀ = 30MPa, Δ = 10MPa)
    nls_nd = nondimensionalize(nls, CharDim)

    @test nls_nd.ξ₀.val == 3
    @test nls_nd.Δ.val == 1

    ds = DecaySoftening(εref = 1.0e-15 / s)
    ds_nd = nondimensionalize(ds, CharDim)

    @test ds_nd.εref.val == 0.1

    # the stored slope must carry units and be rescaled along with min/max values
    @test isdimensional(ls.slope)
    @test ls_nd.slope.val ≈ 0.1

    # evaluation must be consistent across unit systems (regression for unscaled slope)
    char_stress = 1.0e7 # 10 MPa
    @test ls_nd(0.25, ls_nd.max_value.val) ≈ ls(0.25, ustrip(Coh)) / char_stress

    # evaluation for softening ranges other than (lo, hi) = (0, 1)
    @test ls(0.0, 1.0e6) ≈ 1.0e6
    @test ls(0.25, 1.0e6) ≈ 0.75e6
    @test ls(0.5, 1.0e6) ≈ 0.5e6

    # softening nested in a plastic element must survive nondimensionalization
    p_lin = DruckerPrager_regularised(; softening_C = LinearSoftening(5.0e6Pa, 1.0e7Pa, 0.0, 1.0), C = 1.0e7Pa, ϕ = 30.0, η_vp = 1.0e20Pa * s)
    p_lin_nd = nondimensionalize(p_lin, CharDim)
    for EII in (0.0, 0.5, 1.0)
        F_dim = compute_yieldfunction(p_lin; P = 2.0e8, τII = 1.2e8, EII = EII)
        F_nd = compute_yieldfunction(p_lin_nd; P = 2.0e8 / char_stress, τII = 1.2e8 / char_stress, EII = EII)
        @test F_dim / char_stress ≈ F_nd
    end

    p_cap = DruckerPragerCap(; softening_C = NonLinearSoftening(; ξ₀ = 1.5e7Pa, Δ = 7.5e6Pa), C = 1.5e7Pa, ϕ = 30.0, η_vp = 1.0e19Pa * s, pT = -1.5e6Pa, Ψ = 15.0)
    p_cap_nd = nondimensionalize(p_cap, CharDim)
    @test p_cap_nd.softening_C.ξ₀.val ≈ 1.5
    @test p_cap_nd.softening_C.Δ.val ≈ 0.75
    for EII in (0.0, 0.5, 2.0)
        F_dim = compute_yieldfunction(p_cap; P = 1.5e8, τII = 9.0e7, EII = EII)
        F_nd = compute_yieldfunction(p_cap_nd; P = 1.5e8 / char_stress, τII = 9.0e7 / char_stress, EII = EII)
        @test F_dim / char_stress ≈ F_nd
    end

    # MPa-declared parameters must nondimensionalize identically to Pa-declared ones
    p_MPa = DruckerPragerCap(; softening_C = NonLinearSoftening(; ξ₀ = 15.0MPa, Δ = 7.5MPa), C = 15.0MPa, ϕ = 30.0, η_vp = 1.0e19Pa * s, pT = -1.5MPa, Ψ = 15.0)
    p_MPa_nd = nondimensionalize(p_MPa, CharDim)
    @test p_MPa_nd.C.val ≈ p_cap_nd.C.val
    @test p_MPa_nd.pT.val ≈ p_cap_nd.pT.val
    @test p_MPa_nd.softening_C.ξ₀.val ≈ p_cap_nd.softening_C.ξ₀.val
    @test p_MPa_nd.softening_C.Δ.val ≈ p_cap_nd.softening_C.Δ.val
    for EII in (0.0, 0.5, 2.0)
        F_Pa = compute_yieldfunction(p_cap_nd; P = 15.0, τII = 9.0, EII = EII)
        F_MPa = compute_yieldfunction(p_MPa_nd; P = 15.0, τII = 9.0, EII = EII)
        @test F_Pa ≈ F_MPa
    end

    # partial positional construction must not stack-overflow
    nls2 = NonLinearSoftening(30.0, 10.0)
    @test nls2.ξ₀.val == 30.0
    @test nls2.Δ.val == 10.0

end
