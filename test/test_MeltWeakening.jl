using Test
using GeoParams

@testset "MeltWeakening" begin
    diff0 = SetDiffusionCreep(GeoParams.Diffusion.wet_olivine_Mei_2000a)
    dis0 = SetDislocationCreep(GeoParams.Dislocation.wet_olivine_Mei_2000b)
    diff = MeltWeakening(diff0; α = 26)
    dis = MeltWeakening(dis0; α = 31)
    args = (; T = 1473.0, P = 300.0e6, d = 15.0e-6, f = 300.0e6)
    τ, ε = 1.0e6, 1.0e-14

    @test isbits(diff)
    @test diff.law === diff0
    @test startswith(sprint(show, diff), "MeltWeakening: α = 26.0")
    @test param_info(diff).Comment == param_info(diff0).Comment

    @testset "Mei et al. (2002) enhancement factors" begin
        # Abstract: at a given stress, ϕ = 0.12 deforms ≈ 20× (diffusion) and ≈ 40× (dislocation)
        # faster than melt-free samples, and ϕ = 0.02 ≈ 2× (diffusion).
        ratio(law, law0, ϕ) = compute_εII(law, τ; args..., ϕ) / compute_εII(law0, τ; args...)
        @test ratio(diff, diff0, 0.12) ≈ exp(26 * 0.12) rtol = 1.0e-12
        @test 15 < ratio(diff, diff0, 0.12) < 30
        @test 1.5 < ratio(diff, diff0, 0.02) < 2.5
        @test 30 < ratio(dis, dis0, 0.12) < 55
    end

    @testset "no melt leaves the law unchanged" begin
        @test compute_εII(diff, τ; args...) == compute_εII(diff0, τ; args...)
        @test compute_τII(dis, ε; args...) == compute_τII(dis0, ε; args...)
        @test compute_εII(MeltWeakening(dis0; α = 0), τ; args..., ϕ = 0.1) == compute_εII(dis0, τ; args...)
    end

    @testset "ϕ reaches the wrapped law" begin
        nested = MeltWeakening(diff; α = 5)
        @test compute_εII(nested, τ; args..., ϕ = 0.1) / compute_εII(diff0, τ; args...) ≈ exp(31 * 0.1) rtol = 1.0e-12
    end

    @testset "inverse and derivatives" begin
        a = (; args..., ϕ = 0.05)
        for law in (diff, dis)
            @test compute_τII(law, compute_εII(law, τ; a...); a...) ≈ τ rtol = 1.0e-10
            @test dεII_dτII(law, τ; a...) ≈ exp(law.α.val * 0.05) * dεII_dτII(law.law, τ; args...) rtol = 1.0e-12
            τε = compute_τII(law, ε; a...)
            @test dτII_dεII(law, ε; a...) * dεII_dτII(law, τε; a...) ≈ 1 rtol = 1.0e-8
        end
    end

    @testset "composite rheology" begin
        c0 = CompositeRheology(diff0, dis0)
        c = CompositeRheology(diff, dis)
        @test compute_viscosity_εII(c, ε, args) ≈ compute_viscosity_εII(c0, ε, args) rtol = 1.0e-10
        η0 = compute_viscosity_εII(c0, ε, args)
        η = compute_viscosity_εII(c, ε, (; args..., ϕ = 0.05))
        # at fixed strain rate each mechanism weakens by e^{-αϕ/n}; the composite lies in between
        @test exp(-31 * 0.05 / 3.0) > η / η0 > exp(-26 * 0.05 / 1.1) * (1 - 1.0e-8)
    end

    @testset "arrays" begin
        ϕ = [0.0, 0.05, 0.12]
        E = zeros(3)
        compute_εII!(E, diff, fill(τ, 3); T = fill(args.T, 3), P = args.P, d = args.d, f = args.f, ϕ)
        @test E ≈ exp.(26 .* ϕ) .* compute_εII(diff0, τ; args...) rtol = 1.0e-12
        S = zeros(3)
        compute_τII!(S, diff, E; args..., ϕ)
        @test S ≈ fill(τ, 3) rtol = 1.0e-10
    end

    @testset "nondimensional and Float32" begin
        CD = GEO_units()
        nd = nondimensionalize(diff, CD)
        @test nd isa MeltWeakening
        @test isDimensional(diff)
        @test !isDimensional(nd)
        a_nd = (;
            T = nondimensionalize(args.T * K, CD), P = nondimensionalize(args.P * Pa, CD),
            d = nondimensionalize(args.d * m, CD), f = nondimensionalize(args.f * Pa, CD), ϕ = 0.05,
        )
        τ_nd = nondimensionalize(τ * Pa, CD)
        ε_nd = compute_εII(nd, τ_nd; a_nd...)
        @test ustrip(dimensionalize(ε_nd, s^-1, CD)) ≈ compute_εII(diff, τ; args..., ϕ = 0.05) rtol = 1.0e-10
        nd32 = convert_precision(Float32, nd)
        ε32 = compute_εII(nd32, Float32(τ_nd); map(Float32, a_nd)...)
        @test ε32 isa Float32
        @test ε32 ≈ ε_nd rtol = 1.0e-5
    end

    @testset "invalid α" begin
        @test_throws "α must be non-negative" MeltWeakening(diff0; α = -26)
        @test_throws "α must be dimensionless" MeltWeakening(diff0; α = 26s)
    end
end
