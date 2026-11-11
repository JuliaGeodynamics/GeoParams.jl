using Test, GeoParams

@testset "Fragmentation" begin
    @testset "CriticalPorosityFragmentation" begin
        f = CriticalPorosityFragmentation()
        @test f.ϕ_crit == 0.75
        @test compute_fragmentation(f, (; ϕ_gas = 0.7)) < 1
        @test compute_fragmentation(f, (; ϕ_gas = 0.75)) == 1
        @test compute_fragmentation(f, (; ϕ_gas = 0.8)) > 1
        @test compute_fragmentation(f, (; ϕ_gas = 0.6, T = 1.0e3)) == f(; ϕ_gas = 0.6) ≈ 0.8
        @test compute_fragmentation(f, (; ϕ_gas = 0.6f0)) isa Float32
        @test isbits(f)

        @test_throws "ϕ_crit must lie in (0, 1)" CriticalPorosityFragmentation(ϕ_crit = 1)
        @test_throws "ϕ_crit must lie in (0, 1)" CriticalPorosityFragmentation(0)
    end

    @testset "StrainRateFragmentation" begin
        f = StrainRateFragmentation()
        @test f.k == 0.01
        @test UnitValue(f.G) == 1.0e10Pa
        # ratio 1 at εII = (√3/2) k G / η
        η = 1.0e8
        εc = sqrt(3) / 2 * 0.01 * 1.0e10 / η
        @test compute_fragmentation(f, (; εII = 0.99εc, η)) < 1
        @test compute_fragmentation(f, (; εII = εc, η)) ≈ 1
        @test compute_fragmentation(f, (; εII = 1.01εc, η)) > 1
        @test compute_fragmentation(f, (; εII = εc, η, T = 1.0e3)) == f(; εII = εc, η)
        @test compute_fragmentation(f, (; εII = 1.0f0, η = 1.0f8)) isa Float32
        @test isbits(f)

        CharDim = GEO_units()
        f_nd = nondimensionalize(f, CharDim)
        for εII in (0.5εc, 2εc)
            r = compute_fragmentation(f, (; εII, η))
            args_nd = (; εII = nondimensionalize(εII / s, CharDim), η = nondimensionalize(η * Pa * s, CharDim))
            @test compute_fragmentation(f_nd, args_nd) ≈ r
        end
        @test UnitValue(dimensionalize(f_nd, CharDim).G) ≈ 1.0e10Pa

        # equivalent von Mises yield surface: F = 0 exactly where the ratio is 1
        dp = DruckerPrager_regularised(f; η_vp = 1.0e5Pa * s)
        @test dp isa DruckerPrager_regularised
        @test dp.sinϕ.val == 0
        @test UnitValue(dp.C) ≈ sqrt(3) * 0.01 * 1.0e10Pa
        @test compute_yieldfunction(dp, (; τII = 2 * εc * η, P = 1.0e8)) ≈ 0 atol = 1.0e-6

        # G in GPa is converted to Pa; unit-carrying inputs give a dimensionless ratio
        @test compute_fragmentation(StrainRateFragmentation(G = 10GPa), (; εII = εc, η)) ≈ 1
        @test compute_fragmentation(f, (; εII = εc / s, η = η * Pa * s)) ≈ 1

        @test_throws "k must be positive" StrainRateFragmentation(k = 0)
        @test_throws "G must be positive" StrainRateFragmentation(G = -1Pa)
        @test_throws "G must carry stress units" StrainRateFragmentation(G = 1.0e10)
    end
end
