using Test
using GeoParams

@testset "DiffusionCreepGrainSize" begin
    x = DiffusionCreepGrainSize(; Name = "test", d0 = 1mm)

    @test x isa AbstractCreepLaw
    @test Value(x.n) == 1.0
    @test Value(x.r) == 0.0
    @test Value(x.p) == -3.0
    @test Value(x.rg) == 1.38
    @test x.Apparatus == AxialCompression
    @test remove_tensor_correction(x).Apparatus == Invariant

    εII = 1.0e-15
    τII = compute_τII(x, εII; T = 1000.0)
    @test compute_εII(x, τII; T = 1000.0) ≈ εII

    growth_rate = compute_growth_rate(x, τII; T = 1000.0, d = 1.0e-3)
    reduction_rate = compute_grain_size_reduction(x, τII; εII = εII, d = 1.0e-3)
    @test isfinite(growth_rate)
    @test growth_rate > 0
    @test isfinite(reduction_rate)
    @test reduction_rate > 0

    d = 1.0e-3
    dt = 1.5e12
    d_growth = compute_growth_rate(x, τII; T = 1000.0, d = d)
    d_red = compute_grain_size_reduction(x, τII; εII = εII, d = d)
    d_new = (d_growth - d_red) * dt

    @test d_new ≈ 0.225376218760203

    τII_units = compute_τII(x, εII / s; T = 1000K)
    εII_units = compute_εII(x, τII_units; T = 1000K)
    @test ustrip(uconvert(s^-1, εII_units)) ≈ εII
end
