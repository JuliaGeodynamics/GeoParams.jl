@testset "Transformations" begin
    geo = GEO_units(; viscosity = 1.0e19Pas, length = 1000km)
    si = SI_units()

    @test nondimensionalize(10cm / yr, geo) ≈ 0.0031688087814028945
    @test nondimensionalize(10cm / yr, si) ≈ 3.168808781402895e7
    @test nondimensionalize(10MPa / s, geo) ≈ 1.0e12
    @test nondimensionalize(10Pa^1.2 / s, geo) ≈ 39810.71705534975
    @test nondimensionalize(6.3e-2MPa^-3.05 / s, GEO_units()) ≈ 7.068716262102384e14

    velocity = nondimensionalize(3cm / yr, geo)
    @test dimensionalize(velocity, cm / yr, geo) ≈ 3cm / yr
    @test dimensionalize_and_strip(velocity, cm / yr, geo) ≈ 3
    @test udim(velocity, cm / yr, geo) ≈ 3
    @test (@dimstrip velocity cm / yr geo) ≈ 3

    wrapped = GeoUnit([100km, 200km])
    nondimensional = nondimensionalize(wrapped, geo)
    @test nondimensional.val ≈ [0.1, 0.2]
    @test !isdimensional(nondimensional)
    restored = dimensionalize(nondimensional, geo)
    @test restored.val ≈ [100, 200]
    @test isdimensional(restored)

    tuple_result = nondimensionalize((100km, 800 / s), geo)
    @test tuple_result[1] ≈ 0.1
    @test tuple_result[2] ≈ 8.0e14
    @test nondimensionalize((), geo) === ()
    @test_throws ArgumentError nondimensionalize(5.0, geo)
    @test_throws ArgumentError nondimensionalize(Float64[], geo)
    @test (@test_logs (:warn,) nondimensionalize(10cm / yr, nothing)) == 10
    @test dimensionalize(0.1, Pa, nothing) == 0.1
end
