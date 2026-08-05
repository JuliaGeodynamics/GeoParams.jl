@testset "dimensionalize (value, unit, CharDim) forms" begin
    CharDim = GEO_units()
    # 3-arg dimensionalize of a bare nondimensional value back to a unit
    @test dimensionalize(0.1, u"Pa", CharDim) ≈ 1.0e6u"Pa" rtol = 1.0e-12
    @test udim(0.1, u"Pa", CharDim) ≈ 1.0e6 rtol = 1.0e-12
    # no characteristic units -> returned (with a warning) un-dimensionalized
    @test dimensionalize(0.1, u"Pa", nothing) == 0.1
end
