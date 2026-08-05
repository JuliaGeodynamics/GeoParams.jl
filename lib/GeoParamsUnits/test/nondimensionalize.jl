@testset "nondimensionalize dispatch variants" begin
    CharDim = GEO_units()
    # quantity / array / tuple with characteristic units
    @test nondimensionalize(10.0u"Pa", CharDim) isa Number
    @test nondimensionalize([1.0, 2.0]u"Pa", CharDim) isa AbstractArray
    @test nondimensionalize((1.0u"Pa", 2.0u"Pa"), CharDim) isa Tuple
    @test nondimensionalize((), CharDim) === ()
    # no characteristic units -> returned (with a warning) unchanged
    @test nondimensionalize(GeoUnit(1.0u"Pa"), nothing) isa GeoUnit
    @test nondimensionalize([1.0, 2.0], nothing) isa AbstractArray
    # missing units where they are required -> throws
    @test_throws ArgumentError nondimensionalize(1.0, CharDim)
    @test_throws ArgumentError nondimensionalize((1.0, 2.0), CharDim)
    # float array without units where units are required -> throws
    @test_throws ArgumentError nondimensionalize([1.0, 2.0], CharDim)
    # GeoUnit show + indexing
    g = GeoUnit([1.0 2.0; 3.0 4.0]u"Pa")
    @test sprint(show, g) isa String
    @test g[1, 2] isa GeoUnit
end
