@testset "GeoUnit show & dimensionalize" begin
    CharDim = GEO_units()
    # show: dimensional and nondimensional branches + the GeoUnits container
    @test sprint(show, GeoUnit(1.0u"Pa")) isa String   # dimensional branch
    @test sprint(show, GeoUnit(2.0)) isa String        # nondimensional branch
    @test sprint(show, CharDim) isa String
    @test superscript(1 // 2) isa String
    @test which(Unitful.superscript, (Rational{Int64},)).module === Unitful

    # dimensionalize round-trips nondimensionalize for a GeoUnit
    nd = nondimensionalize(GeoUnit(1.0e6u"Pa"), CharDim)
    dim = dimensionalize(nd, CharDim)
    @test ustrip(Value(dim)) ≈ 1.0e6 rtol = 1.0e-9
end
