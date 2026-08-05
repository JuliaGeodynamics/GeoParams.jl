@testset "Characteristic units" begin
    geo = GEO_units(; viscosity = 1.0e19Pas, length = 1000km)
    @test geo.length == 1000km
    @test geo.stress == 10MPa
    @test geo.Pa == 1.0e7u"Pa"
    @test geo.Mass ≈ 1.0e37u"kg"
    @test ustrip(C, geo.temperature) ≈ 1000

    si = SI_units()
    @test si.length == 1000m
    @test si.Pa == 10Pa

    none = NO_units()
    @test none.length == 1
    @test none.Pa == 1

    @test GEO_units(; length = 10).length == 10km
    @test SI_units(; stress = 20).stress == 20Pa
    @test_throws ArgumentError NO_units(; length = 1m)
    @test upgrade_GeoUnits(geo).length == geo.length
    @test sprint(show, geo) isa String
end
