@testset "GeoUnit arithmetic & broadcasting" begin
    G = GeoUnit
    gu = G(2.0u"Pa")
    gu2 = G(3.0)
    # GeoUnit <op> Quantity (returns a Quantity)
    @test gu * 3.0u"Pa" == 6.0u"Pa^2"
    @test 3.0u"Pa" * gu == 6.0u"Pa^2"
    # GeoUnit <op> plain array -> values only (no units)
    @test gu2 * [1.0, 2.0] == [3.0, 6.0]
    @test [1.0, 2.0] * gu2 == [3.0, 6.0]
    # GeoUnit broadcast with a Quantity array -> keeps units
    @test gu2 .* [1.0u"Pa", 2.0u"Pa"] == [3.0u"Pa", 6.0u"Pa"]
    @test [1.0u"Pa", 2.0u"Pa"] .* gu2 == [3.0u"Pa", 6.0u"Pa"]
end
