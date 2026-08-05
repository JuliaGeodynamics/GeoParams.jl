@testset "Arithmetic" begin
    a = GeoUnit(3km)
    b = GeoUnit(2000m)
    @test NumValue(a - b) ≈ 1000
    @test Unit(a - b) == m

    temperature = GeoUnit(400K)
    reference = GeoUnit(293K)
    @test NumValue(temperature + reference) ≈ 693
    @test NumValue(temperature - reference) ≈ 107
    @test Unit(temperature * reference) == K^2
    @test NumValue(temperature / reference) ≈ 400 / 293

    plain = GeoUnit(2.0)
    @test plain .* [1.0, 2.0] == [2.0, 4.0]
    @test plain .* [1m, 2m] == [2m, 4m]
    @test plain .^ [2, 3] == [4.0, 8.0]

    array = GeoUnit([1m, 2m])
    array[2] = 3
    @test array.val == [1, 3]
end
