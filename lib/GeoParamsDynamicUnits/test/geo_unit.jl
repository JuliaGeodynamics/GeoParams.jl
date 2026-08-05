@testset "GeoUnit" begin
    dimensional = GeoUnit(100km)
    @test dimensional isa AbstractGeoUnit
    @test dimensional.val == 100
    @test Unit(dimensional) == km
    @test isdimensional(dimensional)
    @test Value(dimensional) == 100km

    plain = GeoUnit(100)
    @test plain.val == 100.0
    @test Unit(plain) == 1
    @test !isDimensional(plain)

    array = GeoUnit([100km 200km; 300km 400km])
    @test array.val == [100 200; 300 400]
    @test size(array) == (2, 2)
    @test array[1, 2] isa GeoUnit
    @test collect(GeoUnit([1km, 2km])) isa Vector{<:GeoUnit}

    @test GeoUnit(Int32(2)).val isa Float32
    @test GeoUnit(Int64(2)).val isa Float64
    @test_throws ArgumentError GeoUnit(Float64[])
    @test GeoUnit(nothing).val === nothing
    @test GeoUnit(x -> 2x) isa GeoUnit

    @test convert(GeoUnit, Int32[1, 2]).val isa Vector{Float32}
    @test convert(GeoUnit{Float32}, 3.0).val isa Float32
    @test isequal(GeoUnit(5.0), 5.0)
    @test repr("text/plain", dimensional) isa String
end
