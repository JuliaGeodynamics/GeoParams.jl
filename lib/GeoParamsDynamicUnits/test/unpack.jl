@testset "Unpacking" begin
    parameters = (; density = GeoUnit(3000kg / m^3), length = GeoUnit(10km))
    @unpack_val density, length = parameters
    @test density == 3000
    @test length == 10

    @unpack_units density, length = parameters
    @test density == 3000kg / m^3
    @test length == 10km
    @test unpack_vals((parameters.density, parameters.length)) == (3000, 10)
end
