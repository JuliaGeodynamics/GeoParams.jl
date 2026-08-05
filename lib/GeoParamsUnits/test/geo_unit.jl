@testset "GeoUnit conversions & equality" begin
    @test isDimensional(GeoUnit(1.0u"m"))
    @test !isDimensional(GeoUnit(1.0))
    @test !isDimensional(1.0)

    # special constructors
    @test GeoUnit(nothing).val === nothing
    @test GeoUnit(x -> 2x) isa GeoUnit          # function-valued
    @test GeoUnit{Float64}(5.0) isa GeoUnit
    @test GeoUnit(1.0) isa AbstractGeoUnit
    @test GeoUnit{Float64, typeof(NoUnits)} <: AbstractGeoUnit{Float64, typeof(NoUnits)}

    # Empty arrays do not contain enough information to determine a unit.
    @test_throws ArgumentError GeoUnit(Float64[])
    @test_throws ArgumentError GeoUnit(Quantity{Float64}[])
    @test_throws ArgumentError convert(GeoUnit, Int32[])

    # convert from Int → Float
    @test convert(GeoUnit, Int32(5)).val === 5.0f0
    @test convert(GeoUnit, Int32[1, 2, 3]).val == Float32[1, 2, 3]
    @test convert(GeoUnit, Int64[1, 2, 3]).val == Float64[1, 2, 3]
    @test convert(GeoUnit, [1, 2, 3]).val == [1.0, 2.0, 3.0]
    @test GeoUnit(Int32(5)).val === 5.0f0
    @test GeoUnit(Int32[1, 2, 3]).val == Float32[1, 2, 3]
    # Int32 Quantity array
    @test convert(GeoUnit, Float32[1.0, 2.0]u"m").val isa AbstractArray{Float32}

    # equality
    @test isequal(GeoUnit(5.0), 5.0)
    @test isequal(GeoUnit([1.0, 2.0]), [1.0, 2.0])
    @test isequal(GeoUnit(5.0), GeoUnit(5.0))

    # array interface
    ga = convert(GeoUnit, [1.0, 2.0, 3.0]u"m")
    @test size(ga) == (3,)
    @test ga[1] isa GeoUnit                               # getindex(Vararg{Int})
    @test convert(AbstractArray, ga) == ga.val

    # promote_rule + broadcast against a Quantity array
    @test promote_rule(GeoUnit, Quantity) === GeoUnit
    @test GeoUnit(2.0) .* [1.0, 2.0]u"m" == [2.0, 4.0]u"m"
end
