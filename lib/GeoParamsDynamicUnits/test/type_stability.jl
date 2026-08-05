_length(units) = units.length
_velocity(units) = units.velocity

@testset "Type stability" begin
    geo = @inferred GEO_units()
    si = @inferred SI_units()
    none = @inferred NO_units()

    @test @inferred(_length(geo)) == 1000km
    @test @inferred(_velocity(geo)) == geo.m / geo.s
    @test @inferred(_length(si)) == 1000m
    @test @inferred(_length(none)) == 1

    quantity = 3.0cm / yr
    wrapped = @inferred GeoUnit(quantity)
    nondimensional = @inferred nondimensionalize(wrapped, geo)
    @test @inferred(nondimensionalize(quantity, geo)) isa Float64
    @test @inferred(nondimensionalize([1m, 2m], geo)) isa Vector{Float64}
    @test @inferred(compute_units(wrapped, geo)) isa DynamicQuantities.Quantity
    @test @inferred(dimensionalize(nondimensional, geo)) isa GeoUnit
    @test @inferred(dimensionalize(0.1, cm / yr, geo)) isa DynamicQuantities.Quantity
end
