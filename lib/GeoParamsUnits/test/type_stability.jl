using Test
using GeoParamsUnits
using Unitful

_length(units) = units.length
_velocity(units) = units.velocity
_stress(units) = units.stress
_density(units) = units.density

@testset "Type stability" begin
    @testset "Characteristic units" begin
        geo = @inferred GEO_units()
        si = @inferred SI_units()
        none = @inferred NO_units()

        @test @inferred(_length(geo)) == 1000km
        @test @inferred(_velocity(geo)) == geo.m / geo.s
        @test @inferred(_stress(si)) == 10Pa
        @test @inferred(_density(none)) == 1
        @test all(type -> type !== Any, fieldtypes(typeof(getfield(geo, :data))))

        geo_from_numbers = @inferred GEO_units(;
            length = 1000.0,
            temperature = 1000.0,
            stress = 10.0,
            viscosity = 1.0e20,
        )
        @test @inferred(_length(geo_from_numbers)) == 1000.0km
    end

    @testset "Dimensionalization" begin
        geo = GEO_units()
        quantity = 3.0cm / yr
        quantities = [1.0, 2.0]m
        wrapped = @inferred GeoUnit(quantity)

        nondimensional = @inferred nondimensionalize(wrapped, geo)
        @test @inferred(nondimensionalize(quantity, geo)) isa Float64
        @test @inferred(nondimensionalize(quantities, geo)) isa Vector{Float64}
        @test @inferred(compute_units(wrapped, geo)) isa Quantity
        @test @inferred(dimensionalize(nondimensional, geo)) isa typeof(wrapped)
        @test @inferred(dimensionalize(0.1, cm / yr, geo)) isa Quantity
    end
end
