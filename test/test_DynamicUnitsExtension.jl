using Test
using GeoParams
using GeoParamsDynamicUnits

const GPU = GeoParams.Units
const GPD = GeoParamsDynamicUnits
const DQ = GeoParamsDynamicUnits.DynamicQuantities

@testset "GeoParamsDynamicUnits extension" begin
    @test Base.get_extension(GeoParams, :GeoParamsDynamicUnitsExt) !== nothing

    dynamic_density = 2900GPD.kg / GPD.m^3
    density = ConstantDensity(; ρ = dynamic_density)
    @test density.ρ isa GPU.GeoUnit
    @test GPU.ustrip(GPU.upreferred(GPU.Value(density.ρ))) == 2900

    dynamic_geo = GPD.GeoUnit(3GPD.cm / GPD.yr)
    unitful_geo = convert(GPU.GeoUnit, dynamic_geo)
    roundtrip_geo = convert(GPD.GeoUnit, unitful_geo)
    @test unitful_geo.isdimensional
    @test isapprox(
        DQ.uexpand(GPD.Value(roundtrip_geo)),
        DQ.uexpand(GPD.Value(dynamic_geo)),
    )

    fractional_geo = GPD.GeoUnit(2.5GPD.MPa^(-3.05) / GPD.s)
    fractional_roundtrip = convert(GPD.GeoUnit, convert(GPU.GeoUnit, fractional_geo))
    @test isapprox(
        DQ.uexpand(GPD.Value(fractional_roundtrip)),
        DQ.uexpand(GPD.Value(fractional_geo));
        rtol = 1.0e-12,
    )

    celsius = convert(GPU.GeoUnit, GPD.GeoUnit(100GPD.C))
    @test isapprox(GPU.ustrip(GPU.upreferred(GPU.Value(celsius))), 373.15)

    dynamic_units = GPD.GEO_units()
    unitful_units = convert(GPU.GeoUnits, dynamic_units)
    @test unitful_units isa GPU.GeoUnits{GPU.GEO}
    @test isapprox(
        GPU.ustrip(GPU.upreferred(unitful_units.length)),
        DQ.ustrip(DQ.uexpand(dynamic_units.length));
        rtol = 1.0e-12,
    )

    dynamic_value = GeoParams.nondimensionalize(3GPD.cm / GPD.yr, dynamic_units)
    @test isapprox(
        dynamic_value,
        GPU.nondimensionalize(3GPU.cm / GPU.yr, GPU.GEO_units());
        rtol = 1.0e-12,
    )
    @test GeoParams.dimensionalize(1.0, GPD.MPa, dynamic_units) isa DQ.AbstractQuantity

    nondimensional_density = GeoParams.nondimensionalize(density, dynamic_units)
    @test !nondimensional_density.ρ.isdimensional
    dimensional_density = GeoParams.dimensionalize(nondimensional_density, dynamic_units)
    @test dimensional_density.ρ.isdimensional
    @test isapprox(dimensional_density.ρ.val, density.ρ.val; rtol = 1.0e-12)

    unitful_units_roundtrip = convert(GPD.GeoUnits, unitful_units)
    @test unitful_units_roundtrip isa GPD.GeoUnits{GPD.GEO}
end
