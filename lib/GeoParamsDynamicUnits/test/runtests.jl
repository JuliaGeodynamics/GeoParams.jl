using Test
using GeoParamsDynamicUnits
using DynamicQuantities

@testset "GeoParamsDynamicUnits" begin
    include("api.jl")
    include("geo_unit.jl")
    include("characteristic_units.jl")
    include("transformations.jl")
    include("arithmetic.jl")
    include("unpack.jl")
    include("type_stability.jl")
    include("quality.jl")
end
