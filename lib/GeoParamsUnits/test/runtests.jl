using Test
using GeoParamsUnits
using Unitful

@testset "GeoParamsUnits" begin
    include("geo_unit.jl")
    include("units.jl")
    include("constructors.jl")
    include("nondimensionalize.jl")
    include("display.jl")
    include("arithmetic.jl")
    include("dimensionalize.jl")
    include("type_stability.jl")
    include("quality.jl")
end
