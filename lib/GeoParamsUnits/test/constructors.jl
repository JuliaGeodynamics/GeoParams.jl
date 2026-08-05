@testset "GeoUnit integer-quantity constructors" begin
    @test GeoUnit(Int32(2)u"Pa").val isa Float32
    @test GeoUnit(Int64(2)u"Pa").val isa Float64
    @test GeoUnit([Int32(2), Int32(3)]u"Pa").val isa Vector{Float32}
    @test GeoUnit([Int64(2), Int64(3)]u"Pa").val isa Vector{Float64}
end
