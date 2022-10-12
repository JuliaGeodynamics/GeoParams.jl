using Test
using GeoParams
using GLMakie

@testset "StrengthEnvelope.jl" begin
    MatParam  = (SetMaterialParams(Name="UC", Phase=1, Density=ConstantDensity(ρ=2700kg/m^3), CreepLaws = SetDislocationCreep("Wet Quartzite | Ueda et al. (2008)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)),
                 SetMaterialParams(Name="MC", Phase=2, Density=Density=ConstantDensity(ρ=2900kg/m^3), CreepLaws = SetDislocationCreep("Plagioclase An75 | Ji and Zhao (1993)"), Plasticity = DruckerPrager(ϕ=20.0, C=10MPa)),
                 SetMaterialParams(Name="LC", Phase=3, Density=PT_Density(ρ0=2900kg/m^3, α=3e-5/K, β=1e-10/Pa), CreepLaws = SetDislocationCreep("Maryland strong diabase | Mackwell et al. (1998)"), Plasticity = DruckerPrager(ϕ=30.0, C=10MPa)));
    Thickness = [15,10,15]*km;

    τ         = StrengthEnvelope(MatParam, Thickness, LinearTemp(), "test");

    @test isapprox(τ[26], 141.0952, atol=1e-4)
    @test isapprox(τ[38],   8.5801, atol=1e-4)
    @test isapprox(τ[47], 174.5291, atol=1e-4)
    @test isapprox(τ[63],  17.1202, atol=1e-4)
    @test isapprox(τ[75], 418.0952, atol=1e-4)
end