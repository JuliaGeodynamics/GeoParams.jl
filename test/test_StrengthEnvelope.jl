using Test
using GeoParams
import GeoParams.Dislocation

@testset "StrengthEnvelope.jl" begin
    MatParam = (
        SetMaterialParams(;
            Name = "UC",
            Phase = 1,
            Density = ConstantDensity(; ρ = 2700kg / m^3),
            CreepLaws = SetDislocationCreep(Dislocation.wet_quartzite_Ueda_2008),
            Plasticity = DruckerPrager(; ϕ = 30.0, C = 10MPa),
        ),
        SetMaterialParams(;
            Name = "MC",
            Phase = 2,
            Density = Density = ConstantDensity(; ρ = 2900kg / m^3),
            CreepLaws = SetDislocationCreep(Dislocation.plagioclase_An75_Ji_1993),
            Plasticity = DruckerPrager(; ϕ = 20.0, C = 10MPa),
        ),
        SetMaterialParams(;
            Name = "LC",
            Phase = 3,
            Density = PT_Density(; ρ0 = 2900kg / m^3, α = 3.0e-5 / K, β = 1.0e-10 / Pa),
            CreepLaws = SetDislocationCreep(
                Dislocation.strong_diabase_Mackwell_1998
            ),
            Plasticity = DruckerPrager(; ϕ = 30.0, C = 10MPa),
        ),
    )
    Thickness = [15, 10, 15] * km

    # Use default linear gradient and default strain rate
    z, τ, T = StrengthEnvelopeComp(MatParam, Thickness)

    @test isapprox(τ[26].val, 141.0952, atol = 1.0e-4)
    @test isapprox(τ[38].val, 8.5801, atol = 1.0e-4)
    @test isapprox(τ[47].val, 174.5291, atol = 1.0e-4)
    @test isapprox(τ[63].val, 17.1202, atol = 1.0e-4)
    @test isapprox(τ[75].val, 418.0952, atol = 1.0e-4)
    @test isapprox(z[34].val, 13.2, atol = 1.0e-4)
    @test isapprox(T[89].val, 704.0, atol = 1.0e-4)

    # Use HalfspaceCool profile and user-defined strain rate
    z, τ, T = StrengthEnvelopeComp(
        MatParam,
        Thickness,
        HalfspaceCoolTemp(0C, 1350C, 10Myr, 0K / km, 1.0e-6m^2 / s),
        1.0e-10 / s,
    )

    @test isapprox(τ[25].val, 135.7978, atol = 1.0e-4)
    @test isapprox(τ[38].val, 9.32, atol = 1.0e-4)
    @test isapprox(τ[40].val, 113.4964, atol = 1.0e-4)
    @test isapprox(τ[55].val, 19.3755, atol = 1.0e-4)
    @test isapprox(τ[85].val, 23.1451, atol = 1.0e-4)
end
