using Test
using GeoParams

# Exercises the dimensional (`@unpack_units`) input paths and conditional branches
# that the plain-float unit tests skip.
@testset "Units & branching" begin

    # ---- Density: Quantity paths + conduit branches ----
    @testset "Density" begin
        @test compute_density(Compressible_Density(), (; P = 1.0e9Pa)) isa Quantity
        @test compute_density(T_Density(), (; T = 1000.0K)) isa Quantity

        # BubbleFlow: hit both `P < cutoff` (c = a√P) and `P ≥ cutoff` (c = c0) branches
        bf = BubbleFlow_Density(c0 = 0.01)
        cutoff = (NumValue(bf.c0))^2 / (NumValue(bf.a))^2
        @test compute_density(bf, (; P = cutoff / 2)) isa Number          # P < cutoff
        @test compute_density(bf, (; P = cutoff * 2)) isa Number          # P ≥ cutoff
        @test compute_density(bf, (; P = (cutoff / 2)Pa)) isa Number      # Quantity input branch

        # get_α for conduit densities (needs α-bearing sub-densities)
        bfα = BubbleFlow_Density(ρmelt = PT_Density(α = 1.0e-3), ρgas = PT_Density(α = 1.0e-2), c0 = 0.01)
        @test GeoParams.get_α(bfα, (; P = cutoff / 2)) isa Number
        @test GeoParams.get_α(bfα, (; P = cutoff * 2)) isa Number
    end

    # ---- Conductivity: Quantity paths + low/high-T & d≠0 branches ----
    @testset "Conductivity" begin
        wh = T_Conductivity_Whittington()
        @test compute_conductivity(wh, (; T = 500.0)) isa Number          # low-T branch (T ≤ Tcutoff)
        @test compute_conductivity(wh, (; T = 1000.0)) isa Number         # high-T branch
        @test compute_conductivity(wh, (; T = 500.0K)) isa Quantity       # @unpack_units

        whp = T_Conductivity_Whittington_parameterised()
        @test compute_conductivity(whp, (; T = 600.0K)) isa Quantity

        mantle = Set_TP_Conductivity("Mantle")                            # d ≠ 0
        @test compute_conductivity(mantle, (; P = 1.0e9, T = 1000.0)) isa Number
        @test compute_conductivity(mantle, (; P = 1.0e9Pa, T = 1000.0K)) isa Quantity
    end

    # ---- HeatCapacity: low-T branch ----
    @testset "HeatCapacity" begin
        cp = T_HeatCapacity_Whittington()
        @test compute_heatcapacity(cp; T = 500.0) isa Number              # T ≤ Tcutoff branch
        @test compute_heatcapacity(cp; T = 1000.0) isa Number
    end

    # ---- CreepLaw: CorrectionFactor branches ----
    @testset "CreepLaw" begin
        @test GeoParams.CorrectionFactor(GeoParams.SimpleShear) == (2.0, 2.0)
        @test GeoParams.CorrectionFactor(GeoParams.Invariant) == (1.0, 1.0)
        @test_throws ArgumentError GeoParams.CorrectionFactor(99)
    end
end
