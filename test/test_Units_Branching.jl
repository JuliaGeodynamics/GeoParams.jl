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

    # ---- Units: edge-case constructors and helpers ----
    @testset "Units edges" begin
        # upgrade_GeoUnits recreates a GeoUnits from its fields
        g = GEO_units(length = 500km, temperature = 800C, stress = 5MPa, viscosity = 1.0e19Pas)
        g2 = GeoParams.upgrade_GeoUnits(g)
        @test g2.length == g.length

        # GeoUnit setindex! (mutable array path)
        gu = convert(GeoUnit, [1.0, 2.0, 3.0])   # dimensionless array
        gu[2] = 99.0
        @test gu.val[2] == 99.0

        # GeoUnit from Int32 array of Quantities
        gu32 = convert(GeoUnit, convert.(Float32, [1.0, 2.0, 3.0]) * u"m")
        @test gu32.val isa AbstractArray{Float32}

        # broadcasted ^ with GeoUnit (covers Base.broadcasted ^)
        gu_b = GeoUnit(2.0)
        result = gu_b .^ [2, 3]
        @test result ≈ [4.0, 8.0]

        # broadcasted +/- on dimensional GeoUnit arrays (hits the isdimensional_new branch)
        gd1 = convert(GeoUnit, [1.0, 2.0, 3.0] * u"m")
        gd2 = convert(GeoUnit, [4.0, 5.0, 6.0] * u"m")
        sum_d = gd1 .+ gd2
        @test NumValue(sum_d) ≈ [5.0, 7.0, 9.0]
        diff_d = gd2 .- gd1
        @test NumValue(diff_d) ≈ [3.0, 3.0, 3.0]

        # GEO_units / SI_units accept plain (unitless) numbers -> unit-coercion branches
        g_plain = GEO_units(; length = 1000, temperature = 1000, stress = 10, viscosity = 1.0e20)
        @test g_plain.length == 1000km
        s_plain = SI_units(; length = 1000, temperature = 1000, stress = 10, viscosity = 1.0e20)
        @test s_plain.length == 1000m

        # NO_units rejects united arguments (error branches)
        @test_throws ErrorException NO_units(; temperature = 1K)
        @test_throws ErrorException NO_units(; length = 1m)
        @test_throws ErrorException NO_units(; stress = 1Pa)
        @test_throws ErrorException NO_units(; viscosity = 1Pas)

        # superscript override: rational exponent display path
        @test GeoParams.Unitful.superscript(1 // 2) isa String
    end

    # ---- Computations: phase-integer and SVector phase-ratio paths ----
    @testset "Computations paths" begin
        using StaticArrays
        p1 = SetMaterialParams(; Name = "A", Phase = 1, Density = ConstantDensity(ρ = 2900kg / m^3))
        p2 = SetMaterialParams(; Name = "B", Phase = 2, Density = ConstantDensity(ρ = 3100kg / m^3))
        phases_tup = (p1, p2)

        # compute_param with NTuple + integer phase (generated function)
        ρ1 = compute_density(phases_tup, 1, (;))
        @test ρ1 ≈ 2900.0

        # compute_param with SVector phase ratios (generated function)
        ratio = SA[0.4, 0.6]
        ρ_mix = compute_density(phases_tup, ratio, (;))
        @test ρ_mix ≈ 0.4 * 2900.0 + 0.6 * 3100.0

        # nphase fallback: phase index present in no material -> 0.0
        @test compute_density(phases_tup, 99, (;)) == 0.0

        # compute_param! with a PhaseRatio array of the wrong rank -> error guard
        rho2d = zeros(4, 4)                # N = 2 -> PhaseRatio must be 3D
        badPR = zeros(4, 4)                # M = 2 ≠ N+1
        @test_throws ErrorException GeoParams.compute_density!(rho2d, phases_tup, badPR, (;))
    end

    # ---- aliases.jl: error paths ----
    @testset "aliases errors" begin
        # too few arguments
        @test_throws ErrorException GeoParams.checkargs(:GeoParamsAliases)
        # wrong first argument
        @test_throws ErrorException GeoParams.checkargs(:WrongName, :(density = ρ))
        # non-keyword second argument (not an `=` expression)
        @test_throws ErrorException GeoParams.checkargs(:GeoParamsAliases, :density)
        # invalid keyword argument key
        @test_throws ErrorException GeoParams.validate_kwargkeys(
            Dict(:invalid_key => :something), "@use GeoParamsAliases"
        )
        # create_module on an already-defined module name -> error
        @test_throws ErrorException GeoParams.create_module(Main, :Base, :(density = ρ))
    end

    # ---- Density: dimensional get_α path for conduit densities ----
    @testset "get_α Quantity path" begin
        bf = BubbleFlow_Density(ρmelt = PT_Density(α = 1.0e-3), ρgas = PT_Density(α = 1.0e-2), c0 = 0.01)
        @test GeoParams.get_α(bf, (; P = 1.0Pa)) isa Number   # P::Quantity -> @unpack_units branch
    end

    @testset "Density / Viscosity / Conductivity routines" begin
        # ConstantDensity in-place array (default ρ = 2900 kg/m³)
        ρ = zeros(5)
        compute_density!(ρ, ConstantDensity(); P = zeros(5), T = zeros(5))
        @test all(ρ .== 2900.0)
        compute_density!(ρ, ConstantDensity(), (; P = zeros(5), T = zeros(5)))
        @test all(ρ .== 2900.0)

        # Vector_Density show
        vd = Vector_Density(; rho = collect(1.0:5.0))
        @test sprint(show, vd) isa String

        # visco-elastic viscosity helper — regression value
        @test compute_elastoviscosity(ConstantElasticity(), 1.0e20, 1.0e3) ≈ 4.99999750000125e13 rtol = 1.0e-9
        @test compute_elastoviscosity(ConstantElasticity(), 1.0e20, (; dt = 1.0e3)) ≈ 4.99999750000125e13 rtol = 1.0e-9

        # constant conductivity in-place (default k = 3 W/m/K)
        k = zeros(4)
        compute_conductivity!(k, ConstantConductivity())
        @test all(k .== 3.0)
        @test compute_conductivity(ConstantConductivity()) == 3.0
    end

    @testset "conductivity / seismic / heat-capacity routines" begin
        # temperature-dependent conductivity over an array (functor + positional args)
        Tarr = collect(300.0:100.0:800.0)
        for s in (T_Conductivity_Whittington(), T_Conductivity_Whittington_parameterised())
            k = s(Tarr)
            @test length(k) == length(Tarr)
            @test all(>(0), k)
            @test param_info(s) isa MaterialParamsInfo
            @test sprint(show, s) isa String
        end
        # constant seismic velocity show + compute
        sv = ConstantSeismicVelocity()
        @test sprint(show, sv) isa String

        # latent heat capacity construction + show
        lhc = Latent_HeatCapacity()
        @test sprint(show, lhc) isa String
        @test compute_heatcapacity(lhc; T = 1000.0) isa Number
    end

end
