using Test
using GeoParams

@testset "Suspension" begin
    melt = GiordanoMeltViscosity()
    ER = Roscoe1952_Crystals(ϕ_max = 0.6)
    IZ = IshiiZuber1979_Bubbles()
    law = Suspension(melt; crystals = ER, bubbles = IZ)
    T, ϕ_x, ϕ_gas = 1100.0, 0.3, 0.4
    args = (; T, ϕ_x, ϕ_gas)
    ε, τ = 1.0e-3, 1.0e6
    # bubbles fill the fraction ϕ_gas / (1 - ϕ_x) of the bubbly melt (Truby et al., 2015, eq. 3.2)
    ηr = (1 - ϕ_x / 0.6)^-2.5 / (1 - ϕ_gas / (1 - ϕ_x))

    @test isbits(law)
    @test law.law === melt
    @test Roscoe1952_Crystals().ϕ_max == 1 / 1.35
    @test Roscoe1952_Crystals().n == 2.5
    @test Roscoe1952_Crystals(0.6f0, 2) isa Roscoe1952_Crystals{Float32}
    @test startswith(sprint(show, law), "Suspension: crystals = Roscoe1952_Crystals")
    @test param_info(law).Comment == param_info(melt).Comment

    @testset "relative viscosities" begin
        rv = compute_relative_viscosity
        # Einstein (1906) limit: η_r ≈ 1 + 2.5ϕ for dilute spheres
        @test rv(Roscoe1952_Crystals(ϕ_max = 1), 1.0e-6) ≈ 1 + 2.5e-6 rtol = 1.0e-10
        @test rv(IZ, 1.0e-6) ≈ 1 + 1.0e-6 rtol = 1.0e-10
        @test rv(Roscoe1952_Crystals(ϕ_max = 1), 0.2) ≈ 0.8^-2.5 rtol = 1.0e-14
        @test rv(IZ, 0.5) == 2
        # Keller & Suckale (2019): Krieger–Dougherty form with ϕ_max = 0.58, B = 4
        @test rv(Roscoe1952_Crystals(ϕ_max = 0.58, n = 0.58 * 4), 0.29) ≈ 0.5^-2.32 rtol = 1.0e-14
        @test rv(ER, 0.6) == Inf
        @test rv(IZ, 1.0) == Inf
        for ϕ in (-0.1, 0.61)
            @test isnan(rv(ER, ϕ))
        end
        for ϕ in (-0.1, 1.1)
            @test isnan(rv(IZ, ϕ))
        end
        @test rv(IZ, 0) === 1.0
        @test rv(law; ϕ_x, ϕ_gas) ≈ ηr rtol = 1.0e-14
        @test rv(law, (; T, ϕ_x, ϕ_gas)) == rv(law; ϕ_x, ϕ_gas)
        @test rv(Suspension(melt); ϕ_x, ϕ_gas) === 1.0
        @test_throws UndefKeywordError rv(law; ϕ_gas)
        @test_throws UndefKeywordError compute_τII(law, ε; T, ϕ_gas)
    end

    @testset "Newtonian melt: η = η_melt η_r" begin
        η_melt = compute_τII(melt, ε; T) / (2ε)
        @test compute_τII(law, ε; args...) ≈ 2 * η_melt * ηr * ε rtol = 1.0e-12
        @test compute_εII(law, τ; args...) ≈ τ / (2 * η_melt * ηr) rtol = 1.0e-12
        @test compute_τII(law, ε, args) == compute_τII(law, ε; args...)
        @test compute_viscosity_εII(law, ε, args) ≈ η_melt * ηr rtol = 1.0e-12
    end

    @testset "no suspended phases leaves the law unchanged" begin
        @test compute_εII(law, τ; T, ϕ_x = 0, ϕ_gas = 0) == compute_εII(melt, τ; T)
        @test compute_τII(Suspension(melt), ε; args...) == compute_τII(melt, ε; T)
        @test compute_τII(Suspension(melt; bubbles = IZ), ε; args...) ≈ compute_τII(melt, ε; T) / (1 - ϕ_gas / (1 - ϕ_x)) rtol = 1.0e-12
    end

    @testset "nonlinear law: inverse and derivatives" begin
        dis = Suspension(SetDislocationCreep(GeoParams.Dislocation.wet_olivine_Mei_2000b); crystals = ER, bubbles = IZ)
        a = (; T = 1473.0, P = 300.0e6, f = 300.0e6, ϕ_x, ϕ_gas)
        @test compute_τII(dis, compute_εII(dis, τ; a...); a...) ≈ τ rtol = 1.0e-10
        @test dεII_dτII(dis, τ; a...) ≈ dεII_dτII(dis.law, τ; a...) / ηr rtol = 1.0e-12
        ε_dis = 1.0e-14
        τ_dis = compute_τII(dis, ε_dis; a...)
        @test dτII_dεII(dis, ε_dis; a...) * dεII_dτII(dis, τ_dis; a...) ≈ 1 rtol = 1.0e-8
    end

    @testset "composite rheology" begin
        c = CompositeRheology(law, ConstantElasticity())
        η_melt = compute_τII(melt, ε; T) / (2ε)
        @test compute_viscosity_εII(CompositeRheology(law), ε, args) ≈ η_melt * ηr rtol = 1.0e-10
        @test isfinite(compute_τII(c, ε, (; args..., τII_old = 0.0, dt = 1.0)))
    end

    @testset "arrays" begin
        ϕx = [0.0, 0.3, 0.5]
        E = zeros(3)
        compute_εII!(E, law, fill(τ, 3); T, ϕ_x = ϕx, ϕ_gas)
        η_melt = compute_τII(melt, ε; T) / (2ε)
        @test E ≈ τ ./ (2η_melt .* (1 .- ϕx ./ 0.6) .^ -2.5 ./ (1 .- ϕ_gas ./ (1 .- ϕx))) rtol = 1.0e-12
        S = zeros(3)
        compute_τII!(S, law, E; T, ϕ_x = ϕx, ϕ_gas)
        @test S ≈ fill(τ, 3) rtol = 1.0e-10
    end

    @testset "nondimensional and Float32" begin
        CD = GEO_units()
        lin = Suspension(LinearViscous(η = 1.0e5Pa * s); crystals = ER, bubbles = IZ)
        nd = nondimensionalize(lin, CD)
        @test nd isa Suspension
        @test nd.crystals === ER
        ε_nd = nondimensionalize(ε / s, CD)
        τ_nd = compute_τII(nd, ε_nd; ϕ_x, ϕ_gas)
        @test ustrip(dimensionalize(τ_nd, Pa, CD)) ≈ compute_τII(lin, ε; ϕ_x, ϕ_gas) rtol = 1.0e-10
        nd32 = convert_precision(Float32, nd)
        @test nd32.crystals isa Roscoe1952_Crystals{Float32}
        τ32 = compute_τII(nd32, Float32(ε_nd); ϕ_x = Float32(ϕ_x), ϕ_gas = Float32(ϕ_gas))
        @test τ32 isa Float32
        @test τ32 ≈ τ_nd rtol = 1.0e-5
    end

    @testset "capillary-number dependent bubbles" begin
        rv = compute_relative_viscosity
        LM = Mader2013_Bubbles()
        slope(f, ϕ, Ca) = GeoParams.MaterialParameters.ConstitutiveRelationships.relative_viscosity_and_slope(f, ϕ, Ca)[2]
        # dilute limit: Llewellin et al. (2002, eq. 11) with b = 1
        dilute(ϕ, Ca) = 1 + ϕ * (1 - 12 / 5 * Ca^2) / (1 + 36 / 25 * Ca^2)

        @test LM.r.val == 1.0e-4
        @test LM.Γ.val == 0.3
        for Ca in (0.0, 0.3, 1.0, 10.0)
            @test rv(LM, 1.0e-6, Ca) ≈ dilute(1.0e-6, Ca) rtol = 1.0e-10
        end
        let f = LM
            @test rv(f, 0.5, 0) ≈ 2 rtol = 1.0e-14
            @test rv(f, 0.5, 1.0e8) ≈ 0.5^(5 / 3) rtol = 1.0e-10
            @test rv(f, 0, 3.0) ≈ 1 rtol = 1.0e-14
            @test isnan(rv(f, -0.1, 1.0))
            @test isnan(rv(f, 1.1, 1.0))
            @test_throws "depends on the strain rate" rv(f, 0.5)
            # d ln η_r / d ln Ca against central differences
            for Ca in (0.05, 0.4, 0.9, 5.0), ϕ in (0.2, 0.6)
                h = 1.0e-6
                fd = (log(rv(f, ϕ, Ca * exp(h))) - log(rv(f, ϕ, Ca * exp(-h)))) / 2h
                @test slope(f, ϕ, Ca) ≈ fd rtol = 1.0e-6
            end
        end
        η_m, ϕ_xs, ϕ_g = 1.0e5, 0.2, 0.4
        melt = LinearViscous(η = η_m * Pa * s)
        r, Γ = 1.0e-4, 0.3
        let f = LM
            law = Suspension(melt; crystals = ER, bubbles = f)
            a = (; ϕ_x = ϕ_xs, ϕ_gas = ϕ_g)
            ηrx = rv(ER, ϕ_xs)
            ϕ_b = ϕ_g / (1 - ϕ_xs)
            @test isbits(law)
            for ε in (1.0e-5, 1.0e-2, 1.0, 1.0e2)
                # bubbly melt sheared at √(η_r,x/(1 − ϕ_x)) ε̇ (Chateau et al. 2008): Ca = r τ / Γ
                Ca = r * 2 * η_m * sqrt(ηrx / (1 - ϕ_xs)) * ε / Γ
                τ_ε = compute_τII(law, ε; a...)
                @test τ_ε ≈ 2 * η_m * ηrx * rv(f, ϕ_b, Ca) * ε rtol = 1.0e-12
                @test compute_εII(law, τ_ε; a...) ≈ ε rtol = 1.0e-10
                h = 1.0e-6 * ε
                fd = (compute_τII(law, ε + h; a...) - compute_τII(law, ε - h; a...)) / 2h
                @test dτII_dεII(law, ε; a...) ≈ fd rtol = 1.0e-6
                @test dεII_dτII(law, τ_ε; a...) * dτII_dεII(law, ε; a...) ≈ 1 rtol = 1.0e-10
            end
            @test compute_τII(law, 0.0; a...) == 0
            @test compute_εII(law, 0.0; a...) == 0
            @test dτII_dεII(law, 0.0; a...) ≈ 2 * η_m * ηrx / (1 - ϕ_b) rtol = 1.0e-12
            @test isnan(compute_εII(law, τ; ϕ_x = ϕ_xs, ϕ_gas = 1.2))
            @test_throws "depends on the strain rate" rv(law; a...)

            # nonlinear melt
            dis = Suspension(SetDislocationCreep(GeoParams.Dislocation.wet_olivine_Mei_2000b); crystals = ER, bubbles = f)
            ad = (; T = 1473.0, P = 300.0e6, f = 300.0e6, a...)
            @test compute_τII(dis, compute_εII(dis, τ; ad...); ad...) ≈ τ rtol = 1.0e-10

            # arrays and composite rheology
            εs = [1.0e-3, 1.0e-1, 10.0]
            S = zeros(3)
            compute_τII!(S, law, εs; ϕ_x = [0.0, 0.2, 0.4], ϕ_gas = ϕ_g)
            @test S[2] ≈ compute_τII(law, εs[2]; a...) rtol = 1.0e-14
            E = zeros(3)
            compute_εII!(E, law, S; ϕ_x = [0.0, 0.2, 0.4], ϕ_gas = ϕ_g)
            @test E ≈ εs rtol = 1.0e-10
            c = CompositeRheology(law)
            @test compute_τII(c, 1.0, a) ≈ compute_τII(law, 1.0; a...) rtol = 1.0e-8
            @test compute_εII(c, τ, a) ≈ compute_εII(law, τ; a...) rtol = 1.0e-8

            # nondimensional and Float32
            CD = GEO_units()
            nd = nondimensionalize(law, CD)
            ε_nd = nondimensionalize(1.0 / s, CD)
            τ_nd = compute_τII(nd, ε_nd; a...)
            @test ustrip(dimensionalize(τ_nd, Pa, CD)) ≈ compute_τII(law, 1.0; a...) rtol = 1.0e-10
            nd32 = convert_precision(Float32, nd)
            a32 = (; ϕ_x = Float32(ϕ_xs), ϕ_gas = Float32(ϕ_g))
            τ32 = compute_τII(nd32, Float32(ε_nd); a32...)
            @test τ32 isa Float32
            @test τ32 ≈ τ_nd rtol = 1.0e-5
            ε32 = compute_εII(nd32, τ32; a32...)
            @test ε32 isa Float32
            @test ε32 ≈ ε_nd rtol = 1.0e-4
        end
    end

    @testset "Truby et al. (2015) three-phase suspension" begin
        # sample 3P-14 (Table 1): ϕ_p = 0.272, whole-magma gas fraction ϕ*_b = 0.098, η_* = 235 Pa s, μ = 61.7 Pa s,
        # with the Maron–Pierce fit ϕ_m = 0.593 to their bubble-free suspensions
        MP = Suspension(LinearViscous(η = 61.7Pa * s); crystals = Roscoe1952_Crystals(ϕ_max = 0.593, n = 2), bubbles = IZ)
        @test compute_relative_viscosity(MP; ϕ_x = 0.272, ϕ_gas = 0.098) ≈ 235 / 61.7 rtol = 0.05
    end

    @testset "Mueller et al. (2010) maximum packing" begin
        f = Mueller2010_Crystals(aspect_ratio = 7.25)
        @test f isa Roscoe1952_Crystals
        @test f.n == 2
        # ϕ_m that Vona et al. (2011, Table 1) compute with eq. 5.1
        @test f.ϕ_max ≈ 0.3741 rtol = 1.0e-3
        @test Mueller2010_Crystals(aspect_ratio = 5.57).ϕ_max ≈ 0.416 rtol = 1.0e-3
        # Maron–Pierce fits of Mueller et al. (2010, Table 3), to which eq. 5.1 is a regression
        for (rp, ϕm) in ((1, 0.633), (1.75, 0.539), (4.69, 0.433), (9.17, 0.339))
            @test Mueller2010_Crystals(aspect_ratio = rp).ϕ_max ≈ ϕm rtol = 0.06
        end
    end

    @testset "Vona et al. (2011) shear-thinning crystals" begin
        rv = compute_relative_viscosity
        V = Vona2011_Crystals(aspect_ratio = 7.25)
        slope(f, ϕ, x) = GeoParams.MaterialParameters.ConstitutiveRelationships.relative_viscosity_and_slope(f, ϕ, x)[2]
        @test V.ϕ_max == Mueller2010_Crystals(aspect_ratio = 7.25).ϕ_max
        @test V.a == 0.06
        # sample SGP6 (Table 1): ϕ = 1 - 0.728, apparent viscosity at shear rates 0.02 and 0.14 1/s, melt 10^2.54 Pa s
        @test rv(V, 0.272, 0.02) ≈ 6058 / 10^2.54 rtol = 0.01
        @test rv(V, 0.272, 0.14) ≈ 5180 / 10^2.54 rtol = 0.03
        # consistency at 1/s is Maron–Pierce, flow index n = 1 + 2a log10(1 - ϕ/ϕ_m) (eq. 20)
        @test rv(V, 0.2, 1.0) ≈ rv(Mueller2010_Crystals(aspect_ratio = 7.25), 0.2) rtol = 1.0e-14
        @test slope(V, 0.2, 3.0) ≈ 0.12 * log10(1 - 0.2 / V.ϕ_max) rtol = 1.0e-14
        h = 1.0e-6
        @test slope(V, 0.2, 3.0) ≈ (log(rv(V, 0.2, 3exp(h))) - log(rv(V, 0.2, 3exp(-h)))) / 2h rtol = 1.0e-8
        @test rv(V, 0, 0.1) == 1
        @test rv(V, V.ϕ_max, 0.1) == Inf
        @test isnan(rv(V, -0.1, 1.0))
        @test isnan(rv(V, 0.4, 1.0))
        @test_throws "depends on the strain rate" rv(V, 0.2)

        melt = GiordanoMeltViscosity()
        law = Suspension(melt; crystals = V)
        both = Suspension(melt; crystals = V, bubbles = Mader2013_Bubbles())
        a = (; T = 1400.0, ϕ_x = 0.2, ϕ_gas = 0.1)
        @test isbits(law)
        @test_throws "depends on the strain rate" rv(law; a...)
        @test_throws "rate-dependent bubbles must depend on the capillary number" Suspension(melt; bubbles = V)
        for ε in (1.0e-6, 1.0e-2, 1.0)
            # η_r,x is evaluated at the shear rate 2ε̇ of the suspension
            τ_ε = compute_τII(law, ε; a...)
            @test τ_ε ≈ compute_τII(melt, rv(V, 0.2, 2ε) * ε; a...) rtol = 1.0e-12
            for l in (law, both)
                τ_l = compute_τII(l, ε; a...)
                @test compute_εII(l, τ_l; a...) ≈ ε rtol = 1.0e-10
                h = 1.0e-6 * ε
                fd = (compute_τII(l, ε + h; a...) - compute_τII(l, ε - h; a...)) / 2h
                @test dτII_dεII(l, ε; a...) ≈ fd rtol = 1.0e-6
                @test dεII_dτII(l, τ_l; a...) * dτII_dεII(l, ε; a...) ≈ 1 rtol = 1.0e-10
            end
        end
        @test compute_τII(law, 0.0; a...) == 0
        @test compute_εII(law, 0.0; a...) == 0

        # nondimensional and Float32
        lin = Suspension(LinearViscous(η = 1.0e5Pa * s); crystals = V, bubbles = Mader2013_Bubbles())
        al = (; ϕ_x = 0.2, ϕ_gas = 0.1)
        CD = GEO_units()
        nd = nondimensionalize(lin, CD)
        ε_nd = nondimensionalize(0.1 / s, CD)
        τ_nd = compute_τII(nd, ε_nd; al...)
        @test ustrip(dimensionalize(τ_nd, Pa, CD)) ≈ compute_τII(lin, 0.1; al...) rtol = 1.0e-10
        nd32 = convert_precision(Float32, nd)
        @test nd32.crystals isa Vona2011_Crystals{Float32}
        a32 = (; ϕ_x = 0.2f0, ϕ_gas = 0.1f0)
        τ32 = compute_τII(nd32, Float32(ε_nd); a32...)
        @test τ32 isa Float32
        @test τ32 ≈ τ_nd rtol = 1.0e-5
        ε32 = compute_εII(nd32, τ32; a32...)
        @test ε32 isa Float32
        @test ε32 ≈ ε_nd rtol = 1.0e-4
    end

    @testset "negative stress, rate-dependent inverse" begin
        sus = Suspension(LinearViscous(η = 1.0e5Pa * s); crystals = Vona2011_Crystals(), bubbles = Mader2013_Bubbles())
        a = (; ϕ_x = 0.2, ϕ_gas = 0.3)
        @test compute_εII(sus, -1.0e4; a...) ≈ -compute_εII(sus, 1.0e4; a...) rtol = 1.0e-12
        # τ(ε̇) is odd, so compute_τII and compute_εII stay inverses for negative strain rates
        for ε in (1.0e-3, 1.0)
            @test compute_τII(sus, -ε; a...) ≈ -compute_τII(sus, ε; a...) rtol = 1.0e-12
            @test compute_εII(sus, compute_τII(sus, -ε; a...); a...) ≈ -ε rtol = 1.0e-10
        end
    end

    @testset "Vona et al. (2011) at zero shear rate" begin
        V = Vona2011_Crystals()
        rv = compute_relative_viscosity
        # η_r is held below |γ̇| = 1e-30 γ̇0, with the same value in Float32 and Float64
        @test rv(V, 0.2, 0.0) == rv(V, 0.2, 1.0e-40) == rv(V, 0.2, 1.0e-30)
        @test rv(V, 0.2f0, 0.0f0) ≈ rv(V, 0.2, 0.0) rtol = 1.0e-5
        sus = Suspension(LinearViscous(η = 1.0e5Pa * s); crystals = V)
        a = (; ϕ_x = 0.2, ϕ_gas = 0.0)
        @test dτII_dεII(sus, 0.0; a...) ≈ 2 * 1.0e5 * rv(V, 0.2, 0.0) rtol = 1.0e-12
    end

    @testset "invalid parameters" begin
        @test_throws "aspect_ratio must be at least 1" Mueller2010_Crystals(aspect_ratio = 0.5)
        @test_throws "aspect_ratio must be at least 1" Vona2011_Crystals(aspect_ratio = 0.5)
        @test_throws "a must be non-negative" Vona2011_Crystals(a = -0.1)
        @test_throws "ϕ_max must lie in (0, 1]" Roscoe1952_Crystals(ϕ_max = 0)
        @test_throws "ϕ_max must lie in (0, 1]" Roscoe1952_Crystals(ϕ_max = 1.2)
        @test_throws "n must be positive" Roscoe1952_Crystals(n = 0)
        @test_throws "must be relative viscosities or nothing" Suspension(melt; crystals = 0.6)
        @test_throws "crystals cannot depend on the capillary number" Suspension(melt; crystals = Mader2013_Bubbles())
        # fractions are checked even without crystals or bubbles
        @test isnan(compute_relative_viscosity(Suspension(melt; crystals = Roscoe1952_Crystals()); ϕ_x = 0.5, ϕ_gas = 0.9))
        @test isnan(compute_relative_viscosity(Suspension(melt; bubbles = IshiiZuber1979_Bubbles()); ϕ_x = -0.1, ϕ_gas = 0.1))
        @test_throws "bubble radius r must be positive" Mader2013_Bubbles(r = 0m)
        @test_throws "surface tension Γ must be positive" Mader2013_Bubbles(Γ = -1J / m^2)
    end
end
