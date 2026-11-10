using Test
using GeoParams

@testset "Thermal properties of bubbly magma" begin
    CharDim = GEO_units()

    @testset "Moitra et al. (2018) glass/melt" begin
        # Table 3: T [K], cp [J/kg/K], k [W/m/K]; the 973 K point lies off the fits (glass transition).
        data = ((298.0, 812.0, 1.53), (473.0, 1010.0, 1.6), (673.0, 1130.0, 1.67), (1273.0, 1510.0, 2.37))
        cp, k = T_HeatCapacity_Moitra(), T_Conductivity_Moitra()
        for (T, cp_data, k_data) in data
            @test compute_heatcapacity(cp, (; T)) ≈ cp_data rtol = 0.01
            @test compute_conductivity(k, (; T)) ≈ k_data rtol = 0.05
        end
        # eq. 6 at 1273 K, with ρ = 3000 kg/m³ and eq. 4
        @test compute_conductivity(k, (; T = 1273.0)) ≈ 3000 * (2150 - 0.5 * 1273) * (-0.06 + 4.57e-4 * 1273) * 1.0e-6
        @test ustrip(compute_conductivity(k, (; T = 900.0K))) ≈ compute_conductivity(k, (; T = 900.0))
        @test isbits(k)

        Tn = nondimensionalize(900.0K, CharDim)
        kn = nondimensionalize(k, CharDim)
        @test ustrip(dimensionalize(compute_conductivity(kn, (; T = Tn)), Watt / m / K, CharDim)) ≈ compute_conductivity(k, (; T = 900.0))
        cpn = nondimensionalize(cp, CharDim)
        @test ustrip(dimensionalize(compute_heatcapacity(cpn, (; T = Tn)), J / kg / K, CharDim)) ≈ compute_heatcapacity(cp, (; T = 900.0))

        @test compute_conductivity(k, (; T = 900.0f0)) isa Float32
        @test compute_heatcapacity(cp, (; T = 900.0f0)) isa Float32
    end

    @testset "IAPWS 2011 H2O conductivity" begin
        # IAPWS R15-11, Table 4 (λ̄2 = 0 there): T [K], ρ [kg/m³], λ [mW/m/K]
        for (T, ρ, λ) in ((298.15, 998.0, 607.712868), (298.15, 1200.0, 799.038144), (873.15, 0.0, 79.1034659))
            law = IAPWS_H2O_Conductivity(; ρ = ConstantDensity(ρ = ρ * kg / m^3))
            @test compute_conductivity(law, (; P = 1.0e8, T)) ≈ 1.0e-3λ rtol = 1.0e-8
        end
        # λ0·λ1 with the IAPWS-95 density at 850 °C, 200 MPa (iapws: IAPWS95(T, P).rho, _ThCond(rho, T))
        law = IAPWS_H2O_Conductivity(; ρ = ConstantDensity(ρ = 397.66477466060775kg / m^3))
        @test compute_conductivity(law, (; P = 2.0e8, T = 1123.15)) ≈ 0.3664589997022277 rtol = 1.0e-6

        # Full IAPWS 2011 conductivity (iapws: IAPWS95(T, P).k): T [K], P [Pa], k [W/m/K]
        k = IAPWS_H2O_Conductivity()
        data = (
            (873.15, 1.0e5, 0.07917), (873.15, 1.0e8, 0.2945), (873.15, 4.0e8, 0.6862),
            (1123.15, 1.0e5, 0.1127), (1123.15, 3.0e7, 0.1407), (1123.15, 2.0e8, 0.3665), (1123.15, 4.0e8, 0.5573),
            (1473.15, 1.0e5, 0.1626), (1473.15, 1.0e8, 0.2744), (1473.15, 4.0e8, 0.5183),
        )
        for (T, P, k_ref) in data
            @test compute_conductivity(k, (; P, T)) ≈ k_ref rtol = 0.016
        end
        # Dilute-gas limit λ*·λ̄0 (Table 4, 873.15 K)
        @test compute_conductivity(k, (; P = 1.0, T = 873.15)) ≈ 0.0791034659 rtol = 1.0e-6

        args = (; P = 2.0e8, T = 1123.15)
        @test ustrip(compute_conductivity(k, (; P = 2.0e8Pa, T = 1123.15K))) ≈ compute_conductivity(k, args)
        kn = nondimensionalize(k, CharDim)
        args_nd = (; P = nondimensionalize(2.0e8Pa, CharDim), T = nondimensionalize(1123.15K, CharDim))
        @test ustrip(dimensionalize(compute_conductivity(kn, args_nd), Watt / m / K, CharDim)) ≈ compute_conductivity(k, args)
        @test compute_conductivity(k, (; P = 2.0f8, T = 1123.15f0)) isa Float32
        @test isbits(k)
        allocs(law, args) = (compute_conductivity(law, args); @allocated compute_conductivity(law, args))
        @test allocs(k, args) == 0
    end

    @testset "Udoetok (2013) gas mixture" begin
        # Table 2, N2-He at 30 °C (k in 1e-5 cal/°C/cm/s); He is the second gas here.
        kmix = GasMixture_Conductivity(; kh2o = ConstantConductivity(k = 6.116), kco2 = ConstantConductivity(k = 36.35))
        @test compute_conductivity(kmix, (; X_co2 = 0.4681)) ≈ 15.142 rtol = 1.0e-4
        @test compute_conductivity(kmix, (; X_co2 = 0.2893)) ≈ 11.4583 rtol = 1.0e-4
        @test compute_conductivity(kmix, (; X_co2 = 0)) ≈ 6.116
        @test compute_conductivity(kmix, (; X_co2 = 1)) ≈ 36.35
        @test isnan(compute_conductivity(kmix, (; X_co2 = 1.1)))
        @test isnan(compute_conductivity(kmix, (; X_co2 = -0.1)))
    end

    @testset "Heap et al. (2020) Maxwell conductivity" begin
        # Water-saturated Ruapehu andesite R17 (Table 3): ϕ = 0.615, λ = 0.90 W/m/K; λ0 = 1.5, λ_water = 0.6.
        wet = Maxwell_Conductivity(; kmatrix = ConstantConductivity(k = 1.5Watt / m / K), kgas = ConstantConductivity(k = 0.6Watt / m / K))
        @test compute_conductivity(wet, (; ϕ_gas = 0.615)) ≈ 0.9 atol = 0.01
        # Dry samples (λ_air = 0): λ = λ0 · 2(1-ϕ)/(2+ϕ)
        dry = Maxwell_Conductivity(; kmatrix = ConstantConductivity(k = 1.5Watt / m / K), kgas = ConstantConductivity(k = 0.0Watt / m / K))
        @test compute_conductivity(dry, (; ϕ_gas = 0.602)) ≈ 1.5 * 2 * (1 - 0.602) / (2 + 0.602)
        @test compute_conductivity(wet, (; ϕ_gas = 0)) == 1.5
        @test compute_conductivity(wet, (; ϕ_gas = 1)) ≈ 0.6
        # m = 1 is finite
        same = Maxwell_Conductivity(; kmatrix = ConstantConductivity(), kgas = ConstantConductivity())
        @test compute_conductivity(same, (; ϕ_gas = 0.3)) ≈ 3.0
        @test isnan(compute_conductivity(wet, (; ϕ_gas = 1.1)))
        @test isnan(compute_conductivity(wet, (; ϕ_gas = -0.1)))

        # Defaults: Moitra matrix + Udoetok gas, dimensional vs nondimensional
        k = Maxwell_Conductivity()
        @test isbits(k)
        args = (; T = 1200.0, P = 1.0e8, ϕ_gas = 0.3, X_co2 = 0.2)
        k_dim = compute_conductivity(k, args)
        @test k_dim < compute_conductivity(k.kmatrix, args)
        kn = nondimensionalize(k, CharDim)
        args_nd = (; T = nondimensionalize(1200.0K, CharDim), P = nondimensionalize(1.0e8Pa, CharDim), ϕ_gas = 0.3, X_co2 = 0.2)
        @test ustrip(dimensionalize(compute_conductivity(kn, args_nd), Watt / m / K, CharDim)) ≈ k_dim
        @test compute_conductivity(k, (; T = 1200.0f0, P = 1.0f8, ϕ_gas = 0.3f0, X_co2 = 0.2f0)) isa Float32
    end

    @testset "Heap et al. (2020) mass-weighted heat capacity" begin
        # Water-saturated Ruapehu andesite R17 (Table 3): ϕ = 0.615, cp = 2.02 kJ/kg/K, with the paper's
        # ρs = 2750 kg/m³, cp = 750 J/kg/K and water ρf = 1000 kg/m³, cp,f = 4182 J/kg/K.
        ρ = ThreePhase_Density(; ρmelt = ConstantDensity(ρ = 2750kg / m^3), ρgas = ConstantDensity(ρ = 1000kg / m^3))
        wet = ThreePhase_HeatCapacity(; ρ, Cpmelt = ConstantHeatCapacity(Cp = 750J / kg / K), Cpgas = ConstantHeatCapacity(Cp = 4182J / kg / K))
        @test compute_heatcapacity(wet, (; T = 300.0, ϕ_gas = 0.615)) ≈ 2020 rtol = 0.01
        @test compute_heatcapacity(wet, (; T = 300.0, ϕ_gas = 0)) == 750
        @test compute_heatcapacity(wet; T = 300.0, ϕ_gas = 1) ≈ 4182

        # ρ·cp is the volume-weighted sum of the phase values
        Cp = ThreePhase_HeatCapacity()
        args = (; T = 1200.0, P = 1.0e8, ϕ_x = 0.2, ϕ_gas = 0.3)
        ρmix = compute_density(Cp.ρ, args)
        ρcp = sum(
            ϕ * compute_density(ρi, args) * compute_heatcapacity(cpi, args)
                for (ϕ, ρi, cpi) in ((0.5, Cp.ρ.ρmelt, Cp.Cpmelt), (0.2, Cp.ρ.ρsolid, Cp.Cpsolid), (0.3, Cp.ρ.ρgas, Cp.Cpgas))
        )
        @test compute_heatcapacity(Cp, args) ≈ ρcp / ρmix
        # ϕ_gas = 0 never evaluates the gas density (IdealGas_Density at P = 0 would give 0/0 weights otherwise)
        @test compute_heatcapacity(Cp, (; T = 1200.0, ϕ_x = 0.0)) == compute_heatcapacity(Cp.Cpmelt, (; T = 1200.0))
        @test isbits(Cp)
        for bad in ((; ϕ_x = -0.1), (; ϕ_gas = -0.1), (; ϕ_x = 0.6, ϕ_gas = 0.5))
            @test isnan(compute_heatcapacity(Cp, (; T = 1200.0, P = 1.0e8, bad...)))
        end

        Cpn = nondimensionalize(Cp, CharDim)
        args_nd = (; T = nondimensionalize(1200.0K, CharDim), P = nondimensionalize(1.0e8Pa, CharDim), ϕ_x = 0.2, ϕ_gas = 0.3)
        @test ustrip(dimensionalize(compute_heatcapacity(Cpn, args_nd), J / kg / K, CharDim)) ≈ compute_heatcapacity(Cp, args)
        @test compute_heatcapacity(Cp, (; T = 1200.0f0, P = 1.0f8, ϕ_x = 0.2f0, ϕ_gas = 0.3f0)) isa Float32
    end

    @testset "Moitra et al. (2018), smoothed at Tg" begin
        cp, k = T_HeatCapacity_Moitra_smoothed(), T_Conductivity_Moitra_smoothed()
        # Away from Tg = 1000 K (≥ 5ΔTg) the smoothed laws equal the stepped ones.
        for T in (298.0, 473.0, 500.0, 900.0, 1100.0, 1500.0)
            @test compute_heatcapacity(cp, (; T)) ≈ compute_heatcapacity(T_HeatCapacity_Moitra(), (; T)) rtol = 1.0e-6
            @test compute_conductivity(k, (; T)) ≈ compute_conductivity(T_Conductivity_Moitra(), (; T)) rtol = 1.0e-6
        end
        # Across Tg (±2.5ΔTg) the value moves monotonically from the glass branch to the melt branch.
        glass, melt = T_HeatCapacity_Moitra(; Tg = 1.0e4K), T_HeatCapacity_Moitra(; Tg = 0.0K)
        weight(T) = (compute_heatcapacity(cp, (; T)) - compute_heatcapacity(glass, (; T))) / (compute_heatcapacity(melt, (; T)) - compute_heatcapacity(glass, (; T)))
        Ts = range(950.0, 1050.0; length = 1001)
        w = weight.(Ts)
        @test all(diff(w) .> 0) && 0 < first(w) < 0.01 && 0.99 < last(w) < 1
        @test weight(1000.0) ≈ 0.5
        kv = [compute_conductivity(k, (; T)) for T in Ts]
        @test maximum(abs, diff(kv)) < 0.01    # no jump at 0.1 K spacing (the stepped law jumps by 0.33 W/m/K)

        @test ustrip(compute_conductivity(k, (; T = 1000.0K))) ≈ compute_conductivity(k, (; T = 1000.0))
        Tn = nondimensionalize(1000.0K, CharDim)
        kn, cpn = nondimensionalize(k, CharDim), nondimensionalize(cp, CharDim)
        @test ustrip(dimensionalize(compute_conductivity(kn, (; T = Tn)), Watt / m / K, CharDim)) ≈ compute_conductivity(k, (; T = 1000.0))
        @test ustrip(dimensionalize(compute_heatcapacity(cpn, (; T = Tn)), J / kg / K, CharDim)) ≈ compute_heatcapacity(cp, (; T = 1000.0))
        @test compute_conductivity(k, (; T = 1000.0f0)) isa Float32
        @test compute_heatcapacity(cp, (; T = 1000.0f0)) isa Float32
        @test isbits(k) && isbits(cp)
    end

    @testset "Unitful arguments" begin
        args = (; T = 1200.0, P = 1.0e8, X_co2 = 0.2, ϕ_x = 0.2, ϕ_gas = 0.3)
        args_u = (; args..., T = 1200.0K, P = 1.0e8Pa)
        for law in (ConstantConductivity(), GasMixture_Conductivity(), Maxwell_Conductivity(), T_Conductivity_Moitra_smoothed())
            @test ustrip(uconvert(Watt / m / K, compute_conductivity(law, args_u))) ≈ compute_conductivity(law, args)
        end
        for law in (ConstantHeatCapacity(), ThreePhase_HeatCapacity())
            @test ustrip(uconvert(J / kg / K, compute_heatcapacity(law, args_u))) ≈ compute_heatcapacity(law, args)
        end
    end

    @testset "show and param_info" begin
        for law in (T_HeatCapacity_Moitra(), T_HeatCapacity_Moitra_smoothed(), ThreePhase_HeatCapacity(), T_Conductivity_Moitra(), T_Conductivity_Moitra_smoothed(), IAPWS_H2O_Conductivity(), GasMixture_Conductivity(), Maxwell_Conductivity())
            @test param_info(law) isa MaterialParamsInfo
            @test !isempty(sprint(show, law))
        end
    end
end
