using Test, GeoParams, StaticArrays, LaTeXStrings
using ForwardDiff: derivative
const No_MaterialParam = GeoParams.MaterialParameters.No_MaterialParam

@testset "Solubility.jl" begin

    # Liu et al. (2005) silicic -------------------------------------------
    s = Liu2005_Solubility()
    @test isbits(s)
    @test isdimensional(s) === true
    @test param_info(s).Equation isa LaTeXString
    @test occursin("Liu", sprint(show, s))

    # reference value: ~5.7 wt% H2O, no CO2 at X_co2 = 0 (200 MPa, 1200 K)
    mh, mc = compute_dissolved(s, 2.0e8, 1200.0, 0.0)
    @test 0.03 < mh < 0.07
    @test mc == 0.0
    # CO2 dissolves when present in the gas
    @test compute_dissolved(s, 2.0e8, 1200.0, 0.3)[2] > 0.0
    # more soluble at higher pressure
    @test compute_dissolved(s, 4.0e8, 1200.0, 0.0)[1] > mh

    # keyword and NamedTuple call forms
    @test compute_dissolved(s; P = 2.0e8, T = 1200.0, X_co2 = 0.0) == (mh, mc)
    @test compute_dissolved(s, (; P = 2.0e8, T = 1200.0, X_co2 = 0.0)) == (mh, mc)

    # dimensional Quantity input (the @unpack_units branch)
    mhq, mcq = compute_dissolved(s, 200.0MPa, 1200.0K, 0.0)
    @test ustrip(mhq) ≈ mh rtol = 1.0e-6

    # nondimensionalization: mass fractions are dimensionless -> invariant
    CharDim = GEO_units(; viscosity = 1.0e19, length = 1000km)
    s_nd = nondimensionalize(s, CharDim)
    @test isdimensional(s_nd) === false
    mh_nd, mc_nd = compute_dissolved(
        s_nd,
        nondimensionalize(2.0e8Pa, CharDim), nondimensionalize(1200.0K, CharDim), 0.0
    )
    @test mh_nd ≈ mh
    @test mc_nd == 0.0

    # ForwardDiff partials vs central finite differences
    P0, T0, X0 = 2.0e8, 1200.0, 0.3
    dP = ∂dissolved_∂P(s, P0, T0, X0)
    dT = ∂dissolved_∂T(s, P0, T0, X0)
    dX = ∂dissolved_∂Xco2(s, P0, T0, X0)
    @test length(dP) == 2 && length(dT) == 2 && length(dX) == 2
    fd(f, x, h) = (f(x + h) .- f(x - h)) ./ (2h)
    fdP = fd(p -> SVector(compute_dissolved(s, p, T0, X0)...), P0, 1.0e2)
    fdT = fd(t -> SVector(compute_dissolved(s, P0, t, X0)...), T0, 1.0e-3)
    fdX = fd(x -> SVector(compute_dissolved(s, P0, T0, x)...), X0, 1.0e-6)
    @test all(isapprox.(dP, Tuple(fdP); rtol = 1.0e-4))
    @test all(isapprox.(dT, Tuple(fdT); rtol = 1.0e-4))
    @test all(isapprox.(dX, Tuple(fdX); rtol = 1.0e-4))

    # degenerate X_co2 = 0: CO2 output and its pressure sensitivity vanish
    @test ∂dissolved_∂P(s, P0, T0, 0.0)[2] == 0.0

    # find_Xco2: invert compute_dissolved for X_co2 given a target m_h2o
    m_target = compute_dissolved(s, P0, T0, X0)[1]
    @test find_Xco2(s, P0, T0, m_target) ≈ X0 atol = 1.0e-6
    @test find_Xco2(s, P0, T0, compute_dissolved(s, P0, T0, 0.0)[1]) ≈ 0.0 atol = 1.0e-6
    @test find_Xco2(s, P0, T0, compute_dissolved(s, P0, T0, 1.0)[1]) ≈ 1.0 atol = 1.0e-6
    @test_throws ErrorException find_Xco2(s, P0, T0, 10.0) # exceeds achievable range

    # Mafic ----------------------------------------------------------------
    m = Mafic_Solubility()
    @test isbits(m)
    @test isdimensional(m) === true
    @test param_info(m).Equation isa LaTeXString
    @test occursin("Mafic", sprint(show, m))

    mmh, mmc = compute_dissolved(m, 2.0e8, 1200.0, 0.0)
    @test mmh > 0.0
    @test mmc == 0.0
    @test compute_dissolved(m, 2.0e8, 1200.0, 0.3)[2] > 0.0
    # Quantity branch
    @test ustrip(compute_dissolved(m, 200.0MPa, 1200.0K, 0.0)[1]) ≈ mmh rtol = 1.0e-6
    # nondim invariance
    m_nd = nondimensionalize(m, CharDim)
    @test isdimensional(m_nd) === false
    @test compute_dissolved(
        m_nd,
        nondimensionalize(2.0e8Pa, CharDim), nondimensionalize(1200.0K, CharDim), 0.0
    )[1] ≈ mmh
    # derivatives exist and are finite
    @test all(isfinite, ∂dissolved_∂T(m, P0, T0, X0))

    # find_Xco2 for the mafic law too (generic over AbstractSolubility)
    mm_target = compute_dissolved(m, P0, T0, X0)[1]
    @test find_Xco2(m, P0, T0, mm_target) ≈ X0 atol = 1.0e-6
    @test_throws DomainError find_Xco2(m, P0, T0, -0.005)   # not a mass fraction

    # the mafic H2O polynomial has no floor at zero and goes negative outside
    # its calibration; dissolved H2O cannot be negative, so it is floored at
    # zero rather than handed back as an unphysical mass fraction
    @test compute_dissolved(m, 3.0e7, 1000.0, 0.9)[1] == 0.0

    # Sun & Yao (2026) -----------------------------------------------------
    sy = SunYao2026_Solubility()
    @test isbits(sy)
    @test isdimensional(sy) === true
    @test param_info(sy).Equation isa LaTeXString
    @test occursin("Sun & Yao", sprint(show, sy))
    @test sum(sy.X) ≈ 1

    # eq. 7 by hand: X_Na2O = 0.04, X_Al2O3 = 0.087, 2000 bar, 1100 K, pure H2O gas
    # gives 4.0 wt% molecular H2O + 2.0 wt% hydroxyl
    sy_hand = SunYao2026_Solubility((0, 0, 0.087, 0, 0, 0, 0.04, 0, 0), sy.a, sy.b, 1.0e5Pa, 1.0K)
    @test compute_dissolved(sy_hand, 2.0e8, 1100.0, 0.0)[1] ≈ 0.06 rtol = 0.01

    # rhyolite values reported by Sun & Yao (2026) for a similar rhyolite:
    # ~16-17 wt% H2O at 10 kbar and 900 °C (Fig. 5a), reduced by up to ~20 % at 1200 °C (sect. 4.1)
    mh_900 = compute_dissolved(sy, 1.0e9, 1173.15, 0.0)[1]
    @test 0.15 < mh_900 < 0.18
    @test 0.8 < compute_dissolved(sy, 1.0e9, 1473.15, 0.0)[1] / mh_900 < 1
    # pure CO2 gas: < 10 ppm at 10 bar, ~0.4 wt% at 10 kbar (sect. 4.1, Fig. 5c)
    @test compute_dissolved(sy, 1.0e6, 1173.15, 1.0)[2] < 1.0e-5
    @test 0.003 < compute_dissolved(sy, 1.0e9, 1173.15, 1.0)[2] < 0.006
    # 4 wt% H2O + 0.1 wt% CO2 saturates at ~3 kbar with a CO2-rich gas (sect. 4.2)
    let X = find_Xco2(sy, 3.0e8, 1173.15, 0.04)
        @test 0.5 < X < 0.7
        @test compute_dissolved(sy, 3.0e8, 1173.15, X)[2] ≈ 0.001 rtol = 0.2
    end

    # pure H2O gas dissolves no CO2, pure CO2 gas no H2O
    @test compute_dissolved(sy, 2.0e8, 1173.15, 0.0)[2] == 0.0
    @test compute_dissolved(sy, 2.0e8, 1173.15, 1.0)[1] == 0.0
    # and their P and T derivatives are finite, with exact zeros for the absent species
    for ∂ in (∂dissolved_∂P, ∂dissolved_∂T)
        let d0 = ∂(sy, 2.0e8, 1173.15, 0.0), d1 = ∂(sy, 2.0e8, 1173.15, 1.0)
            @test isfinite(d0[1]) && d0[2] == 0.0
            @test d1[1] == 0.0
            @test d1[2] ≈ ∂(sy, 2.0e8, 1173.15, prevfloat(1.0))[2] rtol = 1.0e-6
        end
    end
    # negative partial pressures are outside the domain
    @test_throws DomainError compute_dissolved(sy, 2.0e8, 1173.15, 1.1)
    @test_throws DomainError compute_dissolved(sy, 2.0e8, 1173.15, -0.1)

    sy_ref = compute_dissolved(sy, 2.0e8, 1173.15, 0.3)
    @test all(ustrip.(compute_dissolved(sy, 200.0MPa, 1173.15K, 0.3)) .≈ sy_ref)
    sy_nd = nondimensionalize(sy, CharDim)
    @test isdimensional(sy_nd) === false
    @test all(compute_dissolved(sy_nd, nondimensionalize(2.0e8Pa, CharDim), nondimensionalize(1173.15K, CharDim), 0.3) .≈ sy_ref)
    sy32 = convert_precision(Float32, sy)
    @test @inferred(compute_dissolved(sy32, 2.0f8, 1173.15f0, 0.3f0)) isa Tuple{Float32, Float32}
    @test all(isapprox.(compute_dissolved(sy32, 2.0f8, 1173.15f0, 0.3f0), sy_ref; rtol = 1.0e-5))

    # equilibrium gas fraction of crystal-free rhyolite with an H2O gas of Halbach & Chatterjee (1982)
    let ρ = ThreePhase_Density(; ρgas = HalbachChatterjee_Density()), m_eq = compute_dissolved(sy, 2.0e8, 1173.15, 0.0)[1]
        x_g = 0.02 / (1 - m_eq)
        Vg = x_g / HalbachChatterjee_Density()(; P = 2.0e8, T = 1173.15)
        Vm = (1 - x_g) / 2300
        @test compute_gas_fraction(sy, ρ, (; P = 2.0e8, T = 1173.15, ϕ = 1.0, m_h2o_total = m_eq + 0.02)) ≈ Vg / (Vg + Vm)
        # above ~2.5 GPa the dissolved H2O exceeds a mass fraction of 1, which has no gas fraction
        @test compute_dissolved(sy, 3.0e9, 1173.15, 0.0)[1] > 1
        @test isnan(compute_gas_fraction(sy, ρ, (; P = 3.0e9, T = 1173.15, ϕ = 1.0, m_h2o_total = 0.05)))
    end

    @test_throws "oxd_wt must hold 9 oxide wt%" SunYao2026_Solubility(; oxd_wt = (77.0, 0.1, 12.8))
    @test_throws "oxide wt% must be non-negative" SunYao2026_Solubility(; oxd_wt = (77.0, 0.1, 12.8, 0.8, 0.1, 0.6, -4.0, 4.6, 0.0))
    @test_throws "oxide wt% must not all be zero" SunYao2026_Solubility(; oxd_wt = ntuple(_ -> 0.0, 9))

    # Saturation pressure: inverse of compute_dissolved in (P, X_co2) -----
    for law in (s, m, sy), P in (1.0e8, 2.5e8, 4.0e8), X in (0.0, 0.3, 0.6, 0.9), T in (1123.15, 1373.15)
        Pi, Xi = compute_saturation_pressure(law, compute_dissolved(law, P, T, X)..., T)
        @test Pi ≈ P rtol = 1.0e-8
        @test Xi ≈ X atol = 1.0e-8
    end
    # pure CO2 gas (Mafic_Solubility dissolves H2O even then, so it has no m_h2o = 0 inverse)
    for law in (s, sy)
        @test all(compute_saturation_pressure(law, 0.0, compute_dissolved(law, 2.0e8, 1173.15, 1.0)[2], 1173.15) .≈ (2.0e8, 1.0))
    end
    # 4 wt% H2O + 0.1 wt% CO2 saturates at ~3 kbar with a CO2-rich gas (Sun & Yao 2026, sect. 4.2)
    let (P, X) = compute_saturation_pressure(sy, 0.04, 0.001, 1173.15)
        @test 2.6e8 < P < 2.9e8
        @test 0.5 < X < 0.6
    end
    @test all(compute_saturation_pressure(sy_nd, compute_dissolved(sy, 2.0e8, 1173.15, 0.3)..., nondimensionalize(1173.15K, CharDim)) .≈ (nondimensionalize(2.0e8Pa, CharDim), 0.3))
    @test @inferred(compute_saturation_pressure(sy32, compute_dissolved(sy32, 2.0f8, 1173.15f0, 0.3f0)..., 1173.15f0)) isa Tuple{Float32, Float32}
    @test all(isapprox.(compute_saturation_pressure(sy32, compute_dissolved(sy32, 2.0f8, 1173.15f0, 0.3f0)..., 1173.15f0), (2.0e8, 0.3); rtol = 1.0e-4))
    @test_throws DomainError compute_saturation_pressure(sy, -0.01, 0.001, 1173.15)
    @test_throws DomainError compute_saturation_pressure(sy, 0.0, 0.0, 1173.15)
    # the mafic H2O polynomial exceeds 1 wt% at any pressure at 1100 °C
    @test_throws "no convergence" compute_saturation_pressure(m, 0.005, 0.0, 1373.15)

    # empty-parameter fallback
    @test compute_dissolved(No_MaterialParam()) === (0.0, 0.0)

    # Multiphase: NTuple of phases + a phase array ------------------------
    phases = (
        SetMaterialParams(; Phase = 1, Solubility = Liu2005_Solubility()),
        SetMaterialParams(; Phase = 2, Solubility = Mafic_Solubility()),
    )
    Phasearr = [1, 2, 1]
    n = length(Phasearr)
    m_h2o = zeros(n)
    m_co2 = zeros(n)
    args = (; P = fill(2.0e8, n), T = fill(1200.0, n), X_co2 = fill(0.3, n))
    compute_dissolved!(m_h2o, m_co2, phases, Phasearr, args)
    @test m_h2o[1] ≈ compute_dissolved(Liu2005_Solubility(), 2.0e8, 1200.0, 0.3)[1]
    @test m_h2o[2] ≈ compute_dissolved(Mafic_Solubility(), 2.0e8, 1200.0, 0.3)[1]
    @test all(m_co2 .> 0)

    # single-struct in-place, and generic (view-wrapped) output arrays
    mh_v = view(zeros(n + 2), 2:(n + 1))
    mc_v = view(zeros(n + 2), 2:(n + 1))
    compute_dissolved!(mh_v, mc_v, Liu2005_Solubility(), args)
    @test all(mh_v .≈ compute_dissolved(Liu2005_Solubility(), 2.0e8, 1200.0, 0.3)[1])

    # phase-struct scalar dispatch
    @test compute_dissolved(phases[1], (; P = 2.0e8, T = 1200.0, X_co2 = 0.0))[1] ≈ mh

    # integer-phase selection returns the selected phase's tuple
    argp = (; P = 2.0e8, T = 1200.0, X_co2 = 0.3)
    @test all(isapprox.((compute_dissolved(phases, 2, argp)), compute_dissolved(phases[2], argp); rtol = 1.0e-10))

    # phase-ratio mixing: bulk mass-fraction average, both outputs mixed
    ratios = SA[0.25, 0.75]
    mixed = compute_dissolved_ratio(ratios, phases, argp)
    ref_h = 0.25 * compute_dissolved(phases[1], argp)[1] + 0.75 * compute_dissolved(phases[2], argp)[1]
    ref_c = 0.25 * compute_dissolved(phases[1], argp)[2] + 0.75 * compute_dissolved(phases[2], argp)[2]
    @test mixed[1] ≈ ref_h
    @test mixed[2] ≈ ref_c
    @test compute_dissolved_ratio(Tuple(ratios), phases, argp) == mixed   # NTuple ratios too
    # ratios sum to 1 with equal phases -> the phase value itself
    @test all(compute_dissolved_ratio(SA[1.0, 0.0], phases, argp) .≈ compute_dissolved(phases[1], argp))

    # a phase with no solubility model contributes nothing (no BoundsError)
    phases0 = (phases[1], SetMaterialParams(; Phase = 2))
    @test compute_dissolved(phases0, 2, argp) === (0.0, 0.0)
    @test compute_dissolved_ratio(SA[0.5, 0.5], phases0, argp)[1] ≈ 0.5 * compute_dissolved(phases[1], argp)[1]

    # Float32 structs stay Float32 end-to-end
    s32 = Liu2005_Solubility(; Pref = 1.0f6Pa, Tref = 1.0f0K)
    @test @inferred(compute_dissolved(s32, 2.0f8, 1200.0f0, 0.3f0)) isa Tuple{Float32, Float32}
    m32 = Mafic_Solubility(; T0 = 273.15f0K, Tref = 1.0f0K, Pref = 1.0f6Pa)
    @test @inferred(compute_dissolved(m32, 2.0f8, 1200.0f0, 0.3f0)) isa Tuple{Float32, Float32}

    @test @inferred(compute_dissolved(s, 2.0e8, 1200.0, 0.3)) isa Tuple{Float64, Float64}
    @test @inferred(compute_dissolved_ratio(ratios, phases, argp)) isa Tuple{Float64, Float64}

    # Gas-mixture properties ----------------------------------------------
    g = GasMixture()
    @test isbits(g)
    @test isdimensional(g) === true
    @test param_info(g).Equation isa LaTeXString

    @test effective_molar_mass(g, 0.0) ≈ 18.02e-3
    @test effective_molar_mass(g, 1.0) ≈ 44.01e-3
    @test effective_molar_mass(g, 0.5) ≈ (18.02e-3 + 44.01e-3) / 2
    # c_g = 0 at X_co2 = 0 (reference convention)
    @test compute_gas_heatcapacity(g, 0.0) == 0.0
    # mixed value matches the formula
    let X = 0.3
        m_g = 18.02e-3 * (1 - X) + 44.01e-3 * X
        cg = (18.02e-3 * 3880 * (1 - X) + 44.01e-3 * 1200 * X) / m_g
        @test compute_gas_heatcapacity(g, X) ≈ cg
    end
    # nondimensionalization keeps the struct consistent
    g_nd = nondimensionalize(g, CharDim)
    @test isdimensional(g_nd) === false
    @test compute_gas_heatcapacity(g_nd, 0.0) == 0.0

    # Exsolved gas fraction -----------------------------------------------
    # default mixture: constant 2300 kg/m^3 melt, 2600 kg/m^3 crystals, ideal-gas H2O
    ρ3 = ThreePhase_Density()
    m_eq = compute_dissolved(s, 2.0e8, 1200.0, 0.0)[1]
    argsg = (; P = 2.0e8, T = 1200.0, ϕ = 0.6, m_h2o_total = m_eq + 0.02)
    ϕg = compute_gas_fraction(s, ρ3, argsg)
    # volumes per unit mass of melt + gas; crystals fill 0.4/0.6 of the melt volume
    let x_g = 0.02 / (1 - m_eq), Vg = x_g / (2.0e8 / (461.5 * 1200.0)), Vm = (1 - x_g) / 2300
        @test ϕg ≈ Vg / (Vg + Vm + Vm * 0.4 / 0.6)
        @test compute_gas_fraction(s, ρ3, (; argsg..., ϕ = 1.0)) ≈ Vg / (Vg + Vm)   # no crystals
    end
    # with the H2O gas density of Halbach & Chatterjee (1982)
    ρ3_hc = ThreePhase_Density(; ρgas = HalbachChatterjee_Density())
    let x_g = 0.02 / (1 - m_eq), Vg = x_g / HalbachChatterjee_Density()(; P = 2.0e8, T = 1200.0), Vm = (1 - x_g) / 2300
        @test compute_gas_fraction(s, ρ3_hc, (; argsg..., ϕ = 1.0)) ≈ Vg / (Vg + Vm)
    end
    @test compute_gas_fraction(s, ρ3; argsg...) == ϕg
    @test compute_gas_fraction(s, ρ3, (; argsg..., ϕ = 0.0)) == 0.0                # fully crystalline

    # undersaturated: no gas, and the gas law is not evaluated (Redlich-Kwong is NaN below 273.15 K)
    @test compute_gas_fraction(s, ρ3, (; argsg..., m_h2o_total = m_eq - 0.01)) === 0.0
    ρ3_rk = ThreePhase_Density(; ρgas = RedlichKwong_Density())
    @test compute_gas_fraction(s, ρ3_rk, (; argsg..., T = 250.0, m_h2o_total = 0.05)) === 0.0

    # out-of-range inputs give NaN
    @test isnan(compute_gas_fraction(s, ρ3, (; argsg..., m_h2o_total = 1.0)))
    @test isnan(compute_gas_fraction(s, ρ3, (; argsg..., m_h2o_total = -0.01)))
    @test isnan(compute_gas_fraction(s, ρ3, (; argsg..., ϕ = 1.1)))

    # with ϕ_x = (1 - ϕ_gas)(1 - ϕ), the three-phase mixture holds m_h2o_total of water per mass of melt + gas
    ϕg_rk = compute_gas_fraction(m, ρ3_rk, argsg)
    let ϕ_x = (1 - ϕg_rk) * (1 - argsg.ϕ), ρg = compute_density(ρ3_rk.ρgas, argsg)
        ρmix = compute_density(ρ3_rk, (; argsg..., ϕ_gas = ϕg_rk, ϕ_x))
        water = ϕg_rk * ρg + compute_dissolved(m, 2.0e8, 1200.0, 0.0)[1] * (1 - ϕ_x - ϕg_rk) * 2300
        @test water / (ρmix - ϕ_x * 2600) ≈ argsg.m_h2o_total
    end

    # dimensionless: nondimensionalized laws and inputs give the same fraction
    args_nd = (;
        argsg...,
        P = nondimensionalize(2.0e8Pa, CharDim), T = nondimensionalize(1200.0K, CharDim),
    )
    @test compute_gas_fraction(m_nd, nondimensionalize(ρ3_rk, CharDim), args_nd) ≈ ϕg_rk

    # Float32 in, Float32 out (`m` and `g` are bound to laws above, hence GeoParams.m)
    ρ3_32 = ThreePhase_Density(;
        ρmelt = ConstantDensity(ρ = 2300.0f0kg / GeoParams.m^3),
        ρsolid = ConstantDensity(ρ = 2600.0f0kg / GeoParams.m^3),
        ρgas = IdealGas_Density(; Rs = 461.5f0J / kg / K),
        ρ = 2400.0f0kg / GeoParams.m^3,
    )
    args32 = (; P = 2.0f8, T = 1200.0f0, ϕ = 0.6f0, m_h2o_total = Float32(m_eq + 0.02))
    @test @inferred(compute_gas_fraction(s32, ρ3_32, args32)) isa Float32
    @test compute_gas_fraction(s32, ρ3_32, args32) ≈ ϕg rtol = 1.0e-4
    @test @inferred(compute_gas_fraction(s, ρ3, argsg)) isa Float64
end

@testset "SunYao2026_Solubility pinned values" begin
    # independent recomputation of Sun & Yao (2026) eqs. 7–8 and Table 1 (default rhyolite, 900 °C)
    s, T = SunYao2026_Solubility(), 1173.15
    @test compute_dissolved(s, 1.0e9, T, 0.0)[1] ≈ 0.16884 rtol = 1.0e-4
    @test compute_dissolved(s, 1.0e9, T, 1.0)[2] ≈ 4.5693e-3 rtol = 1.0e-4
    m = compute_dissolved(s, 2.0e8, T, 0.3)
    @test m[1] ≈ 0.044631 rtol = 1.0e-4
    @test m[2] ≈ 4.8731e-4 rtol = 1.0e-4
    P, X = compute_saturation_pressure(s, 0.04, 0.001, T)
    @test P ≈ 2.765e8 rtol = 1.0e-4
    @test X ≈ 0.54339 rtol = 1.0e-4
end
