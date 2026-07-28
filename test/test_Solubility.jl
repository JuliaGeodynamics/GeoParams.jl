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
    @test all(isapprox.((compute_dissolved(phases, 2, argp)),  compute_dissolved(phases[2], argp); rtol = 1e-10))

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
end
