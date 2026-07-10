using Test, GeoParams, StaticArrays, LaTeXStrings
using ForwardDiff: derivative

@testset "GasDensity.jl" begin

    # Redlich-Kwong gas density -------------------------------------------
    rk = RedlichKwong_Density()
    @test isbits(rk)
    @test isdimensional(rk) === true
    @test param_info(rk).Equation isa LaTeXString
    @test sprint(show, rk) isa String
    @test occursin("Redlich-Kwong", sprint(show, rk))

    # reference value: ~360 kg/m³ at 200 MPa, 1200 K
    @test rk(; P = 2.0e8, T = 1200.0) ≈ 3.6e2 rtol = 0.1
    @test compute_density(rk, (; P = 2.0e8, T = 1200.0)) ≈ rk(; P = 2.0e8, T = 1200.0)
    @test rk((; P = 2.0e8, T = 1200.0)) ≈ rk(; P = 2.0e8, T = 1200.0)

    # monotonicity (plan acceptance)
    @test rk(; P = 1.0e8, T = 1123.15) > rk(; P = 1.0e8, T = 1173.15)   # hotter -> lighter
    @test rk(; P = 3.0e8, T = 1123.15) > rk(; P = 1.0e8, T = 1123.15)   # higher P -> denser

    # dimensional Quantity input hits the @unpack_units branch
    @test ustrip(kg / m^3, rk(; P = 200.0MPa, T = 1200.0K)) ≈ rk(; P = 2.0e8, T = 1200.0) rtol = 1.0e-6

    # custom coefficients through the keyword constructor
    rk2 = RedlichKwong_Density(; coeffs = (-100.0, 120.0, 110.0))
    @test rk2.coeffs == (-100.0, 120.0, 110.0)

    # nondimensionalization round-trip: dimensionalized result matches
    CharDim = GEO_units(; viscosity = 1.0e19, length = 1000km)
    rk_nd = nondimensionalize(rk, CharDim)
    @test isdimensional(rk_nd) === false
    ρ_nd = rk_nd(;
        P = nondimensionalize(2.0e8Pa, CharDim),
        T = nondimensionalize(1200.0K, CharDim),
    )
    @test ustrip(dimensionalize(ρ_nd, kg / m^3, CharDim)) ≈ rk(; P = 2.0e8, T = 1200.0)

    # Ideal-gas density ----------------------------------------------------
    ig = IdealGas_Density()
    @test isbits(ig)
    @test isdimensional(ig) === true
    @test param_info(ig).Equation === L"\rho_g = P/(R_s T)"
    @test sprint(show, ig) isa String

    # ρ = P/(Rs T)
    @test ig(; P = 1.0e8, T = 1000.0) ≈ 1.0e8 / (461.5 * 1000.0)
    @test compute_density(ig, (; P = 1.0e8, T = 1000.0)) ≈ 1.0e8 / (461.5 * 1000.0)
    @test ig((; P = 1.0e8, T = 1000.0)) ≈ 1.0e8 / (461.5 * 1000.0)

    # exact analytic derivative ∂ρ/∂P = 1/(Rs T)
    @test derivative(P -> ig(; P, T = 1000.0), 1.0e8) ≈ 1.0 / (461.5 * 1000.0)

    # Quantity branch
    @test ustrip(kg / m^3, ig(; P = 100.0MPa, T = 1000.0K)) ≈ ig(; P = 1.0e8, T = 1000.0) rtol = 1.0e-6

    # nondim round-trip (homogeneous -> exact)
    ig_nd = nondimensionalize(ig, CharDim)
    @test isdimensional(ig_nd) === false
    ρig_nd = ig_nd(;
        P = nondimensionalize(1.0e8Pa, CharDim),
        T = nondimensionalize(1000.0K, CharDim),
    )
    @test ustrip(dimensionalize(ρig_nd, kg / m^3, CharDim)) ≈ ig(; P = 1.0e8, T = 1000.0)

    # Three-phase density --------------------------------------------------
    tp = ThreePhase_Density(;
        ρmelt = ConstantDensity(ρ = 2300kg / m^3),
        ρx = ConstantDensity(ρ = 2600kg / m^3),
        ρgas = ConstantDensity(ρ = 1kg / m^3),
    )
    @test isbits(tp)
    @test isdimensional(tp) === true
    @test param_info(tp).Equation isa LaTeXString
    @test occursin("Three-phase", sprint(show, tp))

    # volume-fraction mixture: 40% crystals, no gas
    @test tp(; ϕ_gas = 0.0, ϕ_x = 0.4) ≈ 0.6 * 2300 + 0.4 * 2600
    @test compute_density(tp, (; ϕ_gas = 0.0, ϕ_x = 0.4)) ≈ 0.6 * 2300 + 0.4 * 2600
    @test tp((; ϕ_gas = 0.0, ϕ_x = 0.4)) ≈ 0.6 * 2300 + 0.4 * 2600

    # adding gas lowers the density
    @test tp(; ϕ_gas = 0.1, ϕ_x = 0.4) < tp(; ϕ_gas = 0.0, ϕ_x = 0.4)

    # default gas is Redlich-Kwong: mixture responds to P, T through the EOS
    tp_rk = ThreePhase_Density()
    ρmix = tp_rk(; ϕ_gas = 0.2, ϕ_x = 0.3, P = 2.0e8, T = 1200.0)
    @test 0 < ρmix < 2600

    # nondim round-trip on the three-phase mixture (constant sub-densities)
    tp_nd = nondimensionalize(tp, CharDim)
    @test isdimensional(tp_nd) === false
    ρtp_nd = tp_nd(; ϕ_gas = 0.0, ϕ_x = 0.4)
    @test ustrip(dimensionalize(ρtp_nd, kg / m^3, CharDim)) ≈ tp(; ϕ_gas = 0.0, ϕ_x = 0.4)

    # SetMaterialParams accepts the new densities in the Density slot
    ph = SetMaterialParams(; Phase = 1, Density = RedlichKwong_Density())
    @test compute_density(ph, (; P = 2.0e8, T = 1200.0)) ≈ rk(; P = 2.0e8, T = 1200.0)
end
