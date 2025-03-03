using Test, GeoParams, StaticArrays, LaTeXStrings

@testset "Permeability.jl" begin

    # Make sure that structs are isbits
    x = ConstantPermeability()
    @test isbits(x)
    @test param_info(x).Equation === L"k = cst"
    @test isdimensional(x) === true
    @test repr("text/plain", x) isa String

    x = HazenPermeability()
    @test isbits(x)
    @test param_info(x).Equation === L"k = C \cdot D_{10}^2"
    @test isdimensional(x) === true

    x = PowerLawPermeability()
    @test isbits(x)
    @test param_info(x).Equation === L"k = c \cdot k_0 \cdot \phi^n"
    @test isdimensional(x) === true

    x = CarmanKozenyPermeability()
    @test isbits(x)
    @test param_info(x).Equation === L"k = c \left(\frac{\phi}{\phi_0}\right)^n"

    # Test the permeability calculations with units
    x1 = ConstantPermeability(; k = 1.0e-12m^2)
    @test x1.k.val == 1.0e-12
    @test GeoParams.get_k(x1) == 1.0e-12
    @test compute_permeability(x1) ≈ 1.0e-12

    x2 = HazenPermeability(; C = 1.0, D10 = 1.0e-4m)
    args = (; ϕ = 0.4)
    @test compute_permeability(x2, args) ≈ 1.0 * (0.1e-3)^2
    @test x2() ≈ 1.0 * (0.1e-3)^2

    x3 = PowerLawPermeability(; c = 1.0, k0 = 1.0e-12m^2, n = 3)
    args = (; ϕ = 0.4)
    @test compute_permeability(x3, args) ≈ 1.0 * 1.0e-12 * (0.4)^3
    @test x3(args) ≈ 1.0 * 1.0e-12 * (0.4)^3

    x4 = CarmanKozenyPermeability(; c = 1.0, ϕ0 = 0.3, n = 3)
    args = (; ϕ = 0.4)
    @test compute_permeability(x4, args) ≈ 1.0 * (0.4 / 0.3)^3
    @test x4(args) ≈ 1.0 * (0.4 / 0.3)^3

    # Test the permeability calculations with non-dimensionalized units
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 1000km)
    x1 = nondimensionalize(x1, CharUnits_GEO)
    @test x1.k.val ≈ 1.0e-24

    x2 = nondimensionalize(x2, CharUnits_GEO)
    @test x2.C.val ≈ 1.0
    @test x2.D10.val ≈ 1.0e-10

    x3 = nondimensionalize(x3, CharUnits_GEO)
    @test x3.k0.val ≈ 1.0e-24

    x4 = nondimensionalize(x4, CharUnits_GEO)
    @test x4.c.val ≈ 1.0


    # Define material parameters with different permeability parameterizations
    Mat_tup = (
        SetMaterialParams(;
            Name = "Mantle", Phase = 1, Permeability = ConstantPermeability(; k = 1.0e-12m^2)
        ),
        SetMaterialParams(; Name = "Crust", Phase = 2, Permeability = HazenPermeability(; C = 1.0, D10 = 1.0e-4m)),
        SetMaterialParams(;
            Name = "UpperCrust",
            Phase = 3,
            Permeability = PowerLawPermeability(; c = 1.0, k0 = 1.0e-12m^2, n = 3, ϕ = 0.1),
            Density = PT_Density(),
        ),
        SetMaterialParams(; Name = "LowerCrust", Phase = 4, Permeability = CarmanKozenyPermeability(; c = 1.0, ϕ0 = 0.3, n = 3), Density = PT_Density()),
    )

    n = 100
    Phases = ones(Int64, n, n, n)
    Phases[:, :, 20:end] .= 2
    Phases[:, :, 50:end] .= 3
    Phases[:, :, 70:end] .= 4

    ϕ = fill(1.0e-2, size(Phases))
    T = fill((800 + 273.15), size(Phases))
    P = fill(10, size(Phases))
    args = (P = P, T = T)

    @test compute_permeability(Mat_tup, Phases[1], args) == 1.0e-12

    compute_permeability!(ϕ, Mat_tup, Phases, args)

    @test sum(ϕ) / n^3 ≈ 1.1484481671481687e-5  # Adjust this value based on expected results

    # Test PhaseRatio and StaticArrays PhaseRatios as input
    args = (P = 0.0, T = 1000.0 + 273.15)
    PhaseRatio = (0.25, 0.25, 0.25, 0.25)
    @test compute_permeability_ratio(PhaseRatio, Mat_tup, args) == 9.26175950925951e-6 # Adjust this value based on expected results

    SvPhaseRatio = SA[0.25, 0.25, 0.25, 0.25]
    @test compute_permeability_ratio(SvPhaseRatio, Mat_tup, args) == 9.26175950925951e-6  # Adjust this value based on expected results

end
