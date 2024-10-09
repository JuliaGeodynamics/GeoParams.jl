using Test, GeoParams, Unitful, StaticArrays, LaTeXStrings

@testset "Permeability.jl" begin

    # Set alias for permeability function
    if !isdefined(Main, :GeoParamsAliases)
        eval(:(@use GeoParamsAliases permeability = k))
    end

    # Make sure that structs are isbits
    x = ConstantPermeability()
    @test isbits(x)
    @test param_info(x).Equation === L"k = cst"
    @test isdimensional(x) === true

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
    x1 = ConstantPermeability(; k=1e-12m^2)
    @test x1.k.val == 1e-12
    @test GeoParams.get_k(x1) == 1e-12

    x2 = HazenPermeability(; C=1.0, D10=1e-4m)
    args = (;ϕ=0.4)
    @test compute_permeability(x2, args) ≈ 1.0 * (0.1e-3)^2
    @test x2() ≈ 1.0 * (0.1e-3)^2

    x3 = PowerLawPermeability(; c=1.0, k0=1e-12m^2, n=3.0)
    args = (;ϕ=0.4)
    @test compute_permeability(x3, args) ≈ 1.0 * 1.0e-12 * (0.4)^3
    @test x3(args) ≈ 1.0 * 1.0e-12 * (0.4)^3

    x4 = CarmanKozenyPermeability(; c=1.0, ϕ0=0.3, n=3.0)
    args = (;ϕ=0.4)
    @test compute_permeability(x4, args) ≈ 1.0 * (0.4 / 0.3)^3
    @test x4(args) ≈ 1.0 * (0.4 / 0.3)^3

    # Test the permeability calculations with non-dimensionalized units
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=1000km)
    x1 = nondimensionalize(x1, CharUnits_GEO)
    @test x1.k.val ≈ 1e-24

    x2 = nondimensionalize(x2, CharUnits_GEO)
    @test x2.C.val ≈ 1.0
    @test x2.D10.val ≈ 1e-10

    x3 = nondimensionalize(x3, CharUnits_GEO)
    @test x3.k0.val ≈ 1e-24

    x4 = nondimensionalize(x4, CharUnits_GEO)
    @test x4.c.val ≈ 1.0

    # # Test allocations
    # k = [0.0]
    # ϕ = 0.4
    # args = (ϕ=ϕ)

    # # Test allocations using k alias
    # k!(k, x3, args)
    # num_alloc = @allocated k!(k, x3, args)
    # @test num_alloc == 0

    # k!(k, x4, args)
    # num_alloc = @allocated k!(k, x4, args)
    # @test num_alloc == 0

end
