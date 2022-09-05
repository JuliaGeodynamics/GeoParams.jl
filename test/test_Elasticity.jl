using Test
using GeoParams

@testset "Elasticity.jl" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=10km)

    # ConstantElasticity ---------
    p = ConstantElasticity()
    info = param_info(p)
    @test isbits(p)
    @test NumValue(p.G) == 5e10

    p_nd = p
    p_nd = nondimensionalize(p_nd, CharUnits_GEO)
    @test p_nd.G.val ≈ 5000.0

    # Compute with dimensional units
    τII = 20e6
    τII_old = 15e6
    dt = 1e6
    args = (τII=τII, τII_old=τII_old, dt=dt)
    @test compute_elastic_shear_strainrate(p, args) ≈ 5.0e-11  # compute

    # Test with arrays
    τII_old_array = ones(10) * 15e6
    τII_array = ones(10) * 20e6
    ε_el_array = similar(τII_array)
    args = (τII=τII_array, τII_old=τII_old_array, dt=dt)
    compute_elastic_shear_strainrate!(ε_el_array, p, args)
    @test ε_el_array[1] ≈ 5.0e-11

    # Check that it works if we give a phase array
    MatParam = (
        SetMaterialParams(; Name="Mantle", Phase=1, Elasticity=ConstantElasticity()),
        SetMaterialParams(;
            Name="Crust", Phase=2, Elasticity=ConstantElasticity(; G=1e10Pa)
        ),
        SetMaterialParams(;
            Name="Crust", Phase=2, HeatCapacity=ConstantHeatCapacity(; cp=1100J / kg / K)
        ),
    )

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    Phases[:, :, 20:end] .= 2
    Phases[:, :, 60:end] .= 2

    τII = ones(size(Phases)) * 20e6
    τII_old = ones(size(Phases)) * 15e6
    ε_el = zero(τII)
    args = (τII=τII, τII_old=τII_old, dt=1e6)
    compute_elastic_shear_strainrate!(ε_el, MatParam, Phases, args)    # computation routine w/out P (not used in most heat capacity formulations)     
    @test maximum(ε_el[1, 1, :]) ≈ 2.5e-10

    # test if we provide phase ratios
    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end

    # Note; using PhaseRatio currently requires an array of timesteps - can probably be fixed to also allow scalars
    dt_arr = ones(size(Phases)) * 1e6     # needs to be an array of timesteps currently
    args = (τII=τII, τII_old=τII_old, dt=dt_arr)
    compute_elastic_shear_strainrate!(ε_el, MatParam, PhaseRatio, args)
    @test maximum(ε_el[1, 1, :]) ≈ 2.5e-10

    args1 = (τII=τII, τII_old=τII_old, dt=1e6)
    compute_elastic_shear_strainrate!(ε_el, MatParam, PhaseRatio, args1)
    @test maximum(ε_el[1, 1, :]) ≈ 2.5e-10

    num_alloc = @allocated compute_elastic_shear_strainrate!(
        ε_el, MatParam, PhaseRatio, args
    )
    @test maximum(ε_el[1, 1, :]) ≈ 2.5e-10
    # @test num_alloc <= 32

    # -----------------------

end
