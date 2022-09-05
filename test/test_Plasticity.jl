using Test
using GeoParams

@testset "Plasticity.jl" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=10km)

    # DruckerPrager ---------
    p = DruckerPrager()
    info = param_info(p)
    @test isbits(p)
    @test NumValue(p.ϕ) == 30

    p_nd = p
    p_nd = nondimensionalize(p_nd, CharUnits_GEO)
    @test p_nd.C.val ≈ 1

    # Compute with dimensional units
    τII = 20e6
    P = 1e6
    args = (P=P, τII=τII)
    @test compute_yieldfunction(p, args) ≈ 1.0839745962155614e7      # no Pfluid

    args_f = (P=P, τII=τII, Pf=0.5e6)
    @test compute_yieldfunction(p, args_f) ≈ 1.1089745962155614e7    # with Pfluid

    # Test with arrays
    P_array = ones(10) * 1e6
    τII_array = ones(10) * 20e6
    F_array = similar(P_array)
    compute_yieldfunction!(F_array, p, (; P=P_array, τII=τII_array))
    @test F_array[1] ≈ 1.0839745962155614e7

    Pf_array = ones(10) * 0.5e6
    Ff_array = similar(P_array)
    compute_yieldfunction!(Ff_array, p, (; P=P_array, τII=τII_array, Pf=Pf_array))
    @test Ff_array[1] ≈ 1.1089745962155614e7

    # Check that it works if we give a phase array
    MatParam = (
        SetMaterialParams(; Name="Mantle", Phase=1, Plasticity=DruckerPrager()),
        SetMaterialParams(; Name="Crust", Phase=2, Plasticity=DruckerPrager(; ϕ=10)),
        SetMaterialParams(;
            Name="Crust", Phase=3, HeatCapacity=ConstantHeatCapacity(; cp=1100J / kg / K)
        ),
    )

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    Phases[:, :, 20:end] .= 2
    Phases[:, :, 60:end] .= 2

    τII = ones(size(Phases)) * 10e6
    P = ones(size(Phases)) * 1e6
    Pf = ones(size(Phases)) * 0.5e6
    F = zero(P)
    args = (P=P, τII=τII)
    compute_yieldfunction!(F, MatParam, Phases, args)    # computation routine w/out P (not used in most heat capacity formulations)     
    @test maximum(F[1, 1, :]) ≈ 839745.962155614

    args_f = (P=P, τII=τII, Pf=Pf)
    Ff = zero(P)
    compute_yieldfunction!(Ff, MatParam, Phases, args_f)    # computation routine w/out P (not used in most heat capacity formulations)     

    # test if we provide phase ratios
    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end
    compute_yieldfunction!(F, MatParam, PhaseRatio, args)
    num_alloc = @allocated compute_yieldfunction!(F, MatParam, PhaseRatio, args)
    @test maximum(F[1, 1, :]) ≈ 839745.962155614
    # @test num_alloc <= 32

    # -----------------------

end
