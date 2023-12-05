using Test, GeoParams, AllocCheck

@testset "Plasticity.jl" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=10km)

    # DruckerPrager ---------
    p = DruckerPrager()
    info = param_info(p)
    @test isbits(p)
    @test NumValue(p.ϕ) == 30
    @test isvolumetric(p) == false

    p_nd = p
    p_nd = nondimensionalize(p_nd, CharUnits_GEO)
    @test p_nd.C.val ≈ 1

    # Compute with dimensional units
    τII = 20e6
    P = 1e6
    args = (P=P, τII=τII)
    args1 = (τII=τII, P=P)
    @test compute_yieldfunction(p, args) ≈ 1.0839745962155614e7      # no Pfluid
    @test compute_yieldfunction(p, args) ≈ compute_yieldfunction(p, args1)    # different order

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
    args_f1 = (Pf=Pf, τII=τII, P=P)

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
    compute_yieldfunction!(F, MatParam, PhaseRatio, args)
    @test maximum(F[1, 1, :]) ≈ 839745.962155614
    num_alloc = check_allocs(compute_yieldfunction!, typeof.((F, MatParam, PhaseRatio, args)))
    @test isempty(num_alloc)

    # Test plastic potential derivatives
    ## 2D
    τij = (1.0, 2.0, 3.0)
    fxx(τij) = 0.5 * τij[1] / second_invariant(τij)
    fyy(τij) = 0.5 * τij[2] / second_invariant(τij)
    fxy(τij) = τij[3] / second_invariant(τij)
    solution2D = [fxx(τij), fyy(τij), fxy(τij)]

    # # using StaticArrays
    # τij_static = @SVector [1.0, 2.0, 3.0]
    # out1 = ∂Q∂τ(p, τij_static)
    # @test out1 == solution2D
    # @test compute_plasticpotentialDerivative(p, τij_static) == ∂Q∂τ(p, τij_static)

    # using tuples
    τij_tuple = (1.0, 2.0, 3.0)
    out2 = ∂Q∂τ(p, τij_tuple)
    @test out2 == Tuple(solution2D)
    @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(p, τij_tuple)

    # using AD
    Q = second_invariant # where second_invariant is a function
    # ad1 = ∂Q∂τ(Q, τij_static)
    # @test out1 == solution2D
    # @test compute_plasticpotentialDerivative(p, τij_static) == ∂Q∂τ(p, τij_static)
    ad2 = ∂Q∂τ(Q, τij_tuple)
    @test out2 == Tuple(solution2D)
    @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(p, τij_tuple)

    ## 3D
    τij = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    gxx(τij) = 0.5 * τij[1] / second_invariant(τij)
    gyy(τij) = 0.5 * τij[2] / second_invariant(τij)
    gzz(τij) = 0.5 * τij[3] / second_invariant(τij)
    gyz(τij) = τij[4] / second_invariant(τij)
    gxz(τij) = τij[5] / second_invariant(τij)
    gxy(τij) = τij[6] / second_invariant(τij)
    solution3D = [gxx(τij), gyy(τij), gzz(τij), gyz(τij), gxz(τij), gxy(τij)]

    # # using StaticArrays
    # τij_static = @SVector [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    # out3 = ∂Q∂τ(p, τij_static)
    # @test out3 == solution3D
    # @test compute_plasticpotentialDerivative(p, τij_static) == ∂Q∂τ(p, τij_static)

    # using tuples
    τij_tuple = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    out4 = ∂Q∂τ(p, τij_tuple)
    @test out4 == Tuple(solution3D)
    @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(p, τij_tuple)

    # using AD
    Q = second_invariant # where second_invariant is a function
    # ad3 = ∂Q∂τ(Q, τij_static)
    # @test out3 == solution3D
    # @test compute_plasticpotentialDerivative(p, τij_static) == ∂Q∂τ(p, τij_static)
    ad4 = ∂Q∂τ(Q, τij_tuple)
    @test out4 == Tuple(solution3D)
    @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(p, τij_tuple)

    # -----------------------

    # composite rheology with plasticity
    η, G = 10, 1
    t_M = η / G
    εII = 1.0
    args = (;)
    pl2 = DruckerPrager(; C=η, ϕ=0)                # plasticity
    c_pl = CompositeRheology(
        LinearViscous(; η=η * Pa * s), ConstantElasticity(; G=G * Pa), pl2
    ) # linear VEP
    c_pl2 = CompositeRheology(ConstantElasticity(; G=G * Pa), pl2) # linear VEP

    # case where old stress is below yield & new stress is above
    args = (τII_old=9.8001101017963, P=0.0, τII=9.8001101017963)
    F_old = compute_yieldfunction(c_pl.elements[3], args)

    # 
    τ1, = local_iterations_εII(c_pl, εII, args; verbose=false, max_iter=10)
    τ2, = compute_τII(c_pl, εII, args; verbose=false)
    @test τ1 == τ2

    args = merge(args, (τII=τ1,))
    F_check = compute_yieldfunction(c_pl.elements[3], args)
    @test abs(F_check) < 1e-12
end
