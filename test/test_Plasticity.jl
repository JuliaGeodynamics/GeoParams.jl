using Test
using GeoParams

@testset "Plasticity.jl" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 10km)
    @testset "DruckerPrager" begin
        # DruckerPrager ---------
        p = DruckerPrager()
        info = param_info(p)
        @test isbits(p)
        @test NumValue(p.ϕ) == 30
        @test isvolumetric(p) == false
        @test repr("text/plain", p) isa String

        p_nd = p
        p_nd = nondimensionalize(p_nd, CharUnits_GEO)
        @test p_nd.C.val ≈ 1


        # Compute with dimensional units
        τII = 20.0e6
        P = 1.0e6
        args = (P = P, τII = τII)
        args1 = (τII = τII, P = P)
        @test compute_yieldfunction(p, args) ≈ 1.0839745962155614e7            # no Pfluid
        @test compute_yieldfunction(p, args) ≈ compute_yieldfunction(p, args1) # different order


        # test cohesion perturbation
        @test compute_yieldfunction(p, (τII = τII, P = P, perturbation_C = 1.0)) ≈ 1.0839745962155614e7  # no perturbation
        @test compute_yieldfunction(p, (τII = τII, P = P, perturbation_C = 0.0)) == 1.95e7                # cohesion = 0.0
        @test compute_yieldfunction(p, (τII = τII, P = P, perturbation_C = 0.5)) ≈ 1.5169872981077807e7  # midway

        args_f = (P = P, τII = τII, Pf = 0.5e6)

        @test compute_yieldfunction(p, args_f) ≈ 1.1089745962155614e7    # with Pfluid

        # Test with arrays
        P_array = fill(1.0e6, 10)
        τII_array = fill(20.0e6, 10)
        F_array = similar(P_array)
        compute_yieldfunction!(F_array, p, (; P = P_array, τII = τII_array))
        @test F_array[1] ≈ 1.0839745962155614e7

        Pf_array = ones(10) * 0.5e6
        Ff_array = similar(P_array)
        compute_yieldfunction!(Ff_array, p, (; P = P_array, τII = τII_array, Pf = Pf_array))
        @test Ff_array[1] ≈ 1.1089745962155614e7

        # Check that it works if we give a phase array
        MatParam = (
            SetMaterialParams(; Name = "Mantle", Phase = 1, Plasticity = DruckerPrager()),
            SetMaterialParams(; Name = "Crust", Phase = 2, Plasticity = DruckerPrager(; ϕ = 10)),
            SetMaterialParams(;
                Name = "Crust", Phase = 3, HeatCapacity = ConstantHeatCapacity(; Cp = 1100J / kg / K)
            ),
        )

        # test computing material properties
        n = 100
        Phases = ones(Int64, n, n, n)
        Phases[:, :, 20:end] .= 2
        Phases[:, :, 60:end] .= 2

        τII = fill(10.0e6, size(Phases)...)
        P = fill(1.0e6, size(Phases)...)
        Pf = fill(0.5e6, size(Phases)...)
        F = zero(P)
        args = (P = P, τII = τII)
        compute_yieldfunction!(F, MatParam, Phases, args)    # computation routine w/out P (not used in most heat capacity formulations)
        @test maximum(F[1, 1, :]) ≈ 839745.962155614

        args_f = (P = P, τII = τII, Pf = Pf)
        args_f1 = (Pf = Pf, τII = τII, P = P)

        Ff = zero(P)
        compute_yieldfunction!(Ff, MatParam, Phases, args_f)    # computation routine w/out P (not used in most heat capacity formulations)
        @test maximum(Ff[1, 1, :]) ≈ 1.089745962155614e6

        # test if we provide phase ratios
        PhaseRatio = zeros(n, n, n, 3)
        for i in CartesianIndices(Phases)
            j = Phases[i]
            I = CartesianIndex(i, j)
            PhaseRatio[I] = 1.0
        end
        compute_yieldfunction!(F, MatParam, PhaseRatio, args)
        num_alloc = @allocated compute_yieldfunction!(F, MatParam, PhaseRatio, args)
        @test maximum(F[1, 1, :]) ≈ 839745.962155614
        @test num_alloc == 0

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
        ad2 = ∂Q∂τ(Q, τij_tuple)
        @test ad2 == (0.5, 0.625, 0.75)

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
        ad4 = ∂Q∂τ(Q, τij_tuple)
        @test ad4 == Tuple(solution3D)
        @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(Q, τij_tuple)

        # -----------------------

        # composite rheology with plasticity
        η, G = 10, 1
        t_M = η / G
        εII = 1.0
        args = (;)
        pl2 = DruckerPrager(; C = η, ϕ = 0)                # plasticity
        c_pl = CompositeRheology(LinearViscous(; η = η * Pa * s), ConstantElasticity(; G = G * Pa), pl2) # linear VEP
        c_pl2 = CompositeRheology(ConstantElasticity(; G = G * Pa), pl2) # linear VEP

        # case where old stress is below yield & new stress is above
        args = (τII_old = 9.8001101017963, P = 0.0, τII = 9.8001101017963)
        F_old = compute_yieldfunction(c_pl.elements[3], args)

        #
        τ1, = local_iterations_εII(c_pl, εII, args; verbose = false, max_iter = 10)
        τ2, = compute_τII(c_pl, εII, args; verbose = false)
        @test τ1 == τ2

        args = merge(args, (τII = τ1,))
        F_check = compute_yieldfunction(c_pl.elements[3], args)
        @test abs(F_check) < 1.0e-12
    end

    @testset "DruckerPrager_regularised" begin

        # DruckerPrager_regularised ---------
        p = DruckerPrager_regularised()
        info = param_info(p)
        @test isbits(p)
        @test NumValue(p.ϕ) == 30
        @test isvolumetric(p) == false
        @test repr("text/plain", p) isa String


        p_nd = p
        p_nd = nondimensionalize(p_nd, CharUnits_GEO)

        @test p_nd.C.val ≈ 1

        # Compute with dimensional units
        τII = 20.0e6
        P = 1.0e6
        args = (P = P, τII = τII)
        args1 = (τII = τII, P = P)
        @test compute_yieldfunction(p, args) ≈ 1.0839745962155614e7            # no Pfluid
        @test compute_yieldfunction(p, args) ≈ compute_yieldfunction(p, args1) # different order

        # test cohesion perturbation
        @test compute_yieldfunction(p, (τII = τII, P = P, perturbation_C = 1.0)) ≈ 1.0839745962155614e7  # no perturbation
        @test compute_yieldfunction(p, (τII = τII, P = P, perturbation_C = 0.0)) == 1.95e7                # cohesion = 0.0
        @test compute_yieldfunction(p, (τII = τII, P = P, perturbation_C = 0.5)) ≈ 1.5169872981077807e7  #

        args_f = (P = P, τII = τII, Pf = 0.5e6)

        @test compute_yieldfunction(p, args_f) ≈ 1.1089745962155614e7    # with Pfluid

        # test with arrays
        P_array = fill(1.0e6, 10)
        τII_array = fill(20.0e6, 10)
        F_array = similar(P_array)
        compute_yieldfunction!(F_array, p, (; P = P_array, τII = τII_array))
        @test F_array[1] ≈ 1.0839745962155614e7

        Pf_array = ones(10) * 0.5e6
        Ff_array = similar(P_array)
        compute_yieldfunction!(Ff_array, p, (; P = P_array, τII = τII_array, Pf = Pf_array))
        @test Ff_array[1] ≈ 1.1089745962155614e7

        # Check that it works if we give a phase array
        MatParam = (
            SetMaterialParams(; Name = "Mantle", Phase = 1, Plasticity = DruckerPrager_regularised()),
            SetMaterialParams(; Name = "Crust", Phase = 2, Plasticity = DruckerPrager_regularised(; ϕ = 10)),
            SetMaterialParams(;
                Name = "Crust", Phase = 3, HeatCapacity = ConstantHeatCapacity(; Cp = 1100J / kg / K)
            ),
        )

        # test computing material properties
        n = 100
        Phases = ones(Int64, n, n, n)
        Phases[:, :, 20:end] .= 2
        Phases[:, :, 60:end] .= 2

        τII = fill(10.0e6, size(Phases)...)
        P = fill(1.0e6, size(Phases)...)
        Pf = fill(0.5e6, size(Phases)...)
        F = zero(P)
        args = (P = P, τII = τII)
        compute_yieldfunction!(F, MatParam, Phases, args)    # computation routine w/out P (not used in most heat capacity formulations)
        @test maximum(F[1, 1, :]) ≈ 839745.962155614

        args_f = (P = P, τII = τII, Pf = Pf)
        args_f1 = (Pf = Pf, τII = τII, P = P)

        Ff = zero(P)
        compute_yieldfunction!(Ff, MatParam, Phases, args_f)
        @test maximum(Ff[1, 1, :]) ≈ 1.089745962155614e6

        # test if we provide phase ratios
        PhaseRatio = zeros(n, n, n, 3)
        for i in CartesianIndices(Phases)
            j = Phases[i]
            I = CartesianIndex(i, j)
            PhaseRatio[I] = 1.0
        end
        compute_yieldfunction!(F, MatParam, PhaseRatio, args)
        num_alloc = @allocated compute_yieldfunction!(F, MatParam, PhaseRatio, args)
        @test maximum(F[1, 1, :]) ≈ 839745.962155614
        # @test num_alloc <= 32

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
        ad2 = ∂Q∂τ(Q, τij_tuple)
        @test ad2 == (0.5, 0.625, 0.75)

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
        ad4 = ∂Q∂τ(Q, τij_tuple)
        @test ad4 == Tuple(solution3D)
        @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(Q, τij_tuple)

        # -----------------------

        # composite rheology with plasticity
        η, G = 10, 1
        t_M = η / G
        εII = 1.0
        args = (;)
        pl2 = DruckerPrager_regularised(; C = η, ϕ = 00)                # plasticity
        c_pl = CompositeRheology(LinearViscous(; η = η * Pa * s), ConstantElasticity(; G = G * Pa), pl2) # linear VEP
        c_pl2 = CompositeRheology(ConstantElasticity(; G = G * Pa), pl2) # linear VEP

        # case where old stress is below yield & new stress is above
        args = (τII_old = 9.0, P = 0.0, τII = 9.0)
        F_old = compute_yieldfunction(c_pl.elements[3], args)

        #
        τ1, = local_iterations_εII(c_pl, εII, args; verbose = false, max_iter = 10)
        τ2, = compute_τII(c_pl, εII, args; verbose = false)
        @test τ1 == τ2

        args = merge(args, (τII = τ1,))
        F_check = compute_yieldfunction(c_pl.elements[3], args)
        @test abs(F_check) < 1.0e-12
    end


    @testset "DruckerPragerCap" begin

        # DruckerPragerCap ---------
        p = DruckerPragerCap()
        info = param_info(p)
        @test isbits(p)
        @test NumValue(p.ϕ) == 30
        @test NumValue(p.pT) == -1.0e5
        @test isvolumetric(p) == false
        @test repr("text/plain", p) isa String

        p_nd = p
        p_nd = nondimensionalize(p_nd, CharUnits_GEO)
        @test p_nd.C.val ≈ 1

        #Compute with dimensional units
        τII = 20.0e6
        P = 1.0e6
        args = (P = P, τII = τII)
        args1 = (τII = τII, P = P)
        @test compute_yieldfunction(p, args) ≈ 1.0991085617509069e7
        @test compute_yieldfunction(p, args) ≈ compute_yieldfunction(p, args1) # different order

        # test cohesion perturbation
        @test compute_yieldfunction(p, (τII = τII, P = P, perturbation_C = 1.0)) ≈ 1.0991085617509069e7  # no perturbation
        @test compute_yieldfunction(p, (τII = τII, P = P, perturbation_C = 0.0)) ≈ 2.2490074976691682e7  # cohesion = 0.0
        @test compute_yieldfunction(p, (τII = τII, P = P, perturbation_C = 0.5)) ≈ 1.5169872981077807e7  # midway

        args_f = (P = P, τII = τII, Pf = 0.5e6)

        @test compute_yieldfunction(p, args_f) ≈ 1.1297073589335589e7    # with Pfluid

        # test with arrays
        P_array = fill(1.0e6, 10)
        τII_array = fill(20.0e6, 10)
        F_array = similar(P_array)
        compute_yieldfunction!(F_array, p, (; P = P_array, τII = τII_array))
        @test F_array[1] ≈ 1.0991085617509069e7

        Pf_array = ones(10) * 0.5e6
        Ff_array = similar(P_array)
        compute_yieldfunction!(Ff_array, p, (; P = P_array, τII = τII_array, Pf = Pf_array))
        @test Ff_array[1] ≈ 1.1297073589335589e7  # with Pfluid

        # Check that it works if we give a phase array
        MatParam = (
            SetMaterialParams(; Name = "Mantle", Phase = 1, Plasticity = DruckerPragerCap()),
            SetMaterialParams(; Name = "Crust", Phase = 2, Plasticity = DruckerPragerCap(; ϕ = 10)),
            SetMaterialParams(;
                Name = "Crust", Phase = 3, HeatCapacity = ConstantHeatCapacity(; Cp = 1100J / kg / K)
            ),
        )

        # test computing material properties
        n = 100
        Phases = ones(Int64, n, n, n)
        Phases[:, :, 20:end] .= 2
        Phases[:, :, 60:end] .= 2

        τII = fill(10.0e6, size(Phases)...)
        P = fill(1.0e6, size(Phases)...)
        Pf = fill(0.5e6, size(Phases)...)
        F = zero(P)
        args = (P = P, τII = τII)
        compute_yieldfunction!(F, MatParam, Phases, args)    # computation routine w/out P (not used in most heat capacity formulations)
        @test maximum(F[1, 1, :]) ≈ 2.919742711754192e6

        args_f = (P = P, τII = τII, Pf = Pf)
        args_f1 = (Pf = Pf, τII = τII, P = P)

        Ff = zero(P)
        compute_yieldfunction!(Ff, MatParam, Phases, args_f)
        @test maximum(Ff[1, 1, :]) ≈ 3.2926429126483365e6

        # test if we provide phase ratios
        PhaseRatio = zeros(n, n, n, 3)
        for i in CartesianIndices(Phases)
            j = Phases[i]
            I = CartesianIndex(i, j)
            PhaseRatio[I] = 1.0
        end
        compute_yieldfunction!(F, MatParam, PhaseRatio, args)
        num_alloc = @allocated compute_yieldfunction!(F, MatParam, PhaseRatio, args)
        @test maximum(F[1, 1, :]) ≈ 2.919742711754192e6
        @test num_alloc == 0

        # Test plastic potential derivatives
        ## 2D
        τij = (1.0, 2.0, 3.0)
        τII_test = second_invariant(τij)
        dQτII = ∂Q∂τII(p, τII_test)
        fxx(τij) = dQτII * τij[1] / second_invariant(τij)
        fyy(τij) = dQτII * τij[2] / second_invariant(τij)
        fxy(τij) = 2 * dQτII * τij[3] / second_invariant(τij)
        solution2D = [fxx(τij), fyy(τij), fxy(τij)]

        # using tuples
        τij_tuple = (1.0, 2.0, 3.0)
        out2 = ∂Q∂τ(p, τij_tuple)
        @test out2 == Tuple(solution2D)
        @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(p, τij_tuple)

        # using AD
        Q = second_invariant # where second_invariant is a function
        ad2 = ∂Q∂τ(Q, τij_tuple)
        @test ad2 == (0.5, 0.625, 0.75)

        ## 3D
        τij = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
        dQτII = ∂Q∂τII(p, second_invariant(τij))
        gxx(τij) = dQτII * τij[1] / second_invariant(τij)
        gyy(τij) = dQτII * τij[2] / second_invariant(τij)
        gzz(τij) = dQτII * τij[3] / second_invariant(τij)
        gyz(τij) = 2 * dQτII * τij[4] / second_invariant(τij)
        gxz(τij) = 2 * dQτII * τij[5] / second_invariant(τij)
        gxy(τij) = 2 * dQτII * τij[6] / second_invariant(τij)
        solution3D = [gxx(τij), gyy(τij), gzz(τij), gyz(τij), gxz(τij), gxy(τij)]

        # using tuples
        τij_tuple = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
        out4 = ∂Q∂τ(p, τij_tuple)
        @test out4 == Tuple(solution3D)
        @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(p, τij_tuple)

        # using AD
        Q = second_invariant # where second_invariant is a function
        ad4 = ∂Q∂τ(Q, τij_tuple)
        @test all(isapprox.(2 .* dQτII .* ad4, Tuple(solution3D); rtol = 1.0e-5))
        @test all(isapprox.(compute_plasticpotentialDerivative(p, τij_tuple), 2 .* dQτII .* ad4; rtol = 1.0e-5))

        # -----------------------

        # composite rheology with plasticity
        η, G = 10, 1
        t_M = η / G
        εII = 1.0
        args = (;)
        pl2 = DruckerPragerCap(C = η, ϕ = 00, pT = -1)                # plasticity
        c_pl = CompositeRheology(LinearViscous(; η = η * Pa * s), ConstantElasticity(; G = G * Pa), pl2) # linear VEP
        c_pl2 = CompositeRheology(ConstantElasticity(; G = G * Pa), pl2) # linear VEP

        # case where old stress is below yield & new stress is above
        args = (τII_old = 9.8001101017963, P = 0.0, τII = 20.8001101017963)
        F_old = compute_yieldfunction(c_pl.elements[3], args)
        #
        τ1, = local_iterations_εII(c_pl, εII, args; verbose = false, max_iter = 10)
        τ2, = compute_τII(c_pl, εII, args; verbose = false)
        @test τ1 == τ2

        args = merge(args, (τII = τ1,))
        F_check = compute_yieldfunction(c_pl.elements[3], args)
        @test abs(F_check) < 1.0e-12
    end

end
