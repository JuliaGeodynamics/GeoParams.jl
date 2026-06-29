using Test
using GeoParams
using StaticArrays

@testset "Plasticity.jl" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 10km)
    @testset "DruckerPrager" begin
        # DruckerPrager ---------
        p = DruckerPrager()
        info = param_info(p)
        @test info.Equation === L"$F = \\tau_{II} - \\cos(ϕ)C - \\sin(ϕ)(P-P_f); Q=\\tau_{II} - \\sin(Ψ)(P-P_f)$"
        @test sprint(show, p) == "Drucker-Prager plasticity with: C = 1.0e7 Pa, ϕ = 30.0ᵒ, Ψ = 0.0ᵒ"
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
        ## 2D — constrained gradient: τzz = -τxx - τyy is a function of stored components
        τij = (1.0, 2.0, 3.0)
        fxx(τij) = (τij[1] + 0.5 * τij[2]) / second_invariant(τij)
        fyy(τij) = (0.5 * τij[1] + τij[2]) / second_invariant(τij)
        fxy(τij) = τij[3] / second_invariant(τij)
        solution2D = [fxx(τij), fyy(τij), fxy(τij)]

        # using StaticArrays (SVector{3} dispatch of ∂Q∂τ)
        out1 = ∂Q∂τ(p, @SVector [1.0, 2.0, 3.0])
        @test out1 isa SVector && collect(out1) ≈ solution2D

        # using tuples
        τij_tuple = (1.0, 2.0, 3.0)
        out2 = ∂Q∂τ(p, τij_tuple)
        @test out2 == Tuple(solution2D)
        @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(p, τij_tuple)

        # using AD — must agree with the analytical constrained gradient
        Q = second_invariant # where second_invariant is a function
        ad2 = ∂Q∂τ(Q, τij_tuple)
        @test ad2 == (0.5, 0.625, 0.75)
        @test collect(out1) ≈ collect(ad2)

        ## 3D
        τij = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
        gxx(τij) = 0.5 * τij[1] / second_invariant(τij)
        gyy(τij) = 0.5 * τij[2] / second_invariant(τij)
        gzz(τij) = 0.5 * τij[3] / second_invariant(τij)
        gyz(τij) = τij[4] / second_invariant(τij)
        gxz(τij) = τij[5] / second_invariant(τij)
        gxy(τij) = τij[6] / second_invariant(τij)
        solution3D = [gxx(τij), gyy(τij), gzz(τij), gyz(τij), gxz(τij), gxy(τij)]

        # using StaticArrays (SVector{6} dispatch of ∂Q∂τ)
        out3 = ∂Q∂τ(p, @SVector [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        @test out3 isa SVector && collect(out3) ≈ solution3D

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
        @test info.Equation === L"$F = \\tau_{II} - \\cos(ϕ)C - \\sin(ϕ)(P-P_f) - 2η_vpε̇II_pl ; Q=\\tau_{II} - \\sin(Ψ)(P-P_f)$"
        @test sprint(show, p) == "Regularized Drucker-Prager plasticity with: C = 1.0e7 Pa, ϕ = 30.0ᵒ, Ψ = 0.0ᵒ, η_vp=1.0e20 Pa s"
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
        ## 2D — constrained gradient: τzz = -τxx - τyy is a function of stored components
        τij = (1.0, 2.0, 3.0)
        fxx(τij) = (τij[1] + 0.5 * τij[2]) / second_invariant(τij)
        fyy(τij) = (0.5 * τij[1] + τij[2]) / second_invariant(τij)
        fxy(τij) = τij[3] / second_invariant(τij)
        solution2D = [fxx(τij), fyy(τij), fxy(τij)]

        # using StaticArrays (SVector{3} dispatch of ∂Q∂τ)
        out1 = ∂Q∂τ(p, @SVector [1.0, 2.0, 3.0])
        @test out1 isa SVector && collect(out1) ≈ solution2D

        # using tuples
        τij_tuple = (1.0, 2.0, 3.0)
        out2 = ∂Q∂τ(p, τij_tuple)
        @test out2 == Tuple(solution2D)
        @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(p, τij_tuple)

        # using AD — must agree with the analytical constrained gradient
        Q = second_invariant # where second_invariant is a function
        ad2 = ∂Q∂τ(Q, τij_tuple)
        @test ad2 == (0.5, 0.625, 0.75)
        @test collect(out1) ≈ collect(ad2)

        ## 3D
        τij = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
        gxx(τij) = 0.5 * τij[1] / second_invariant(τij)
        gyy(τij) = 0.5 * τij[2] / second_invariant(τij)
        gzz(τij) = 0.5 * τij[3] / second_invariant(τij)
        gyz(τij) = τij[4] / second_invariant(τij)
        gxz(τij) = τij[5] / second_invariant(τij)
        gxy(τij) = τij[6] / second_invariant(τij)
        solution3D = [gxx(τij), gyy(τij), gzz(τij), gyz(τij), gxz(τij), gxy(τij)]

        # using StaticArrays (SVector{6} dispatch of ∂Q∂τ)
        out3 = ∂Q∂τ(p, @SVector [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        @test out3 isa SVector && collect(out3) ≈ solution3D

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
        @test info.Equation === L"$F = \tau_{II} - kP - c \;\;\mathrm{or}\;\; a(\sqrt{\tau_{II}^2 + (P-p_y)^2} - R_y),\; Q = \tau_{II} - k_q P - \mathrm{const} \;\mathrm{or}\;\; b(\sqrt{\tau_{II}^2 + (P-p_q)^2} - R_f)$"
        @test sprint(show, p) == "DruckerPragerCap(ϕ=30.0, Ψ=0.0, C=1.0e7 Pa, η_vp=1.0e20 Pa s, pT=-100000.0 Pa"
        @test isbits(p)

        # mixed-units constructor: C/pT unit promotion warning branches
        cp_mix1 = DruckerPragerCap(; C = 10.0e6, pT = -1.0e5Pa)   # C unitless, pT united
        @test cp_mix1 isa DruckerPragerCap
        cp_mix2 = DruckerPragerCap(; C = 10.0e6Pa, pT = -1.0e5)   # C united, pT unitless
        @test cp_mix2 isa DruckerPragerCap
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
        ## 2D — constrained gradient: τzz = -τxx - τyy is a function of stored components
        τij = (1.0, 2.0, 3.0)
        τII_test = second_invariant(τij)
        dQτII = ∂Q∂τII(p, τII_test)
        fxx(τij) = dQτII * (τij[1] + 0.5 * τij[2]) / second_invariant(τij)
        fyy(τij) = dQτII * (0.5 * τij[1] + τij[2]) / second_invariant(τij)
        fxy(τij) = 2 * dQτII * τij[3] / second_invariant(τij)
        solution2D = [fxx(τij), fyy(τij), fxy(τij)]

        # using tuples
        τij_tuple = (1.0, 2.0, 3.0)
        out2 = ∂Q∂τ(p, τij_tuple)
        @test out2 == Tuple(solution2D)
        @test compute_plasticpotentialDerivative(p, τij_tuple) == ∂Q∂τ(p, τij_tuple)

        # using AD — must agree with analytical constrained gradient scaled by Aτ
        Q = second_invariant # where second_invariant is a function
        ad2 = ∂Q∂τ(Q, τij_tuple)
        @test ad2 == (0.5, 0.625, 0.75)
        @test all(isapprox.(dQτII .* collect(ad2), [out2[1], out2[2], out2[3] / 2]; rtol = 1.0e-5))

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

    @testset "no functor ambiguities" begin
        # The specialized yield-function functors must each be a strict subtype of the
        # generic one (free softening type vars bounded by AbstractSoftening), otherwise
        # `compute_yieldfunction` throws a MethodError on real calls. Guard against the
        # regression where the specializations were either dead or mutually ambiguous.
        mod = GeoParams.MaterialParameters.ConstitutiveRelationships
        amb = Test.detect_ambiguities(mod; recursive = false)
        plast = filter(
            m -> occursin("Plasticity", string(m[1].file)) || occursin("Plasticity", string(m[2].file)),
            amb,
        )
        @test isempty(plast)

        # every softening combination must dispatch and apply softening at EII > 0
        sC = LinearSoftening(0.0e0, 1.0e7, 0.0e0, 1.0e0)
        sϕ = LinearSoftening(15.0e0, 30.0e0, 0.0e0, 1.0e0)
        sΨ = LinearSoftening(0.0e0, 20.0e0, 0.0e0, 1.0e0)
        for p in (
                DruckerPrager(; softening_C = sC),
                DruckerPrager(; softening_ϕ = sϕ),
                DruckerPrager(; softening_ϕ = sϕ, softening_C = sC),
                DruckerPrager_regularised(; softening_C = sC),
                DruckerPrager_regularised(; softening_ϕ = sϕ),
                DruckerPragerCap(; Ψ = 20.0, softening_C = sC),
                DruckerPragerCap(; Ψ = 20.0, softening_ϕ = sϕ, softening_Ψ = sΨ),
                DruckerPragerCap(; Ψ = 20.0, softening_ϕ = sϕ, softening_C = sC, softening_Ψ = sΨ),
            )
            F0 = compute_yieldfunction(p; P = 1.0e6, τII = 2.0e7, EII = 0.0)
            F1 = compute_yieldfunction(p; P = 1.0e6, τII = 2.0e7, EII = 1.0)
            @test isfinite(F0) && isfinite(F1)
        end
    end

    @testset "plastic_strain / lambda" begin
        mod = GeoParams.MaterialParameters.ConstitutiveRelationships
        p = DruckerPrager(; ϕ = 30.0, Ψ = 10.0, C = 1.0e7Pa)
        εvp = mod.plastic_strain(p, (1.0, 1.0, 1.0), 1.0e-15)
        @test isfinite(εvp) && εvp > 0
        # accumulating form: plastic_strain(εvp, p, τij, λ̇, dt)
        εvp2 = mod.plastic_strain(0.0, p, (1.0, 1.0, 1.0), 1.0e-15, 1.0e3)
        @test isfinite(εvp2) && εvp2 > 0
        @test isfinite(mod.lambda(1.0e6, p, 1.0e20, 1.0e19))
        @test isfinite(mod.lambda(1.0e6, p, 1.0e20, 1.0e19; K = 2.0e10, dt = 1.0e3, h = 1.0e5, τij = (1.0e6, 1.0e6, 1.0e6)))

        # ∂Q∂τ with a NamedTuple argument (dispatches on args.τij)
        @test ∂Q∂τ(p, (; τij = (1.0, 2.0, 3.0))) == ∂Q∂τ(p, (1.0, 2.0, 3.0))
    end

    @testset "dimensionalize round-trip" begin
        CharDim = GEO_units(; viscosity = 1.0e19, length = 10km)
        @test dimensionalize(NoSoftening(), CharDim) === NoSoftening()
        for p in (
                DruckerPrager(; C = 1.0e7Pa, ϕ = 30.0),
                DruckerPrager_regularised(; C = 1.0e7Pa, η_vp = 1.0e20Pa * s),
                DruckerPragerCap(; C = 1.0e7Pa, pT = -5.0e5Pa, η_vp = 1.0e20Pa * s, Ψ = 10.0),
            )
            pd = dimensionalize(nondimensionalize(p, CharDim), CharDim)
            @test pd.C.val ≈ p.C.val
        end
    end

    @testset "DruckerPragerCap flow potential & derivatives" begin
        p = DruckerPragerCap()
        P, τII = 1.0e6, 20.0e6
        # regression values for the cap flow potential and its derivatives
        @test compute_flowpotential(p; P = P, τII = τII) ≈ 1.3360009139380395e7 rtol = 1.0e-9
        @test compute_flowpotential(p, (; P = P, τII = τII)) ≈ 1.3360009139380395e7 rtol = 1.0e-9
        @test ∂Q∂P(p, P; τII = τII) ≈ -0.31343077600453645 rtol = 1.0e-9
        @test ∂Q∂τII(p, P; τII = τII) ≈ 0.06521675174605653 rtol = 1.0e-9
        @test ∂F∂P(p, P; τII = τII) ≈ -0.6037375326226843 rtol = 1.0e-9
        @test ∂F∂τII(p, P; τII = τII) ≈ 0.08062094643606207 rtol = 1.0e-9

        # array variant matches the scalar result
        Pa = fill(P, 5)
        τa = fill(τII, 5)
        Q = similar(Pa)
        compute_flowpotential!(Q, p; P = Pa, τII = τa)
        @test Q[1] ≈ 1.3360009139380395e7 rtol = 1.0e-9

        # softening functor variants (EII > 0 path)
        sC = LinearSoftening(0.0e0, 1.0e7, 0.0e0, 1.0e0)
        sϕ = LinearSoftening(15.0e0, 30.0e0, 0.0e0, 1.0e0)
        sΨ = LinearSoftening(0.0e0, 20.0e0, 0.0e0, 1.0e0)
        for ps in (
                DruckerPragerCap(; Ψ = 20.0, softening_C = sC),
                DruckerPragerCap(; Ψ = 20.0, softening_ϕ = sϕ),
                DruckerPragerCap(; Ψ = 20.0, softening_Ψ = sΨ),
                DruckerPragerCap(; Ψ = 20.0, softening_ϕ = sϕ, softening_C = sC),
                DruckerPragerCap(; Ψ = 20.0, softening_ϕ = sϕ, softening_Ψ = sΨ),
                DruckerPragerCap(; Ψ = 20.0, softening_C = sC, softening_Ψ = sΨ),
                DruckerPragerCap(; Ψ = 20.0, softening_ϕ = sϕ, softening_C = sC, softening_Ψ = sΨ),
            )
            @test compute_yieldfunction(ps; P = P, τII = τII, EII = 0.5) isa Number
            @test compute_flowpotential(ps; P = P, τII = τII, EII = 0.5) isa Number
        end
    end

    @testset "DruckerPrager generic (args, kwargs) interfaces" begin
        p = DruckerPrager(; C = 1.0e6)
        a = (; P = 1.0e6, τII = 2.0e7)
        # functor + generic positional-args forms (regression values)
        @test p(a) ≈ 1.863397459621556e7 rtol = 1.0e-9
        @test compute_yieldfunction(p, a) ≈ 1.863397459621556e7 rtol = 1.0e-9
        @test ∂Q∂P(p, a, (;)) == -0.0   # associated DP with Ψ = 0 -> ∂Q/∂P = 0

        # ∂Q∂τ on a user function (AD gradient over a Vector)
        Q(τ) = sqrt(sum(abs2, τ))
        @test ∂Q∂τ(Q, [1.0, 2.0, 2.0]) ≈ [1 / 3, 2 / 3, 2 / 3] rtol = 1.0e-9
        # 3-arg (args, kwargs) form + the MaterialParams plasticity dispatch
        τ = @SVector [1.0, 2.0, 3.0]
        @test ∂Q∂τ(p, τ, (;)) ≈ [0.5, 0.625, 0.75] rtol = 1.0e-9
        mp = SetMaterialParams(; Phase = 1, Plasticity = DruckerPrager(; C = 1.0e6))
        @test ∂Q∂τ(mp, τ) ≈ [0.5, 0.625, 0.75] rtol = 1.0e-9
    end

end
