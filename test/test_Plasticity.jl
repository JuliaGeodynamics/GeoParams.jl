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
    num_alloc = @allocated compute_yieldfunction!(F, MatParam, PhaseRatio, args)
    @test maximum(F[1, 1, :]) ≈ 839745.962155614
    # @test num_alloc <= 32

    # Test plastic potential derivatives
    ## 2D
    τij = (1.0, 2.0, 3.0)
    fxx(τij) = 0.5 * τij[1]/second_invariant(τij)
    fyy(τij) = 0.5 * τij[2]/second_invariant(τij)
    fxy(τij) = τij[3]/second_invariant(τij)
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
    gxx(τij) = 0.5 * τij[1]/second_invariant(τij)
    gyy(τij) = 0.5 * τij[2]/second_invariant(τij)
    gzz(τij) = 0.5 * τij[3]/second_invariant(τij)
    gyz(τij) = τij[4]/second_invariant(τij)
    gxz(τij) = τij[5]/second_invariant(τij)
    gxy(τij) = τij[6]/second_invariant(τij)
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


    # try elastic element
    e1 = ConstantElasticity(G=10Pa)           # elasticity
    pl1 = DruckerPrager(ϕ=0, C=15)                # plasticity

    
    τ_new  = compute_τII(e1,1, args)
    
    args = (τII_old=5.0, dt=1.0, τII=5.0)
    F_old = compute_yieldfunction(pl1, args)

    args = (τII_old=5.0, dt=1.0, τII=τ_new)
    F_new = compute_yieldfunction(pl1, args)

    ε_pl =  compute_εII(pl1, 1.0, τ_new, args)  # plastic strian rate
    
    
    J[1,1] =  GeoParams.MaterialParameters.ConstitutiveRelationships.dεII_dτII_elements(e1,τ,args);

    # ===
    # Local Iterations
    τ_initial = 5.0
    εII_total = 1.0

    x = [τ_new; 0.0]
    
    verbose=true
    tol=1e-6
    iter = 0
    ϵ = 2 * tol
    τII_prev = τ_initial
    τ_parallel = 0
    max_iter = 10
    J = zeros(2,2)
    r = zeros(2)
    while (ϵ > tol) && (iter < max_iter)
        iter += 1

        τ   = x[1]
        λ̇   = x[2]

        args = merge(args, (τII=τ,))
        
        # Update part of jacobian related to serial elements
        c = e1
        r[1]   = εII_total - GeoParams.MaterialParameters.ConstitutiveRelationships.compute_εII(c,τ,args)
        J[1,1] = GeoParams.MaterialParameters.ConstitutiveRelationships.dεII_dτII(c,x[1],args);
        
        r[1]   -=  λ̇*∂Q∂τII(pl1, τ)         # add plastis strainrate

        # Add contributions from || elements
        F = compute_yieldfunction(pl1,args);
#        J[2,2] = 1

        J[1,2] = ∂Q∂τII(pl1, τ)     
        J[2,1] = ∂F∂τII(pl1, τ)    
        #J[2,2] = F

        r[2] =  -F 
#        r[2] = 0.0
        
        #fill_J_parallel!(J, r, x, c, τ, args)
    
        # update solution
        dx  = J\r 
        x .+= dx   

        ϵ    = sum(abs.(dx)./(abs.(x .+ 1e-9)))
        verbose && println(" iter $(iter) $ϵ")
        @show x
    end
    verbose && println("---")
    # ===
    

end