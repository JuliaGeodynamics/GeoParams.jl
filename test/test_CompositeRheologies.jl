using Test
using GeoParams, ForwardDiff

@testset "CompositeRheologies" begin

    # Diffusion & dislocation creep in series
    v1 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
    v2 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
    t1 = (v1, v2)
    args = (T=900.0, d=100e-6)
    εII = 1e-15

    
    τII = compute_τII(t1, εII, args)
    @test τII ≈ 8.892678156850036e8
    η = τII / (2εII)

    # Combined creeplaws
    c1 = CompositeRheology(v1,v2)
    τII = compute_τII(c1, εII, args) # Same computation as 
    @test τII ≈ 8.892678156850036e8

    # Same with arrays:    
    εII_array = ones(10) * 1e-5
    τII_array = similar(εII_array)
    compute_τII!(τII_array, t1, εII_array, args)
    @test τII_array[1] ≈ 1.918028581543394e12

    # Add elasticity in the mix
    e1 = ConstantElasticity()
    t2 = (v1, v2, e1)
    dt_maxwell = η / NumValue(e1.G)
    args = (T=900.0, d=100e-6, τII_old=0.0, dt=dt_maxwell / 20)
    εII = 1e-15

    τII_vec = [0.0]
    t = [0.0]
    for i in 1:100
        τII = compute_τII(t2, εII, args)
        args = merge(args, (; τII_old=τII))
        push!(τII_vec, τII)
        push!(t, i * args.dt)
    end
    SecYear = 3600 * 24 * 365.25
    @test sum(τII_vec) ≈ 7.840307351918251e10 # this was the original number
   # SecYear = 3600 * 24 * 365.25
  
    # put the different rheological elements in a composite rheology structure
    v1 = LinearViscous()
    v2 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
    v3 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
    e1 = ConstantElasticity()
    p1 = DruckerPrager()

    # Check that constructing them works
    c2 = CompositeRheology( (v1, v2, v3, e1, Parallel(p1, v1, v2),v2, Parallel(p1, v1), v2,v3 ))   
    @test isa(c2.elements[1], AbstractCreepLaw)

    c3 = CompositeRheology( (v1, v2, v3, e1, Parallel(p1, e1, Parallel( (v1, v2), v3) ), v2,v3) )   
    @test isa(c3.elements[3], AbstractCreepLaw)
    
    c4 = Parallel((v2,v3,e1, Parallel(v2,v3),v2, Parallel(v2,(v3, v2))),v3,(e1,p1))
    @test isa(c4.elements[2], AbstractCreepLaw)

    # Perform computations for different complexites

    # Linear viscous
    v1  = (LinearViscous(η=1e21Pas), )
    c5  =  CompositeRheology(v1)

    args = (T=900.0, d=100e-6, τII_old=0.0)
    εII = 1e-15
    τII   = compute_τII(v1, εII, args)
    τII_2 = compute_τII(c5, εII, args)
    @test τII ≈ τII_2 ≈ 2*1e21*1e-15

    εII   = compute_εII(v1, τII, args)
    @test εII ≈ 1e-15


    # Linear viscoelastic 0D rheology test
    η  =  1e21;
    G  =  1e10;
    η,G  =  10, 1;
    
    t_M=  η/G
    εII = 1.;
    args=(;)
    c6  =  CompositeRheology((LinearViscous(η=η*Pa*s),ConstantElasticity(G=G*Pa)))
    t_vec, τ_vec =   time_τII_0D(c6, εII, args; t=(0.,t_M*4), nt=100, verbose=false)

    analytical_sol = @. 2.0*η*(1.0-exp(-t_vec/t_M))*εII
    
    err = sum(abs.(τ_vec .- analytical_sol))/length(t_vec)
    @test err ≈   0.0900844333898483


    # Case with nonlinear viscosity
    v1 = LinearViscous(η=1e23Pas)
    c6 =  CompositeRheology(v1,v3,ConstantElasticity())
    c7 =  CompositeRheology(v1,v3)
    y = (c6,c7)

    args = (T=1100.0, d=100e-6, τII_old=0.0)
    εII = 1e-15
    SecYear = 3600*24*365
    t_vec, τ_vec =   time_τII_0D(c6, εII, args; t=(0.,0.01*1e6*SecYear), nt=100, verbose=false)
    @test sum(τ_vec) ≈ 4.41429358183189e8


    # Parallel element
    pa1 = Parallel(LinearViscous(η=1e23Pas), LinearViscous(η=5e22Pas))  # put elements in parallel
    τII  = compute_τII(pa1, εII, args)
    τII_check = 2*εII*( NumValue(pa1.elements[1].η) +  NumValue(pa1.elements[2].η))
    @test τII ≈ τII_check

    # Parallel element with one element that has two components 
    c8    = (LinearViscous(η=1e23Pas), v2);
    v1    = LinearViscous(η=5e22Pas)
    pa2   = Parallel(v1, v2)  # put elements in parallel
    τII_1 = compute_τII(v1, εII, args)
    τII_2 = compute_τII(v2, εII, args)
    
    τII   = compute_τII(pa2, εII, args)       # this allocates!
    @test τII == (τII_1+τII_2)

    # Local iterations for a Parallel object for a given stress
    pa3 = Parallel(LinearViscous(η=1e23Pas), LinearViscous(η=5e22Pas))    # parallel object with 1 component/level
    τII = 1e6
    εII = compute_εII(pa3, τII, args, verbose=false)       # the strainrate of the parallel object (constant for all elements)
    τII_check = compute_τII(pa3, εII, args)
    @test τII == τII_check


    pa4 = Parallel( (LinearViscous(η=1e23Pas), v2), LinearViscous(η=5e22Pas))    # parallel object with 2 
    pa5 = Parallel( (LinearViscous(η=1e23Pas), v2), (LinearViscous(η=5e22Pas), Parallel(v1,v2)))    # parallel object with 2 entries but another parallel object 

    # test stress circuit implementations
    τII = GeoParams.MaterialParameters.ConstitutiveRelationships.stress_circuit(pa3, εII, args)
    @test τII == 1e6

    τII1 = GeoParams.MaterialParameters.ConstitutiveRelationships.stress_circuit(pa4, εII, args)
    @test τII1 ≈ 345408.7183763343

    # Check derivatives with FD
    pa6 = Parallel(v1,v2)
    Δτ  = τII*1e-6;
    ε0  = compute_εII(pa6, τII, args, verbose=false)
    ε1  = compute_εII(pa6, τII + Δτ, args, verbose=false)
    dε_dτ_FD = (ε1-ε0)/Δτ
    dεII_dτII(pa6, τII, args)

    # 
    ε = 1e-15
    Δε = 1e-6*ε
    τ0 = compute_τII(pa6, ε, args, verbose=false)
    τ1 = compute_τII(pa6, ε+Δε, args, verbose=false)
    dτ_dε_FD = (τ1-τ0)/Δε

    dτII_dεII(pa6, ε, args)
    
    # use AD to check the individual derivatives of the separate creeplaws 
    pa7 = Parallel(v1,v2,v3)
    for i=1:length(pa7.elements)
        f(ε) = compute_τII(pa7.elements[i], ε, args)
        der_AD = ForwardDiff.derivative(f,1e-15)
        @test der_AD ≈ dτII_dεII(pa7.elements[i], 1e-15, args)

        g(τ) = compute_εII(pa7.elements[i], τ, args)
        der_AD = ForwardDiff.derivative(g, 1e6)
        @test der_AD ≈ dεII_dτII(pa7.elements[i], 1e6, args)
    end

    # test derivatives of the combined creep laws
    f2(ε) = compute_τII(pa7, ε, args)
    der_AD2 = ForwardDiff.derivative(f2,1e-15)
    @test der_AD2 ≈ dτII_dεII(pa7, 1e-15, args)   # correct; we need this derivative for the jacobian


    c9 = CompositeRheology(v2,v3,Parallel(v1,v2))
    der = GeoParams.MaterialParameters.ConstitutiveRelationships.dεII_dτII_elements(c9,1e6,args)         # this should compute dεII/dτII for serial elements that are NOT Parallel blocks
    εII = GeoParams.MaterialParameters.ConstitutiveRelationships.compute_εII_elements(c9,1e6,args)       # sums εII for serial elements that are NOT Parallel blocks

    # Check that the sum & derivatives only apply to non-parallel elements
    ε,dε = 0.0, 0.0
    for i=1:length(c9.elements)
        if !isa(c9.elements[i], Parallel)
            ε += compute_εII(c9.elements[i],1e6,args)
            dε += dεII_dτII(c9.elements[i],1e6,args)
        end
    end
    @test ε ≈ εII
    @test dε ≈ der

    εII= 3e-15
    v1 = LinearViscous()
    v2 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
    v3 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
    c  = CompositeRheology(v2,v1,Parallel(v2,v1))
    τ  = local_iterations_εII(c, εII, args, verbose=true)

    @test τ ≈ 569147.233065495

end

