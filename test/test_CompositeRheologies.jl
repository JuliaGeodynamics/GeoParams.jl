using Test
using GeoParams #, ForwardDiff

@testset "CompositeRheologies" begin

    # Diffusion & dislocation creep in series
    pp = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
    pp1 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
    v = (pp, pp1)
    args = (T=900.0, d=100e-6)
    εII = 1e-15

    
    τII = compute_τII(v, εII, args)
    @test τII ≈ 8.892678156850036e8
    η = τII / (2εII)

    # Combined creeplaws
    v1 = CompositeRheology(pp,pp1)
    τII = compute_τII(v, εII, args) # Same computation as 
    @test τII ≈ 8.892678156850036e8

    # Same with arrays:    
    εII_array = ones(10) * 1e-5
    τII_array = similar(εII_array)
    compute_τII!(τII_array, v, εII_array, args)
    @test τII_array[1] ≈ 1.918028581543394e12

    # compute strainrate given stress  [not working yet!]
    # εII1  = compute_εII(v,τII, args) 

    # Add elasticity in the mix
    el = ConstantElasticity()
    v = (pp, pp1, el)
    dt_maxwell = η / NumValue(el.G)
    args = (T=900.0, d=100e-6, τII_old=0.0, dt=dt_maxwell / 20)
    εII = 1e-15

    τII_vec = [0.0]
    t = [0.0]
    for i in 1:100
        τII = compute_τII(v, εII, args)
        args = merge(args, (; τII_old=τII))
        push!(τII_vec, τII)
        push!(t, i * args.dt)
    end
    SecYear = 3600 * 24 * 365.25
    @test sum(τII_vec) ≈ 7.840307351918251e10 # this was the original number
   # SecYear = 3600 * 24 * 365.25
  
    # put the different rheological elements in a composite rheology structure
    pp0 = LinearViscous()
    pp1 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
    pp2 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
    
    pp3 = ConstantElasticity()
    pp4 = DruckerPrager()

    # Check that constructing them works
    a = CompositeRheology( (pp0, pp1, pp2, pp3, Parallel(pp4, pp0, pp1),pp1, Parallel(pp4, pp0), pp1,pp2 ))   
    @test isa(a.elements[1], AbstractCreepLaw)

    b = CompositeRheology( (pp0, pp1, pp2, pp3, Parallel(pp4, pp3, Parallel( (pp1, pp0), pp2) ), pp1,pp2) )   
    @test isa(b.elements[3], AbstractCreepLaw)
    
    a=Parallel((pp1,pp2,pp3, Parallel(pp1,pp2),pp1, Parallel(pp1,(pp2, pp1))),pp2,(pp3,pp4))

    # Perform computations for different complexites

    # Linear viscous
    v1  = (LinearViscous(η=1e21Pas), )
    v2  =  CompositeRheology(v1)

    args = (T=900.0, d=100e-6, τII_old=0.0)
    εII = 1e-15
    τII   = compute_τII(v1, εII, args)
    τII_2 = compute_τII(v2, εII, args)
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
    x  =  CompositeRheology((LinearViscous(η=η*Pa*s),ConstantElasticity(G=G*Pa)))
    t_vec, τ_vec =   time_τII_0D(x, εII, args; t=(0.,t_M*4), nt=100, verbose=false)

    analytical_sol = @. 2.0*η*(1.0-exp(-t_vec/t_M))*εII
    
    err = sum(abs.(τ_vec .- analytical_sol))/length(t_vec)
    @test err ≈   0.0900844333898483


    # Case with nonlinear viscosity
    pp0= LinearViscous(η=1e23Pas)
    x  =  CompositeRheology(pp0,pp2,ConstantElasticity())
    x1 =  CompositeRheology(pp0,pp2)
    y = (x,x1)

    args = (T=1100.0, d=100e-6, τII_old=0.0)
    εII = 1e-15
    SecYear = 3600*24*365
    t_vec, τ_vec =   time_τII_0D(x, εII, args; t=(0.,0.01*1e6*SecYear), nt=100, verbose=false)
    @test sum(τ_vec) ≈ 4.41429358183189e8


    # Parallel element
    x = Parallel(LinearViscous(η=1e23Pas), LinearViscous(η=5e22Pas))  # put elements in parallel
    τII   = compute_τII(x, εII, args)
    τII_check = 2*εII*( NumValue(x.elements[1].η) +  NumValue(x.elements[2].η))
    @test τII ≈ τII_check

    # Parallel element with one element that has two components 
    v1    = (LinearViscous(η=1e23Pas), pp1);
    v2    = LinearViscous(η=5e22Pas)
    x     = Parallel(v1, v2)  # put elements in parallel
    τII_1 = compute_τII(v1, εII, args)
    τII_2 = compute_τII(v2, εII, args)
    
    τII   = compute_τII(x, εII, args)       # this allocates!
    @test τII == (τII_1+τII_2)

    # Local iterations for a Parallel object for a given stress
    x = Parallel(LinearViscous(η=1e23Pas), LinearViscous(η=5e22Pas))    # parallel object with 1 component/level
    τII = 1e6
    εII = compute_εII(x, τII, args, verbose=true)       # the strainrate of the parallel object (constant for all elements)
    τII_check = compute_τII(x, εII, args)
    @test τII == τII_check


    x1 = Parallel( (LinearViscous(η=1e23Pas), pp1), LinearViscous(η=5e22Pas))    # parallel object with 2 
    x2 = Parallel( (LinearViscous(η=1e23Pas), pp1), (LinearViscous(η=5e22Pas), Parallel(pp0,pp1)))    # parallel object with 2 entries but another parallel object 

    

    # test stress circuit implementations
    τII = GeoParams.MaterialParameters.ConstitutiveRelationships.stress_circuit(x, εII, args)
    @test τII == 1e6

    τII1 = GeoParams.MaterialParameters.ConstitutiveRelationships.stress_circuit(x1, εII, args)
    @test τII1 ≈ 345408.7183763343

    # Check derivatives with FD
    xt = Parallel(pp0,pp1)
    Δτ = τII*1e-6;
    ε0 = compute_εII(xt, τII, args, verbose=false)
    ε1 = compute_εII(xt, τII + Δτ, args, verbose=false)
    dε_dτ_FD = (ε1-ε0)/Δτ
    dεII_dτII(x1, τII, args)

    # 
    ε = 1e-15
    Δε = 1e-6*ε
    τ0 = compute_τII(xt, ε, args, verbose=false)
    τ1 = compute_τII(xt, ε+Δε, args, verbose=false)
    dτ_dε_FD = (τ1-τ0)/Δε

    dτII_dεII(x1, ε, args)
    
    #=
    # use AD to test derivatives
    xt = Parallel(pp0,pp1, pp2)
    for i=1:length(xt.elements)
        @show i
        f(ε) = compute_τII(xt.elements[i], ε, args)
        der_AD = ForwardDiff.derivative(f,1e-15)
        @test der_AD ≈ dτII_dεII(xt.elements[i], 1e-15, args)

        g(τ) = compute_εII(xt.elements[i], τ, args)
        der_AD = ForwardDiff.derivative(g,1e-15)
        @test der_AD ≈ dεII_dτII(xt.elements[i], 1e6, args)

    end

=#

    εII = compute_εII(x, τII, args, verbose=true) 

end

