using Test
using GeoParams

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
    @test isa(a.rheology_chain[1], AbstractCreepLaw)

    b = CompositeRheology( (pp0, pp1, pp2, pp3, Parallel(pp4, pp3, Parallel( (pp1, pp0), pp2) ), pp1,pp2) )   
    @test isa(b.rheology_chain[3], AbstractCreepLaw)
    
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




end
