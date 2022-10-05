using Test
using GeoParams, ForwardDiff

@testset "CompositeRheologies" begin

    # Define a range of rheological components
    v1 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
    v2 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
    v3 = LinearViscous()
    v4 = LinearViscous(η=1e22Pa*s)
    e1 = ConstantElasticity()           # elasticity
    pl1= DruckerPrager()                # plasticity

    # Parallel elements
    p1 = Parallel(v3,v4)                # linear elements
    p2 = Parallel(v1,v2)                # includes nonlinear viscous elements
    p3 = Parallel(v1,v2,v3)             # includes nonlinear viscous elements

    # CompositeRheologies
    c1 = CompositeRheology(v1,v2)
    c2 = CompositeRheology(v3,v4)       # two linear rheologies
    c3 = CompositeRheology(v1,v2, e1)   # with elasticity
    c4 = CompositeRheology(v1,v3, p1)   # with linear || element
    c5 = CompositeRheology(v1,v4, p2)   # with nonlinear || element
    c6 = CompositeRheology(v1,v4,p1,p2) # with 2 || elements
    c7 = CompositeRheology(v3,e1,pl1)   # with plastic element
    p4 = Parallel(c3,v3)                # Parallel element with composite one as well    

    # Check that we can construct complicated rheological elements
    c = CompositeRheology( (v1, v2, v3, e1, Parallel(p1, v1, v2),v2, Parallel(p1, v1), v2,v3 ))   
    @test isa(c.elements[1], AbstractCreepLaw)

    c = CompositeRheology( (v1, v2, v3, e1, Parallel(p1, e1, Parallel( (v1, v2), v3) ), v2,v3) )   
    @test isa(c.elements[3], AbstractCreepLaw)
    
    c = Parallel((v2,v3,e1, Parallel(v2,v3),v2, Parallel(v2,(v3, v2))),v3,(e1,p1))
    @test isa(c.elements[2], AbstractCreepLaw)
    
    args = (T=900.0, d=100e-6, τII_old=1e6, dt=1e8)
    εII, τII = 2e-15, 2e6

    # Check derivatives 
    vec1 = [c1 c2 c3 c4 c5 p1 p2 p3] 
    for v in vec1   # problem with c3 (elasticity) and c4 (that has || elements) 
        τII = 1e9
        Δτ  = τII*1e-6;
        ε0  = compute_εII(v, τII, args, verbose=false);
        ε1  = compute_εII(v, τII + Δτ, args, verbose=false);
        dε_dτ_FD = (ε1-ε0)/Δτ
        dε_dτ_AD = dεII_dτII_AD(v, τII, args)
        dε_dτ    = dεII_dτII(v, τII, args)       
        @test dε_dτ ≈ dε_dτ_AD  
        @test dε_dτ ≈ dε_dτ_FD  rtol=1e-6 
        @test dε_dτ_AD ≈ dε_dτ_FD  rtol=1e-6 
        
        εII = 1e-12 #2e-18
        Δε  = εII*1e-6;
        τ0  = compute_τII(v, εII, args, verbose=false);
        τ1  = compute_τII(v, εII + Δε, args, verbose=false);
        dτ_dε_FD = (τ1-τ0)/Δε
        dτ_dε_AD = dτII_dεII_AD(v, εII, args)
        dτ_dε    = dτII_dεII(v, εII, args)       
        @test dτ_dε ≈ dτ_dε_AD  
        @test dτ_dε ≈ dτ_dε_AD  rtol=1e-5 
        @test dτ_dε_AD ≈ dτ_dε_FD  rtol=1e-5 
        
    end


    # Check computational routines for strainrate & stress
    for v in vec1

        # Check computational routines if total strainrate is given
        τ    = compute_τII(v, εII, args)   
        τ_AD = compute_τII_AD(v, εII, args)   
        ε    = compute_εII(v, τ,   args)   
        ε_AD = compute_εII_AD(v, τ,   args)   
        
        @test ε ≈ εII ≈ ε_AD 
        @test τ ≈ τ_AD

        if isa(v,Parallel)
            # check; for || elements, ε is constant and τ is the sum
            τ_check = 0;
            for i=1:length(v.elements)
                τ_check += compute_τII(v.elements[i], εII, args)   
            end
            @test τ_check ≈ τ
        end
        
        # Check computational routines if stress is given
        ε    = compute_εII(v, τII, args)   
        ε_AD = compute_εII_AD(v, τII, args)   
        τ    = compute_τII(v, ε,   args)   
        τ_AD = compute_τII_AD(v, ε,   args)   
        
        @test τ ≈ τII ≈ τ_AD
        @test ε ≈ ε_AD
        
        if isa(v,CompositeRheology)
            # check; for serial elements, ε is constant and τ is the sum
            ε_check = 0;
            for i=1:length(v.elements)
                ε_check += compute_εII(v.elements[i], τII, args)   
            end
            @test ε_check ≈ ε
        end

        # Viscosity calculation routines
        η    =  computeViscosity_εII(v, εII, args)
        η_AD =  computeViscosity_εII_AD(v, εII, args)
        @test η ≈ η_AD


        args_dim = (T=900.0K, d=100e-6m, τII_old=1e6Pa)
        εII_dim, τII_dim = 2e-15/s, 2e6Pa

        # Do calculations with values that have units (currently only works for single creep laws)
#        τ_dim = compute_τII(v, εII_dim, args_dim)
#        ε_dim = compute_εII(v, τ_dim,   args_dim)
#        @test ε_dim ≈ εII_dim
#        η_dim = computeViscosity_εII(v,εII_dim,args_dim)


    end

    # CompositeRheology cases with parallel elements
    εII =  3e-15
    for v in [c4 c5 c6]    

        τ_AD = compute_τII_AD(v, εII, args)     # using AD
       
        # check result. For a parallel element we should satisfy the following equations:
        #   τ_parallel == τ_AD 
        #   sum(ε) = εII
        ε_vec = [compute_εII(v.elements[i], τ_AD, args) for i=1:length(v.elements)]
        @test sum(ε_vec) ≈ εII

        ε_parallel,τ_parallel = 0.0, 0.0
        for i=1:length(v.elements)
            # Check that the stress of each || element is the same as the total one
            if isa(v.elements[i], Parallel)
                ε_parallel = ε_vec[i]      # parallel strainrate
                τ_parallel = compute_τII_AD(v.elements[i], ε_parallel, args) 
                @test τ_parallel  ≈ τ_AD
            end
        end

        @test εII ≈ compute_εII(v,τ_AD,args)    # sum of ε

        # Check with analytical jacobian
        τ    = compute_τII(v, εII, args)        # using analytical jacobians (expanded for || elements)
        @test τ_AD ≈ τ

    end

    # Test Parallel rheology with a CompositeRheology branch
    εII =  3e-15
    for v in [p4]
        τ_AD = compute_τII_AD(v, εII, args)     
        τ    = compute_τII(v, εII, args)     

        τ_par = 0.
        for i=1:length(v.elements)
            τ_par  += compute_τII(v.elements[i], εII, args)  
        end

        @test τ_AD ≈ τ ≈ τ_par
        
    end

    # Composite cases with plasticity 
    for v in [c7]
        #τ_AD = compute_τII_AD(v, εII, args)     
        τ    = compute_τII(v, εII, args)   
        args = merge(args, (τII=τ,P=1e9))  
        F    = compute_yieldfunction(v.elements[3],args)

        # Check that F==0    
        F = 0.
        for i=1:length(v.elements)
            if isa(v.elements,AbstractPlasticity)
                F  += compute_τII(v.elements[i], args)  
            end
        end
        @test F ≈ 0.0

    end

    # 0D rheology tests
    η,G  =  10, 1;
    t_M  =  η/G
    εII  =  1.;
    args = (;)
    pl2  =  DruckerPrager(C=η, ϕ=0)                # plasticity

    c_lin =  CompositeRheology(LinearViscous(η=η*Pa*s),ConstantElasticity(G=G*Pa)) # linear VE
    t_vec, τ_vec    =   time_τII_0D(c_lin, εII, args; t=(0.,t_M*4), nt=100, verbose=true)
    analytical_sol  =   @. 2.0*η*(1.0-exp(-t_vec/t_M))*εII
    err             =   sum(abs.(τ_vec .- analytical_sol))/length(t_vec)
    @test err ≈ 0.0900844333898483

    # Plasticity calculation
    F = zero(τ_vec)
    P = 0.0
    for (i,τ) = enumerate(τ_vec)
        args = (P=P, τII=τ)
        F[i] = compute_yieldfunction(pl2,args)
    end

    # VEP calculation
    c_pl  =  CompositeRheology(LinearViscous(η=η*Pa*s),ConstantElasticity(G=G*Pa),pl2) # linear VEP
    t_vec, τ_vec    =   time_τII_0D(c_pl, εII, args; t=(0.,t_M*4), nt=100, verbose=false)
    

    # test non-dimensionalisation
    CharDim = GEO_units()
    for (ind,v) in enumerate(vec1)
        p_nd  = nondimensionalize(v,CharDim)
        p_dim = dimensionalize(p_nd,CharDim)

        # check
        for (i,v_local) in enumerate(v.elements)
            names = fieldnames(typeof(v_local))
            p_local = p_dim[i]
            for name in names
                d_dim = getfield(p_local, name)
                v_dim = getfield(v_local, name)
                if isa(d_dim,GeoUnit)
                    if !isnan(NumValue(v_dim))
                        @test NumValue(v_dim) ≈ NumValue(d_dim)
                    end
                end
            end
        end

    end

    # Specify a composite rheology in the MaterialParam struct 
    MatParam_nd = SetMaterialParams(Name="Viscous Matrix", Phase=2,
                Density   = ConstantDensity(),
                CompositeRheology   = c5, CharDim=CharDim)
    @test NumValue(MatParam_nd.CompositeRheology[1][2].η) ≈ 100

    MatParam   = SetMaterialParams(Name="Viscous Matrix", Phase=2,
                Density   = ConstantDensity(),
                CompositeRheology   = c5)
    
    @test NumValue(MatParam.CompositeRheology[1][2].η) ≈ 1e22

end
