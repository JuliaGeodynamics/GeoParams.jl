using Test
using GeoParams, ForwardDiff

@testset "CompositeRheologies" begin

    # Define a range of rheological components
    v1 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
    v2 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
    v3 = LinearViscous()
    v4 = LinearViscous(η=1e22Pa*s)
    e1 = ConstantElasticity()           # elasticity
    e2 = SetConstantElasticity(; G=5e10, Kb=1e11)
    #pl1= DruckerPrager(C=1e6)                # plasticity
    pl1= DruckerPrager(C=1e6/cosd(30))        # plasticity which ends up with the same yield stress as pl3
    pl2= DruckerPrager(C=1e6, ϕ=0, Ψ=10)      # plasticity
    pl3= DruckerPrager(C=1e6, ϕ=0)            # plasticity
    
    # Parallel elements
    p1 = Parallel(v3,v4)                # linear elements
    p2 = Parallel(v1,v2)                # includes nonlinear viscous elements
    p3 = Parallel(v1,v2,v3)             # includes nonlinear viscous elements
    p4 = Parallel(pl1, LinearViscous(η=1e20Pa*s)) # viscoplastic regularisation

    # CompositeRheologies
    c1 = CompositeRheology(v1,v2)
    c2 = CompositeRheology(v3,v4)       # two linear rheologies
    c3 = CompositeRheology(v1,v2, e1)   # with elasticity
    c4 = CompositeRheology(v1,v3, p1)   # with linear || element
    c5 = CompositeRheology(v1,v4, p2)   # with nonlinear || element
    c6 = CompositeRheology(v1,v4,p1,p2) # with 2 || elements
    c7 = CompositeRheology(v4,e1)       # viscoelastic with linear viscosity
    c8 = CompositeRheology(v4,e1,pl1)   # with plastic element
    c9 = CompositeRheology(v4,e1,p4)    # with visco-plastic parallel element
    c10= CompositeRheology(e1,pl3)      # elastoplastic
    c11= CompositeRheology(e1,Parallel(pl3,LinearViscous(η=1e19Pa*s)))      # elasto-viscoplastic
    
    c12= CompositeRheology(e2,v3)       # viscoelasticity with volumetric elasticity
    c13= CompositeRheology(e2,pl2)      # volumetric elastoplastic
    
    c14= CompositeRheology(SetConstantElasticity(G=1e10, Kb=2e11), LinearViscous(η=1e20), DruckerPrager(C=3e5, Ψ=10))   # case A
    c15= CompositeRheology(SetConstantElasticity(G=1e10, Kb=2e11), LinearViscous(η=1e20), Parallel(DruckerPrager(C=3e5, Ψ=10),LinearViscous(η=1e19Pa*s)))   # case A
    c16= CompositeRheology(SetConstantElasticity(G=1e10, Kb=2e11), LinearViscous(η=1e20), DruckerPrager_regularised(C=3e5, Ψ=10, η_vp=1e19))   # case A

    p4 = Parallel(c3,v3)                # Parallel element with composite one as well    

    # Check that we can construct complicated rheological elements
    c = CompositeRheology( (v1, v2, v3, e1, Parallel(p1, v1, v2),v2, Parallel(p1, v1), v2,v3 ))   
    @test isa(c.elements[1], AbstractCreepLaw)

    c = CompositeRheology( (v1, v2, v3, e1, Parallel(p1, e1, Parallel( CompositeRheology(v1, v2), v3) ), v2,v3) )   
    @test isa(c.elements[3], AbstractCreepLaw)
    
    c = Parallel(CompositeRheology(v2,v3,e1, Parallel(v2,v3),v2, Parallel(v2,CompositeRheology(v3, v2))),v3,CompositeRheology(e1,p1))
    @test isa(c.elements[2], AbstractCreepLaw)
    
    args = (T=900.0, d=100e-6, τII_old=1e6, dt=1e8)
    εII, τII = 2e-15, 2e6

    # test volumetric parts
    test_vec = [v1 v2 v3 v4 e1 e2 pl1 pl2]
    sol      = [false false false false false true false true]
    for i = 1:length(sol)
        @test isvolumetric(test_vec[i]) == sol[i]
    end

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
        τ0, = compute_τII(v, εII, args, verbose=false);
        τ1, = compute_τII(v, εII + Δε, args, verbose=false);
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
        τ,    = compute_τII(v, εII, args)   
        τ_AD, = compute_τII_AD(v, εII, args)   
        ε     = compute_εII(v, τ,   args)   
        ε_AD  = compute_εII_AD(v, τ,   args)   
        
        @test ε ≈ εII ≈ ε_AD 
        @test τ ≈ τ_AD

        if isa(v,Parallel)
            # check; for || elements, ε is constant and τ is the sum
            τ_check = 0;
            for i=1:length(v.elements)
                τ_check += compute_τII(v.elements[i], εII, args)[1]   
            end
            @test τ_check ≈ τ
        end
        
        # Check computational routines if stress is given
        ε,    = compute_εII(v, τII, args)   
        ε_AD, = compute_εII_AD(v, τII, args)   
        τ,    = compute_τII(v, ε,   args)   
        τ_AD, = compute_τII_AD(v, ε,   args)   
        
        @test τ ≈ τII ≈ τ_AD
        @test ε ≈ ε_AD
        
        if isa(v,CompositeRheology)
            # check; for serial elements, ε is constant and τ is the sum
            ε_check = 0;
            for i=1:length(v.elements)
                ε_check += compute_εII(v.elements[i], τII, args)[1]   
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
    εII =  1e-15
#    for v in [c4 c5 c6]    
    for v in [c4 c5]    


        τ_AD, = compute_τII_AD(v, εII, args)     # using AD
       
        # check result. For a parallel element we should satisfy the following equations:
        #   τ_parallel == τ_AD 
        #   sum(ε) = εII
        ε_vec = [compute_εII(v.elements[i], τ_AD, args) for i=1:length(v.elements)]
        @test sum(ε_vec) ≈ εII

        ε_parallel,τ_parallel = 0.0, 0.0
        for i=1:length(v.elements)
            # Check that the stress of each || element is the same as the total one
            if isa(v.elements[i], Parallel)
                ε_parallel  = ε_vec[i]      # parallel strainrate
                τ_parallel, = compute_τII_AD(v.elements[i], ε_parallel, args) 
                @test τ_parallel  ≈ τ_AD
            end
        end

        @test εII ≈ compute_εII(v,τ_AD,args)    # sum of ε

        # Check with analytical jacobian
        τ,   = compute_τII(v, εII, args)        # using analytical jacobians (expanded for || elements)
        @test τ_AD ≈ τ

    end

    # Test Parallel rheology with a CompositeRheology branch
    εII =  3e-15
    for v in [p4]
        τ_AD, = compute_τII_AD(v, εII, args)     
        τ,    = compute_τII(v, εII, args)     

        τ_par = 0.
        for i=1:length(v.elements)
            τ_par  += compute_τII(v.elements[i], εII, args)[1]  
        end

        @test τ_AD ≈ τ ≈ τ_par
        
    end

    # Composite cases with (non-parallel) plasticity 
    εII =  3e-15
    args = merge(args, (τII_old=7e5,P=0.0, dt=8e8))
    for v in [c10 c8]
       # τ_AD, = compute_τII_AD(v, εII, args)     
        τ,    = compute_τII(v, εII, args, verbose=false)   
        
        args_old = merge(args, (τII=args.τII_old,))  
        Fold = compute_yieldfunction(c8.elements[3],args_old)

        args = merge(args, (τII=τ,))  
        F    = compute_yieldfunction(c8.elements[3],args)

        # Check that F==0    
        F = 0.
        for i=1:length(v.elements)
            if isa(v.elements,AbstractPlasticity)
                F  += compute_τII(v.elements[i], args)[1]  
            end
        end
        @test F ≈ 0.0

    end

    # cases with plastic elements that have parallel elements as well  
    εII = 1e-15  
    args = (T = 900.0, d = 0.0001, τII_old = 700000.0, dt = 8.0e9, P = 0.0)
    for v in [c8, c9, c10, c11]     

        v_pl = v[length(v.elements)];   # assuming it is the last element
        if isa(v_pl,Parallel)
            v_pl = v_pl[1]
            τ,λ,τ_plastic = compute_τII(v, εII, args, verbose=false)   
        else
            τ,λ           = compute_τII(v, εII, args, verbose=false)  
            τ_plastic = τ
        end

        # check that the sum of strainrates is correct 
        ε_check = 0.0
        for i=1:length(v.elements)
            if !isplastic(v[i])
                ε_check += compute_εII(v[i],τ,args)
            else
                # plastic element
                ε_check += ∂Q∂τII(v_pl, τ_plastic)*λ
            end
        end
        @test ε_check ≈ εII

        # in case we have a || plastic element, check that the sum(τ_parallel)==τ
        if length(v.elements)>2
            if  isa(v[3],Parallel)
                τ_check = 0.0
                v_par = v[3]
                ε_pl = ∂Q∂τII(v_pl, τ_plastic)*λ        # plastic strainrate of || element
                for i=1:length(v_par.elements)
                    if !isplastic(v_par[i])
        
                        τ_check += compute_τII(v_par[i],ε_pl,args)[1]
                    else
                        # plastic element
                        τ_check += τ_plastic
                    end
                end
                @test τ_check ≈ τ       # sum of stress of || element should be the same as 

            end
        end

        args_new = merge(args, (τII=τ_plastic,))
        
        F = compute_yieldfunction(v_pl,args_new)
        @test abs(F)<1e-10

        # check that if we start @ the yield stress, we stay @ the same stress
        
        # perform a 0D time-dependent test and make sure that we remain on the yield stress
        t_vec, τ_vec    =   time_τII_0D(v, εII, args; t=(0.,args.dt*10), nt=50, verbose=false)

        @test τ_vec[end-3] ≈ τ_vec[end]

    end
      
       

    # 0D rheology tests
    η,G  =  10, 1;
    t_M  =  η/G
    εII  =  1.;
    args = (;)
    pl2  =  DruckerPrager(C=η, ϕ=0)                # plasticity

    c_lin =  CompositeRheology(LinearViscous(η=η*Pa*s),ConstantElasticity(G=G*Pa)) # linear VE
    t_vec, τ_vec    =   time_τII_0D(c_lin, εII, args; t=(0.,t_M*4), nt=100, verbose=false)
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
    
    # Compare regularised DP vs. doing the same with a Parallel element
    εII = 1e-15;  
    args = (T = 900.0, d = 0.0001, τII_old = 700000.0, dt = 8.0e9, P = 0.0);
    t_max = 8e11;
    
    # Test elasto-plastic combinations with regularisations 
    # (either done explicitly using a Parallel element or doing the same with DruckerPrager_regularised)
    c_ep   = CompositeRheology(ConstantElasticity(),DruckerPrager() ); 
    c_ep_reg = CompositeRheology(ConstantElasticity(),DruckerPrager_regularised(η_vp=1e20) ); 
    c_e_vp = CompositeRheology(ConstantElasticity(),Parallel(DruckerPrager(),LinearViscous(η=1e20)) ); 
    
    _, τ_vec  =   time_τII_0D(c_ep,     εII, args; t=(0.,t_max), verbose=false)
    _, τ_vec1 =   time_τII_0D(c_e_vp,   εII, args; t=(0.,t_max), verbose=false)
    _, τ_vec2 =   time_τII_0D(c_ep_reg, εII, args; t=(0.,t_max), verbose=false)
    @test τ_vec1[end] ≈  τ_vec2[end]
    @test τ_vec1[end] > τ_vec[end] 

    # Test this with various flavors of viscoelastoplasticity if adding a parallel element
    c_vep = CompositeRheology(LinearViscous(η=1e22),ConstantElasticity(),DruckerPrager());  
    c_vep_reg = CompositeRheology(LinearViscous(η=1e22),ConstantElasticity(),DruckerPrager_regularised(η_vp=1e20) ); 
    c_ve_vp = CompositeRheology(LinearViscous(η=1e22),ConstantElasticity(),Parallel(DruckerPrager(),LinearViscous(η=1e20)) ); 
    
    _, τ_vec  =   time_τII_0D(c_vep, εII, args; t=(0.,t_max), nt=10, verbose=false)
    _, τ_vec1 =   time_τII_0D(c_ve_vp, εII, args; t=(0.,t_max), nt=10, verbose=false)
    _, τ_vec2 =   time_τII_0D(c_vep_reg, εII, args; t=(0.,t_max), nt=10, verbose=false)
    @test τ_vec1[end] ≈  τ_vec2[end]
    @test τ_vec1[end] > τ_vec[end] 

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


    
    # Test cases with/without volumetric elements
    # Note that the only 'special' case implemented one where we have 
    # volumetric plasticity, which requires iterations. 
    # In all other cases we assume the coupling  
#    εII = 1e-15  
#    εvol = -1e-18;

#    εxx,εzz  = 6.8e-15, -7e-15
    εxx,εzz  = 7e-15, -6.8e-15

    εII  = sqrt(0.5*(εxx^2 + εzz^2))
    εvol = εxx + εzz

    args = (T = 900.0, d = 0.0001, τII_old = 700000.0, dt = 8.0e9, P = 0.0, P_old = 1e6)
    for v in [c8 c3 c12 c13 c14]     

        p, τII = compute_p_τII(v, εII, εvol, args, verbose=false)
        if !isvolumetric(v)
            @test p == args.P_old
        end

        # 0D rheology functions p
        t_max = args.dt*2;
        t_vec, P_vec, τ_vec    =   time_p_τII_0D(v, εII, εvol, args; t=(0.,t_max), nt=20, verbose=false)
        if !isvolumetric(v)
            @test sum(P_vec) == 0.0
        
        elseif isvolumetric(v) && !isvolumetricplastic(v)
            # 'uncoupled' volumetric deformation
            Kb_computed = (P_vec[end] - P_vec[1])/(-t_max*εvol)
            if isa(v[1],AbstractElasticity)
                @test Kb_computed  ≈  NumValue(v[1].Kb)
            end

        end

    end

    # case with dilatant plasticity
    εxx,εzz  = 7e-15, -6.8e-15
    εII  = sqrt(0.5*(εxx^2 + εzz^2))
    εvol = εxx + εzz
    args = (T = 900.0, d = 0.0001, τII_old = 700000.0, dt = 8.0e9, P = 0.0, P_old = 1e6)
    
    t_max = args.dt*2;
    _, P_vec, τ_vec    =   time_p_τII_0D(c14, εII, εvol, args; t=(0.,t_max), nt=20, verbose=false)
    #@test  sum(P_vec) ≈ 1.1803339314328427e6
    #@test  sum(τ_vec) ≈ 4.742438875602647e6
    
    # dilatant plasticity in || with viscous element
    _, P_vec1, τ_vec1  =   time_p_τII_0D(c15, εII, εvol, args; t=(0.,t_max), nt=20, verbose=false)
    _, P_vec2, τ_vec2  =   time_p_τII_0D(c16, εII, εvol, args; t=(0.,t_max), nt=20, verbose=false)

    @test τ_vec1[end] ≈  τ_vec2[end]
    @test τ_vec1[end] > τ_vec[end] 

end
