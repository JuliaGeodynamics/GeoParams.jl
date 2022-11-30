using StaticArrays, GeoParams, ForwardDiff

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
c = CompositeRheology( (v1, v2, v3, e1, Parallel(p1, e1, Parallel( CompositeRheology(v1, v2), v3) ), v2,v3) )   
c = Parallel(CompositeRheology(v2,v3,e1, Parallel(v2,v3),v2, Parallel(v2,CompositeRheology(v3, v2))),v3,CompositeRheology(e1,p1))
    
args = (T=900.0, d=100e-6, τII_old=1e6, dt=1e8)
εII, τII = 2e-15, 2e6
# compute_τII(v, εII, args, verbose=false);

c = CompositeRheology(e1,v4)       # two linear rheologies
ε
εij, τij, τij_old  = ntuple(x -> rand(3), 3) 
compute_ε(e1, τij, τij_old, args)
compute_τ(e1, εij, τij_old, args)
compute_dτdε(e1, εij, τij_old, args)
compute_dεdτ(e1, εij, τij_old, args)

compute_ε(v4, τij, args)
compute_τ(v4, εij, args)
compute_dτdε(v4, εij, args)
compute_dεdτ(v4, εij, args)

εij, τij, τij_old  = ntuple(x ->Tuple(rand(3)), 3) 
compute_ε(e1, τij, τij_old, args)
compute_τ(e1, εij, τij_old, args)
compute_dτdε(e1, εij, τij_old, args)
compute_dεdτ(e1, εij, τij_old, args)

compute_ε(v4, τij, args)
compute_τ(v4, εij, args)
compute_dτdε(v4, εij, args)
compute_dεdτ(v4, εij, args)

εij, τij, τij_old  = ntuple(x -> SVector{3,Float64}(rand(3)), 3) 
compute_ε(e1, τij, τij_old, args)
compute_τ(e1, εij, τij_old, args)
compute_dτdε(e1, εij, τij_old, args)
compute_dεdτ(e1, εij, τij_old, args)

compute_ε(v4, τij, args)
compute_τ(v4, εij, args)
compute_dτdε(v4, εij, args)
compute_dεdτ(v4, εij, args)
