using Test
using GeoParams
    
@testset "CompositeRheologies" begin

    # Diffusion & dislocation creep
    pp   = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
    pp1  = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
    v    = (pp,pp1)
    args = (T=900.0, d=100e-6)
    εII  = 1e-15
    
    τII  = compute_τII(v,εII, args) 
    @test τII ≈ 8.892678156850036e8
    η = τII/(2εII)

    # Same with arrays:    
    εII_array       =   ones(10)*1e-5
    τII_array       =   similar(εII_array)
    compute_τII!(τII_array, v,εII_array, args) 
    @test τII_array[1] ≈ 1.918028581543394e12

    # Add elasticity in the mix
    el      = ConstantElasticity();
    v         = (pp,pp1, el)
    dt_maxwell= η/NumValue(el.G);  
    args    = (T=900.0, d=100e-6,  τII_old = 0.0, dt=dt_maxwell/20)
    εII     = 1e-15
    
    τII_vec = [0.0];
    t       = [0.0];
    for i=1:100
        τII  = compute_τII(v,εII, args)
        args = merge(args, (;τII_old=τII))
        τII_vec = [τII_vec; τII]
        t       = [t; i*args.dt]
    end
    SecYear = 3600*24*365.25
    @test sum(τII_vec) ≈ 8.016289174454965e9


end
