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

    
    εII_array       =   ones(10)*1e-5
    τII_array       =   similar(εII_array)
    
    args_array = (;T=T_array )
    compute_τII!(τII_array, v,εII_array, args) 



end
