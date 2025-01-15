using Test
using GeoParams

@testset "ChemicalDiffusion" begin

    # test the diffusion data structure
    x1 = GeoParams.DiffusionData()
    @test isbits(x1)
    @test Value(x1.D0) == 0.0m^2/s

    # test Rutile Hf data
    Hf_Rt_perp = GeoParams.Rutile.Rt_Hf_Cherniak2007_⊥c()
    # calculate D and test here

    # test Rutile Hf data
    Hf_Rt_para = GeoParams.Rutile.Rt_Hf_Cherniak2007_Ξc()
    # calculate D and test here


end