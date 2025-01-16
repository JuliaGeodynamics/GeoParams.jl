using Test
using GeoParams

@testset "ChemicalDiffusion" begin

    # test the diffusion data structure
    x1 = DiffusionData()
    @test isbits(x1)

    # check auto unit conversion
    Hf_Rt_perp = Rutile.Rt_Hf_Cherniak2007_⊥c;
    Hf_Rt_perp = SetChemicalDiffusion(Hf_Rt_perp; D0=10km^2/s)
    @test Hf_Rt_perp.D0.val == 1.0e7

    # test the diffusion parameter calculation
    D = ustrip(compute_D(x1))
    @test D == 0.0

    # test the diffusion parameter calculation with arrays
    Hf_Rt_para = Rutile.Rt_Hf_Cherniak2007_Ξc;
    Hf_Rt_para = SetChemicalDiffusion(Hf_Rt_para)

    # with unit
    T = ones(10) * 1273.15K
    D = zeros(10)m^2/s
    compute_D!(D, Hf_Rt_para;T=T)
    @test ustrip(D[1]) ≈ 1.06039e-21 atol = 1e-24

    # without unit
    D = zeros(10)
    T = ones(10) * 1273.15
    compute_D!(D, Hf_Rt_para;T=T, P=zeros(size(T)))
    @test D[1] ≈ 1.06039e-21 atol = 1e-24

    # test experimental data with literature values

    # Benchmark Rutile Hf data from Cherniak 2007 (HD 15/01/25)
    Hf_Rt_para = Rutile.Rt_Hf_Cherniak2007_Ξc;
    Hf_Rt_para = SetChemicalDiffusion(Hf_Rt_para)
    D =  ustrip(compute_D(Hf_Rt_para, T=1273.15K))
    @test D ≈ 1.06039e-21 atol = 1e-24

    Hf_Rt_perp = Rutile.Rt_Hf_Cherniak2007_⊥c;
    Hf_Rt_perp = SetChemicalDiffusion(Hf_Rt_perp)
    D =  ustrip(compute_D(Hf_Rt_perp, T=1273.15K))
    @test D ≈ 1.21560e-21 atol = 1e-24

    # Hf_Rt_para = Rutile.Rt_Hf_Cherniak2007_Ξc;
    # Hf_Rt_para = SetChemicalDiffusion(Hf_Rt_para)

    # CharUnits_GEO   =   GEO_units(length=10u"cm")
    # SetMaterialParams(Name="Viscous Matrix",
    #                   ChemDiffusion   = Hf_Rt_para,
    #                   CharUnits_GEO)
end



