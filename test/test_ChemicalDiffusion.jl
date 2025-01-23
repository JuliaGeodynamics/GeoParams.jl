using Test
using GeoParams

@testset "ChemicalDiffusion" begin

    # test the diffusion data structure
    x1 = DiffusionData()
    @test isbits(x1)

    # check auto unit conversion
    Hf_Rt_perp = Rutile.Rt_Hf_Cherniak2007_perp_c
    Hf_Rt_perp = SetChemicalDiffusion(Hf_Rt_perp; D0 = 10km^2 / s)
    @test Hf_Rt_perp.D0.val == 1.0e7

    # test the diffusion parameter calculation
    D = ustrip(compute_D(x1))
    @test D == 0.0

    # test the diffusion parameter calculation with arrays
    Hf_Rt_para = Rutile.Rt_Hf_Cherniak2007_para_c
    Hf_Rt_para = SetChemicalDiffusion(Hf_Rt_para)

    # with unit
    T = ones(10) * 1273.15K
    D = zeros(10)m^2 / s
    compute_D!(D, Hf_Rt_para; T = T)
    @test ustrip(D[1]) ≈ 1.06039e-21 atol = 1.0e-24

    # without unit
    D = zeros(10)
    T = ones(10) * 1273.15
    compute_D!(D, Hf_Rt_para; T = T, P = zeros(size(T)))
    @test D[1] ≈ 1.06039e-21 atol = 1.0e-24

    # test SetMaterialParams
    phase = SetMaterialParams(
        Name = "Chemical Diffusion",
        ChemDiffusion = Hf_Rt_para
    )

    @test phase.ChemDiffusion[1].D0.val ≈ 9.1e-15

    # test nondimensionalisation
    CharUnits_GEO = GEO_units(length = 10cm)
    phase = SetMaterialParams(
        Name = "Chemical Diffusion",
        ChemDiffusion = Hf_Rt_para,
        CharDim = CharUnits_GEO
    )

    @test phase.ChemDiffusion[1].D0.val ≈ 9.099999999999998

    # test experimental data with literature values

    # -------------------------- Rutile --------------------------

    # Benchmark Hf data from Cherniak 2007 (HD 15/01/25)
    Hf_Rt_para = Rutile.Rt_Hf_Cherniak2007_para_c
    Hf_Rt_para = SetChemicalDiffusion(Hf_Rt_para)
    D = ustrip(compute_D(Hf_Rt_para, T = 1273.15K))
    @test D ≈ 1.06039e-21 atol = 1.0e-24

    Hf_Rt_perp = Rutile.Rt_Hf_Cherniak2007_perp_c
    Hf_Rt_perp = SetChemicalDiffusion(Hf_Rt_perp)
    D = ustrip(compute_D(Hf_Rt_perp, T = 1273.15K))
    @test D ≈ 1.2156e-21 atol = 1.0e-24

    # Benchmark Zr data from Cherniak 2007 (HD 20/01/25)
    Zr_Rt = Rutile.Rt_Zr_Cherniak2007_para_c
    Zr_Rt = SetChemicalDiffusion(Zr_Rt)
    D = ustrip(compute_D(Zr_Rt, T = 1273.15K))
    @test D ≈ 1.0390187e-21 atol = 1.0e-24

    # -------------------------- Garnet --------------------------

    # Benchmark Fe data from Chakraborty 1992 (HD 18/01/25)
    Fe_Grt = Garnet.Grt_Fe_Chakraborty1992
    Fe_Grt = SetChemicalDiffusion(Fe_Grt)
    D = ustrip(uconvert(cm^2 / s, compute_D(Fe_Grt, T = 1373.15K, P = 1GPa)))
    @test D ≈ 1.308812e-14 atol = 1.0e-18

    # Benchmark Mg data from Chakraborty 1992 (HD 18/01/25)
    Mg_Grt = Garnet.Grt_Mg_Chakraborty1992
    Mg_Grt = SetChemicalDiffusion(Mg_Grt)
    D = ustrip(uconvert(cm^2 / s, compute_D(Mg_Grt, T = 1373.15K, P = 1GPa)))
    @test D ≈ 1.041487e-14 atol = 1.0e-18

    # Benchmark Mn data from Chakraborty 1992 (HD 18/01/25)
    Mn_Grt = Garnet.Grt_Mn_Chakraborty1992
    Mn_Grt = SetChemicalDiffusion(Mn_Grt)
    D = ustrip(uconvert(cm^2 / s, compute_D(Mn_Grt, T = 1373.15K, P = 1GPa)))
    @test D ≈ 6.909072e-14 atol = 1.0e-18

    # Benchmark REE data from Bloch 2020 (HD 23/01/25)
    REE_Grt_fast = Garnet.Grt_REE_Bloch2020_fast
    REE_Grt_fast = SetChemicalDiffusion(REE_Grt_fast)
    D = ustrip(compute_D(REE_Grt_fast, T = 1323.15K, P = 1.0GPa))
end

# REE_Grt_fast = Garnet.Grt_REE_Bloch2020_fast
# REE_Grt_fast = SetChemicalDiffusion(REE_Grt_fast)


# D = ustrip(compute_D(REE_Grt_fast, T = 1323.15K, P = 1.0GPa))


