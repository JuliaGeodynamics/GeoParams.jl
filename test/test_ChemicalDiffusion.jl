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

    # Benchmark Hf data from Cherniak (2007) (HD 15/01/25)
    Hf_Rt_para = Rutile.Rt_Hf_Cherniak2007_para_c
    Hf_Rt_para = SetChemicalDiffusion(Hf_Rt_para)
    D = ustrip(compute_D(Hf_Rt_para, T = 1273.15K))
    @test D ≈ 1.06039e-21 atol = 1.0e-24

    Hf_Rt_perp = Rutile.Rt_Hf_Cherniak2007_perp_c
    Hf_Rt_perp = SetChemicalDiffusion(Hf_Rt_perp)
    D = ustrip(compute_D(Hf_Rt_perp, T = 1273.15K))
    @test D ≈ 1.2156e-21 atol = 1.0e-24

    # Benchmark Zr data from Cherniak (2007) (HD 20/01/25)
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
    REE_Grt_slow = Garnet.Grt_REE_Bloch2020_slow
    REE_Grt_slow = SetChemicalDiffusion(REE_Grt_slow)
    D = ustrip(compute_D(REE_Grt_slow, T = 1323.15K))

    D_paper = exp(-10.24 - (221057)/ (2.303*8.31446261815324 * 1323.15))
    @test D ≈ D_paper

    # Benchmark REE data from Bloch 2020 (HD 23/01/25)
    REE_Grt_fast = Garnet.Grt_REE_Bloch2020_fast
    REE_Grt_fast = SetChemicalDiffusion(REE_Grt_fast)
    D = ustrip(compute_D(REE_Grt_fast, T = 1323.15K, P = 1.0GPa))

    D_paper = exp(-9.28 - (265200 + 10800*1) / (2.303*8.31446261815324 * 1323.15))
    @test D ≈ D_paper

    # -------------------------- Rhyolitic Melt --------------------------

    # Benchmark Sc data from Holycross and Watson 2018 (HD 24/01/25)
    Sc_Melt_highH2O = Melt.Melt_Sc_Holycross2018_rhyolitic_highH2O
    Sc_Melt_highH2O = SetChemicalDiffusion(Sc_Melt_highH2O)

    D = ustrip(compute_D(Sc_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 8.492199e-14 atol = 1.0e-18

    Sc_Melt_mediumH2O = Melt.Melt_Sc_Holycross2018_rhyolitic_mediumH2O
    Sc_Melt_mediumH2O = SetChemicalDiffusion(Sc_Melt_mediumH2O)

    D = ustrip(compute_D(Sc_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈ 1.556627e-13 atol = 1.0e-18

    # Benchmark V data from Holycross and Watson 2018 (HD 24/01/25)
    V_Melt_highH2O = Melt.Melt_V_Holycross2018_rhyolitic_highH2O
    V_Melt_highH2O = SetChemicalDiffusion(V_Melt_highH2O)

    D = ustrip(compute_D(V_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 9.746426609e-14 atol = 1.0e-18

    V_Melt_mediumH2O = Melt.Melt_V_Holycross2018_rhyolitic_mediumH2O
    V_Melt_mediumH2O = SetChemicalDiffusion(V_Melt_mediumH2O)

    D = ustrip(compute_D(V_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈ 3.23892920e-13 atol = 1.0e-18

    # Benchmark Y data from Holycross and Watson 2018 (HD 24/01/25)
    Y_Melt_highH2O = Melt.Melt_Y_Holycross2018_rhyolitic_highH2O
    Y_Melt_highH2O = SetChemicalDiffusion(Y_Melt_highH2O)

    D = ustrip(compute_D(Y_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 1.592339304e-13 atol = 1.0e-18

    Y_Melt_mediumH2O = Melt.Melt_Y_Holycross2018_rhyolitic_mediumH2O
    Y_Melt_mediumH2O = SetChemicalDiffusion(Y_Melt_mediumH2O)

    D = ustrip(compute_D(Y_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈ 4.11331120e-13 atol = 1.0e-18

    # Benchmark Zr data from Holycross and Watson 2018 (HD 24/01/25)
    Zr_Melt_highH2O = Melt.Melt_Zr_Holycross2018_rhyolitic_highH2O
    Zr_Melt_highH2O = SetChemicalDiffusion(Zr_Melt_highH2O)

    D = ustrip(compute_D(Zr_Melt_highH2O, T = 1163.15K))
    @test  D ≈ 4.7013667e-14 atol = 1.0e-18

    Zr_Melt_mediumH2O = Melt.Melt_Zr_Holycross2018_rhyolitic_mediumH2O
    Zr_Melt_mediumH2O = SetChemicalDiffusion(Zr_Melt_mediumH2O)

    D = ustrip(compute_D(Zr_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈ 1.1681056405e-13 atol = 1.0e-18

    # Benchmark Hf data from Holycross and Watson 2018 (HD 24/01/25)
    Hf_Melt_mediumH2O = Melt.Melt_Hf_Holycross2018_rhyolitic_mediumH2O
    Hf_Melt_mediumH2O = SetChemicalDiffusion(Hf_Melt_mediumH2O)

    D = ustrip(compute_D(Hf_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈ 9.97822548e-14 atol = 1.0e-18

    # Benchmark Nb data from Holycross and Watson 2018 (HD 24/01/25)
    Nb_Melt_highH2O = Melt.Melt_Nb_Holycross2018_rhyolitic_highH2O
    Nb_Melt_highH2O = SetChemicalDiffusion(Nb_Melt_highH2O)

    D = ustrip(compute_D(Nb_Melt_highH2O, T = 1163.15K))
    @test  D ≈ 1.07505876683e-13 atol = 1.0e-18

    Nb_Melt_mediumH2O = Melt.Melt_Nb_Holycross2018_rhyolitic_mediumH2O
    Nb_Melt_mediumH2O = SetChemicalDiffusion(Nb_Melt_mediumH2O)

    D = ustrip(compute_D(Nb_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈ 1.473878454e-13 atol = 1.0e-18

    # Benchmark La data from Holycross and Watson 2018 (HD 24/01/25)
    La_Melt_highH2O = Melt.Melt_La_Holycross2018_rhyolitic_highH2O
    La_Melt_highH2O = SetChemicalDiffusion(La_Melt_highH2O)

    D = ustrip(compute_D(La_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 2.5040159e-13 atol = 1.0e-18

    La_Melt_mediumH2O = Melt.Melt_La_Holycross2018_rhyolitic_mediumH2O
    La_Melt_mediumH2O = SetChemicalDiffusion(La_Melt_mediumH2O)

    D = ustrip(compute_D(La_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈ 5.8301422e-13 atol = 1.0e-18

    # Benchmark Ce data from Holycross and Watson 2018 (HD 24/01/25)
    Ce_Melt_highH2O = Melt.Melt_Ce_Holycross2018_rhyolitic_highH2O
    Ce_Melt_highH2O = SetChemicalDiffusion(Ce_Melt_highH2O)

    D = ustrip(compute_D(Ce_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 2.4168221e-13 atol = 1.0e-18

    Ce_Melt_mediumH2O = Melt.Melt_Ce_Holycross2018_rhyolitic_mediumH2O
    Ce_Melt_mediumH2O = SetChemicalDiffusion(Ce_Melt_mediumH2O)

    D = ustrip(compute_D(Ce_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈ 5.0854503e-13 atol = 1.0e-18

    # Benchmark Pr data from Holycross and Watson 2018 (HD 24/01/25)
    Pr_Melt_highH2O = Melt.Melt_Pr_Holycross2018_rhyolitic_highH2O
    Pr_Melt_highH2O = SetChemicalDiffusion(Pr_Melt_highH2O)

    D = ustrip(compute_D(Pr_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 2.5893475e-13 atol = 1.0e-18

    Pr_Melt_mediumH2O = Melt.Melt_Pr_Holycross2018_rhyolitic_mediumH2O
    Pr_Melt_mediumH2O = SetChemicalDiffusion(Pr_Melt_mediumH2O)

    D = ustrip(compute_D(Pr_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  5.448527e-13 atol = 1.0e-18

    # Benchmark Nd data from Holycross and Watson 2018 (HD 24/01/25)
    Nd_Melt_highH2O = Melt.Melt_Nd_Holycross2018_rhyolitic_highH2O
    Nd_Melt_highH2O = SetChemicalDiffusion(Nd_Melt_highH2O)

    D = ustrip(compute_D(Nd_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 2.120037e-13 atol = 1.0e-18

    Nd_Melt_mediumH2O = Melt.Melt_Nd_Holycross2018_rhyolitic_mediumH2O
    Nd_Melt_mediumH2O = SetChemicalDiffusion(Nd_Melt_mediumH2O)

    D = ustrip(compute_D(Nd_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  5.6794099e-13 atol = 1.0e-18

    # Benchmark Sm data from Holycross and Watson 2018 (HD 24/01/25)
    Sm_Melt_highH2O = Melt.Melt_Sm_Holycross2018_rhyolitic_highH2O
    Sm_Melt_highH2O = SetChemicalDiffusion(Sm_Melt_highH2O)

    D = ustrip(compute_D(Sm_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 2.016882e-13 atol = 1.0e-18

    Sm_Melt_mediumH2O = Melt.Melt_Sm_Holycross2018_rhyolitic_mediumH2O
    Sm_Melt_mediumH2O = SetChemicalDiffusion(Sm_Melt_mediumH2O)

    D = ustrip(compute_D(Sm_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  5.112780e-13 atol = 1.0e-18

    # Benchmark Eu data from Holycross and Watson 2018 (HD 24/01/25)

    Eu_Melt_mediumH2O = Melt.Melt_Eu_Holycross2018_rhyolitic_mediumH2O
    Eu_Melt_mediumH2O = SetChemicalDiffusion(Eu_Melt_mediumH2O)

    D = ustrip(compute_D(Eu_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  1.423718e-12 atol = 1.0e-18

    # Benchmark Gd data from Holycross and Watson 2018 (HD 24/01/25)
    Gd_Melt_highH2O = Melt.Melt_Gd_Holycross2018_rhyolitic_highH2O
    Gd_Melt_highH2O = SetChemicalDiffusion(Gd_Melt_highH2O)

    D = ustrip(compute_D(Gd_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 1.704840e-13 atol = 1.0e-18

    Gd_Melt_mediumH2O = Melt.Melt_Gd_Holycross2018_rhyolitic_mediumH2O
    Gd_Melt_mediumH2O = SetChemicalDiffusion(Gd_Melt_mediumH2O)

    D = ustrip(compute_D(Gd_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  5.094053e-13 atol = 1.0e-18

    # Benchmark Tb data from Holycross and Watson 2018 (HD 24/01/25)
    Tb_Melt_highH2O = Melt.Melt_Tb_Holycross2018_rhyolitic_highH2O
    Tb_Melt_highH2O = SetChemicalDiffusion(Tb_Melt_highH2O)

    D = ustrip(compute_D(Tb_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 1.700945e-13 atol = 1.0e-18

    Tb_Melt_mediumH2O = Melt.Melt_Tb_Holycross2018_rhyolitic_mediumH2O
    Tb_Melt_mediumH2O = SetChemicalDiffusion(Tb_Melt_mediumH2O)

    D = ustrip(compute_D(Tb_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  4.284803e-13 atol = 1.0e-18

    # Benchmark Dy data from Holycross and Watson 2018 (HD 24/01/25)
    Dy_Melt_highH2O = Melt.Melt_Dy_Holycross2018_rhyolitic_highH2O
    Dy_Melt_highH2O = SetChemicalDiffusion(Tb_Melt_highH2O)

    D = ustrip(compute_D(Dy_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 1.7009448e-13 atol = 1.0e-18

    Dy_Melt_mediumH2O = Melt.Melt_Dy_Holycross2018_rhyolitic_mediumH2O
    Dy_Melt_mediumH2O = SetChemicalDiffusion(Dy_Melt_mediumH2O)

    D = ustrip(compute_D(Dy_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  4.628963e-13 atol = 1.0e-18

    # Benchmark Ho data from Holycross and Watson 2018 (HD 24/01/25)
    Ho_Melt_highH2O = Melt.Melt_Ho_Holycross2018_rhyolitic_highH2O
    Ho_Melt_highH2O = SetChemicalDiffusion(Ho_Melt_highH2O)

    D = ustrip(compute_D(Ho_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 1.463340e-13 atol = 1.0e-18

    Ho_Melt_mediumH2O = Melt.Melt_Ho_Holycross2018_rhyolitic_mediumH2O
    Ho_Melt_mediumH2O = SetChemicalDiffusion(Ho_Melt_mediumH2O)

    D = ustrip(compute_D(Ho_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  4.072542e-13 atol = 1.0e-18

    # Benchmark Er data from Holycross and Watson 2018 (HD 24/01/25)
    Er_Melt_highH2O = Melt.Melt_Er_Holycross2018_rhyolitic_highH2O
    Er_Melt_highH2O = SetChemicalDiffusion(Er_Melt_highH2O)

    D = ustrip(compute_D(Er_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 1.160297e-13 atol = 1.0e-18

    Er_Melt_mediumH2O = Melt.Melt_Er_Holycross2018_rhyolitic_mediumH2O
    Er_Melt_mediumH2O = SetChemicalDiffusion(Er_Melt_mediumH2O)

    D = ustrip(compute_D(Er_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  4.387794e-13 atol = 1.0e-18

    # Benchmark Yb data from Holycross and Watson 2018 (HD 24/01/25)
    Yb_Melt_highH2O = Melt.Melt_Yb_Holycross2018_rhyolitic_highH2O
    Yb_Melt_highH2O = SetChemicalDiffusion(Er_Melt_highH2O)

    D = ustrip(compute_D(Yb_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 1.160297256e-13 atol = 1.0e-18

    Yb_Melt_mediumH2O = Melt.Melt_Yb_Holycross2018_rhyolitic_mediumH2O
    Yb_Melt_mediumH2O = SetChemicalDiffusion(Yb_Melt_mediumH2O)

    D = ustrip(compute_D(Yb_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  4.219549e-13 atol = 1.0e-18

    # Benchmark Lu data from Holycross and Watson 2018 (HD 24/01/25)
    Lu_Melt_highH2O = Melt.Melt_Lu_Holycross2018_rhyolitic_highH2O
    Lu_Melt_highH2O = SetChemicalDiffusion(Lu_Melt_highH2O)

    D = ustrip(compute_D(Lu_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 1.01538143e-13 atol = 1.0e-18

    Lu_Melt_mediumH2O = Melt.Melt_Lu_Holycross2018_rhyolitic_mediumH2O
    Lu_Melt_mediumH2O = SetChemicalDiffusion(Lu_Melt_mediumH2O)

    D = ustrip(compute_D(Lu_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  3.775836e-13 atol = 1.0e-18

    # Benchmark Th data from Holycross and Watson 2018 (HD 24/01/25)
    Th_Melt_highH2O = Melt.Melt_La_Holycross2018_rhyolitic_highH2O
    Th_Melt_highH2O = SetChemicalDiffusion(Th_Melt_highH2O)

    D = ustrip(compute_D(Th_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 2.50401591567e-13 atol = 1.0e-18

    # Benchmark U data from Holycross and Watson 2018 (HD 24/01/25)
    U_Melt_highH2O = Melt.Melt_U_Holycross2018_rhyolitic_highH2O
    U_Melt_highH2O = SetChemicalDiffusion(U_Melt_highH2O)

    D = ustrip(compute_D(U_Melt_highH2O, T = 1123.15K))
    @test  D ≈ 3.0722503094e-14 atol = 1.0e-18

    U_Melt_mediumH2O = Melt.Melt_U_Holycross2018_rhyolitic_mediumH2O
    U_Melt_mediumH2O = SetChemicalDiffusion(U_Melt_mediumH2O)

    D = ustrip(compute_D(U_Melt_mediumH2O, T = 1273.15K))
    @test  D ≈  1.556627767e-13 atol = 1.0e-18

end
