using Test
using GeoParams

# measure inside a function: on Julia 1.10, @allocated at test scope counts boxing of the result
allocs(f::F, args::Vararg{Any, N}) where {F, N} = (f(args...); @allocated f(args...))

@testset "WaterFugacity" begin

    # Fugacity of H2O [MPa] from the online calculator of A. Withers
    # (https://www.esci.umn.edu/people/researchers/withe012/fugacity.htm), which implements
    # Pitzer & Sterner (1994). Rows: P = 100:100:2800 MPa, columns: T = 100:100:1500 °C.
    fH_Withers = [
        0.1805 2.4001 10.5506 26.2335 46.0866 64.0398 77.0565 85.7693 91.5531 95.4524 98.1368 100.0236 101 102 103
        0.3197 3.8892 16.4169 40.3115 72.1575 105.717 135.5804 159.2991 176.9813 189.7775 198.936 205.479 210 214 216
        0.5573 6.1846 24.905 59.5835 105.5426 155.1193 201.62 241.0616 272.2636 295.8816 313.296 325.937 335 341 346
        0.9588 9.695 37.1336 86.1191 149.721 218.05 282.9415 339.3973 385.4246 421.2079 448.099 467.81 482 492 499
        1.6311 15.0232 54.6422 122.5294 208.3012 299.0737 384.972 460.181 522.2919 571.2836 608.53 635.988 656 669 678
        2.7484 23.0544 79.5544 172.2363 285.7386 403.38 513.4012 609.3546 688.7611 751.7224 799.827 835.339 861 878 889
        4.5915 35.083 114.7944 239.7427 387.6674 537.2665 674.8434 793.601 891.3743 968.7492 1027.79 1071.22 1102 1122 1135
        7.6122 52.9943 164.3747 330.9664 521.2406 708.4863 877.2081 1020.711 1137.6272 1229.441 1299.01 1349.72 1385 1408 1421
        12.532 79.5228 233.7834 453.6668 695.5345 926.6135 1130.032 1299.884 1436.2226 1541.924 1620.1 1677.7 1716 1740 1752
        20.502 118.62 330.5029 617.9914 922.0459 1203.467 1444.839 1642.027 1797.3189 1915.566 2002.32 2063.01 2102 2125 2135
        33.345 175.977 464.7042 837.1759 1215.307 1553.62 1835.541 2060.086 2232.7949 2361.19 2452.82 2514.39 2552 2571 2575
        53.936 259.761 650.1716 1128.442 1593.649 1995.007 2318.929 2569.418 2756.5453 2891.307 2983.73 3042.41 3074 3086 3081
        86.79 381.66 905.5 1512.14 2080 2549.7 2915 3188 3384 3520 3608 3658 3679 3678 3660
        139.02 558.33 1255.9 2023.2 2703 3245 3649 3938 4137 4265 4340 4375 4378 4359 4322
        221.66 813.488 1735.08 2693.07 3501 4113 4549 4844 5034 5145 5197 5207 5184 5139 5077
        351.92 1180.75 2388.34 3571.93 4517 5195 5650 5937 6103 6181 6198 6171 6112 6032 5937
        556.48 1707.67 3276.3 4721.84 5809 6541 6995 7252 7373 7401 7365 7286 7179 7052 6913
        876.5 2461 4480 6222 7447 8209 8634 8831 8880 8831 8722 8574 8401 8215 8020
        1376 3536 6107 8176 9520 10274 10627 10724 10663 10509 10299 10058 9802 9538 9274
        2152 5065 8302 10713 12138 12826 13045 12987 12771 12471 12127 11767 11402 11042 10692
        3354 7233 11254 14001 15436 15970 15973 15688 15258 14762 14244 13729 13229 12749 12292
        5212 10301 15216 18252 19582 19838 19514 18908 18186 17432 16691 15981 15311 14683 14097
        8073 14631 20521 23737 24786 24589 23787 22740 21629 20542 19514 18560 17680 16872 16129
        12468 20728 27612 30802 31304 30412 28936 27292 25672 24156 22768 21511 20374 19346 18413
        19201 29291 37065 39883 39456 37540 35130 32691 30412 28351 26514 24882 23432 22138 20979
        29488 41295 49648 51535 49632 46250 42571 39088 35961 33215 30820 28730 26899 25286 23857
        45166 58085 66362 66458 62315 56877 51497 46654 42450 38847 35764 33115 30826 28832 27084
        69002 81519 88523 85540 78096 69824 62189 55592 50028 45360 41433 38108 35268 32821 30694
    ]
    P_MPa = 100:100:2800
    T_C = 100:100:1500

    fH = [compute_water_fugacity(P * 1.0e6, T + 273.15) / 1.0e6 for P in P_MPa, T in T_C]
    max_dev = maximum(abs.(fH ./ fH_Withers .- 1))
    @test max_dev < 5.0e-3   # the table is rounded to integers at high T; max. deviation is ~0.4%

    # Float32 and Unitful versions
    f64 = compute_water_fugacity(1.0e9, 1073.15)
    f32 = compute_water_fugacity(1.0f9, 1073.15f0)
    @test f32 isa Float32
    @test f32 ≈ f64 rtol = 1.0e-4
    @test compute_water_fugacity(1.0e9, 1073.15) ≈ fH_Withers[10, 8] * 1.0e6 rtol = 1.0e-4
    fq = compute_water_fugacity(1GPa, 800C)
    @test unit(fq) == Pa
    @test ustrip(fq) ≈ f64
    @test ustrip(MPa, compute_water_fugacity(1000MPa, 1073.15K)) ≈ f64 / 1.0e6

    # no allocations
    compute_water_fugacity(1.0e9, 1073.15)
    @test allocs(compute_water_fugacity, 1.0e9, 1073.15) == 0

    # approaches an ideal gas at low pressure and high temperature
    @test compute_water_fugacity(1.0e5, 2273.15) ≈ 1.0e5 rtol = 1.0e-3
    # stable (vapour) root below the saturation pressure, liquid above it (Psat(300 °C) ≈ 8.6 MPa)
    @test compute_water_fugacity(5.0e6, 573.15) / 5.0e6 > 0.85
    @test compute_water_fugacity(2.0e7, 573.15) < 1.0e7

    # Shinevar et al. (2015) crustal approximation, within a factor 1.5 of PS94
    # for 0.2-1 GPa and 600-1000 K
    for P in (0.2e9, 0.4e9, 0.6e9, 0.8e9, 1.0e9), T in 600.0:50.0:1000.0
        ratio = compute_water_fugacity_Shinevar(P, T) / compute_water_fugacity(P, T)
        @test 1 / 1.5 < ratio < 1.5
    end
    @test compute_water_fugacity_Shinevar(1.0e9, 800.0) ≈ 5.521e9 * exp(-(31.28e3 - 1.0e9 * 2.009e-5) / (8.3145 * 800.0))
    @test compute_water_fugacity_Shinevar(1.0f9, 800.0f0) isa Float32
    @test ustrip(Pa, compute_water_fugacity_Shinevar(1GPa, 800K)) ≈ compute_water_fugacity_Shinevar(1.0e9, 800.0)

    # Water content of olivine (Zhao et al. 2004), C_OH in H/10⁶ Si with fH2O in MPa
    P, T = 2.0e9, 1573.15
    fH2O = compute_water_fugacity(P, T)
    C_OH = compute_water_content_olivine(fH2O, P, T)
    @test C_OH ≈ 26 * exp(-(40.0e3 + P * 10.0e-6) / (8.3145 * T)) * fH2O / 1.0e6
    @test 1000 < C_OH < 10000
    @test compute_water_content_olivine(fH2O * Pa, 2GPa, 1300C) ≈ C_OH
    @test compute_water_content_olivine(Float32(fH2O), 2.0f9, 1573.15f0) isa Float32
end
