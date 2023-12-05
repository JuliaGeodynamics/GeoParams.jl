using Test
using GeoParams

@testset "EnergyParameters.jl" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=10km)

    # Heat capacity ---------

    # Constant heat capacity
    cp1 = ConstantHeatCapacity()
    info = param_info(cp1)
    @test isbits(cp1)
    @test cp1.cp.val == 1050.0

    cp1_nd = cp1
    cp1_nd = nondimensionalize(cp1_nd, CharUnits_GEO)
    @test cp1_nd.cp.val ≈ 1.3368075000000002e22
    @test compute_heatcapacity(cp1_nd) ≈ 1.3368075000000002e22  # compute
    @test cp1_nd() ≈ 1.3368075000000002e22  # compute
    @test compute_heatcapacity(cp1_nd; random_name=1) ≈ 1.3368075000000002e22  # compute

    # Temperature-dependent heat capacity
    # dimensional
    T = 250.0:100:1250
    cp2 = T_HeatCapacity_Whittington()
    Cp = similar(T)
    @test isbits(cp2)
    args = (; T=T)
    compute_heatcapacity!(Cp, cp2, args)
    @test sum(Cp) ≈ 11667.035717418683

    # nondimensional
    cp2_nd = T_HeatCapacity_Whittington()
    cp2_nd = nondimensionalize(cp2_nd, CharUnits_GEO)
    T_nd = Float64.(T * K / CharUnits_GEO.Temperature)
    Cp_nd = similar(T)
    args = (; T=T_nd)
    compute_heatcapacity!(Cp_nd, cp2_nd, args)
    @test sum(Cp_nd) ≈ 1.4853886523631602e23

    # Dimensionalize again and double-check the results
    @test sum(abs.(ustrip.(Cp_nd * CharUnits_GEO.heatcapacity) - Cp)) < 1e-11

    # Test with arrays
    T_array = T * ones(10)'
    Cp_array = similar(T_array)
    compute_heatcapacity!(Cp_array, cp1, (;))
    @test Cp_array[1] ≈ 1050

    Cp_array = similar(T_array)
    compute_heatcapacity!(Cp_array, cp2, (; T=T_array))
    @test sum(Cp_array[:, 1]) ≈ 11667.035717418683

    T_array = T * ones(10)'
    Cp_array = zeros(size(T_array))
    compute_heatcapacity!(Cp_array, cp2, (; T=T_array))
    @test sum(Cp_array[:, 1]) ≈ 11667.035717418683

    # Check that it works if we give a phase array
    MatParam = Array{MaterialParams,1}(undef, 2)
    MatParam[1] = SetMaterialParams(;
        Name="Mantle", Phase=1, HeatCapacity=ConstantHeatCapacity()
    )

    MatParam[2] = SetMaterialParams(;
        Name="Crust", Phase=2, HeatCapacity=T_HeatCapacity_Whittington()
    )

    Mat_tup = Tuple(MatParam)

    Mat_tup1 = (
        SetMaterialParams(; Name="Mantle", Phase=1, HeatCapacity=ConstantHeatCapacity()),
        SetMaterialParams(;
            Name="Crust", Phase=2, HeatCapacity=ConstantHeatCapacity(; cp=1100J / kg / K)
        ),
    )

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    Phases[:, :, 20:end] .= 2

    Cp = zeros(size(Phases))
    T = ones(size(Phases)) * 1500
    P = zeros(size(Phases))

    args = (; T=T)
    compute_heatcapacity!(Cp, Mat_tup, Phases, args)    # computation routine w/out P (not used in most heat capacity formulations)
    @test sum(Cp[1, 1, :]) ≈ 121399.0486067196

    # check with array of constant properties (and no required input args)
    args1 = (;)
    compute_heatcapacity!(Cp, Mat_tup1, Phases, args1)    # computation routine w/out P (not used in most heat capacity formulations)
    @test sum(Cp[1, 1, :]) ≈ 109050.0

    num_alloc = @allocated compute_heatcapacity!(Cp, Mat_tup, Phases, args)
    @test sum(Cp[1, 1, :]) ≈ 121399.0486067196
    # @test num_alloc <= 32

    # test if we provide phase ratios
    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end
    compute_heatcapacity!(Cp, Mat_tup, PhaseRatio, args)
    num_alloc = @allocated compute_heatcapacity!(Cp, Mat_tup, PhaseRatio, args)
    @test sum(Cp[1, 1, :]) ≈ 121399.0486067196
    # @test num_alloc <= 32

    # -----------------------

    # Conductivity ----------

    # Constant

    # Constant conductivity
    cond = ConstantConductivity()
    @test isbits(cond)
    @test NumValue(cond.k) == 3.0
    @test cond.k.unit == u"W" / K / m

    cond = nondimensionalize(cond, CharUnits_GEO)
    @test NumValue(cond.k) ≈ 3.8194500000000007

    @test compute_conductivity(cond; T=100.0) ≈ 3.8194500000000007 # compute

    # Temperature-dependent conductivity
    # dimensional
    T = Vector{Float64}(250:100:1250)
    cond2 = T_Conductivity_Whittington()
    k = compute_conductivity(cond2, T)
    @test isbits(cond2)
    @test sum(k) ≈ 27.503366436682285

    # nondimensional
    cond2_nd = T_Conductivity_Whittington()
    cond2_nd = nondimensionalize(cond2_nd, CharUnits_GEO)
    T_nd = Float64.(ustrip.(T / CharUnits_GEO.Temperature))
    k_nd = compute_conductivity(cond2_nd, T_nd)
    @test sum(k_nd) ≈ 35.01591097886205

    k1 = zeros(size(T))
    args = (; T=ustrip.(T))
    compute_conductivity!(k1, cond2, args)
    @test sum(abs.(k - k1)) < 1e-13

    # Dimensionalize again and double-check the results
    @test sum(abs.(ustrip.(k_nd * CharUnits_GEO.conductivity) - k)) < 1e-11

    # Temperature-dependent parameterised conductivity
    # dimensional
    T = Vector{Float64}(250:100:1250)
    cond2 = T_Conductivity_Whittington_parameterised()
    k = compute_conductivity(cond2, T)
    @test isbits(cond2)
    @test sum(k) ≈ 27.553653387829254

    # nondimensional
    cond2_nd = T_Conductivity_Whittington_parameterised()
    cond2_nd = nondimensionalize(cond2_nd, CharUnits_GEO)
    T_nd = Float64.(ustrip.(T / CharUnits_GEO.Temperature))
    k_nd = compute_conductivity(cond2_nd, T_nd)
    @test sum(k_nd) ≈ 35.079933810714806

    # Dimensionalize again and double-check the results
    @test sum(abs.(ustrip.(k_nd * CharUnits_GEO.conductivity) - k)) < 1e-11

    # Check if we use arrays
    T_array = ustrip.(T) * ones(100)'
    k_array = copy(T_array)
    P_array = copy(T_array)

    args = (T=T_array, P=P_array)
    compute_conductivity!(k_array, cond, args)
    @test k_array[1] ≈ 3.8194500000000007

    compute_conductivity!(k_array, cond2, (; T=T_array))
    @test sum(k_array) ≈ 2755.3653387829254

    k_TP = Set_TP_Conductivity("LowerCrust")
    compute_conductivity!(k_array, k_TP, (P=P_array, T=T_array))
    @test sum(k_array) ≈ 2055.7129327367625

    # Check that it works if we give a phase array
    MatParam = Array{MaterialParams,1}(undef, 3)
    MatParam[1] = SetMaterialParams(;
        Name="Mantle", Phase=1, Conductivity=ConstantConductivity()
    )

    MatParam[2] = SetMaterialParams(;
        Name="Crust", Phase=2, Conductivity=T_Conductivity_Whittington()
    )

    MatParam[3] = SetMaterialParams(;
        Name="MantleLithosphere", Phase=3, Conductivity=Set_TP_Conductivity("Mantle")
    )

    Mat_tup = Tuple(MatParam)

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    Phases[:, :, 20:end] .= 2
    Phases[:, :, 60:end] .= 3

    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end

    k = zeros(size(Phases))
    T = ones(size(Phases)) * 1500
    P = zeros(size(Phases))
    args = (P=P, T=T)

    compute_conductivity!(k, Mat_tup, Phases, args)
    @test sum(k) ≈ 1.9216938849389635e6
    # num_alloc = @allocated compute_conductivity!(k, Mat_tup, Phases, args)
    # @test num_alloc <= 32

    compute_conductivity!(k, Mat_tup, PhaseRatio, args)
    @test sum(k) ≈ 1.9216938849389635e6
    # num_alloc = @allocated compute_conductivity!(k, Mat_tup, PhaseRatio, args)
    # @test num_alloc <= 32

    ######

    # TP-dependent conductivity for different predefines cases
    T = Vector{Float64}(250:100:1250)
    P = 1e6 * ones(size(T)) / ustrip(uconvert(Pa, 1MPa))  # must be in MPa!
    List = ["LowerCrust" "Mantle" "OceanicCrust" "UpperCrust"]
    Sol_kT = [20.55712932736763 28.700405819019323 20.55712932736763 19.940302462417037]
    for i in 1:length(List)
        k_TP = Set_TP_Conductivity(List[i])
        k = compute_conductivity(k_TP, P, T)         # note that P must be in MPa
        @test sum(k) ≈ Sol_kT[i]

        k_TP_nd = deepcopy(k_TP)
        k_TP_nd = nondimensionalize(k_TP_nd, CharUnits_GEO)
        T_nd = Float64.(ustrip.(T / CharUnits_GEO.Temperature))
        P_nd = Float64.(ustrip(P / CharUnits_GEO.stress))
        k_nd = compute_conductivity(k_TP_nd, P_nd, T_nd)

        @test ustrip(sum(abs.(ustrip.(k_nd * CharUnits_GEO.conductivity) - k))) < 1e-11
    end

    T = [200.0 300.0; 400.0 500.0]
    k1 = compute_conductivity(cond2, T)
    # -----------------------

    # Latent heat -----------
    a = ConstantLatentHeat()
    Q_L = compute_latent_heat(a)
    @test isbits(a)
    @test Q_L == 400

    a = nondimensionalize(a, CharUnits_GEO)
    Q_L = compute_latent_heat(a)
    @test Q_L ≈ 4e21

    # Check that it works if we give a phase array (including with an empty field)
    Mat_tup = (
        SetMaterialParams(; Name="Mantle", Phase=1, LatentHeat=ConstantLatentHeat()),
        SetMaterialParams(;
            Name="Crust", Phase=2, LatentHeat=ConstantLatentHeat(; Q_L=153kJ / kg)
        ),
        SetMaterialParams(; Name="MantleLithosphere", Phase=3),
    )

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    Phases[:, :, 20:end] .= 2
    Phases[:, :, 60:end] .= 3

    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end

    Hl = zeros(size(Phases))
    z = ones(size(Phases)) * 10e3
    args = (;)

    compute_latent_heat!(Hl, Mat_tup, Phases, args)
    @test minimum(Hl) ≈ 0.0
    @test maximum(Hl) ≈ 400
    @test Hl[50, 50, 50] ≈ 153.0

    compute_latent_heat!(Hl, Mat_tup, PhaseRatio, args)
    @test sum(Hl) ≈ 1.372e8

    # -----------------------

    # Radioactive heat ------
    a = ConstantRadioactiveHeat()
    H_r = compute_radioactive_heat(a)
    @test isbits(a)
    @test H_r ≈ 1.0e-6

    a = nondimensionalize(a, CharUnits_GEO)
    H_r = compute_radioactive_heat(a)
    @test H_r ≈ 0.1

    # depth-dependent radioactive heating:
    a = ExpDepthDependentRadioactiveHeat()
    z = 10e3
    H_r = compute_radioactive_heat(a, (z = z))
    @test isbits(a)
    @test H_r ≈ 3.678794411714423e-7

    Nx, Nz = 101, 101
    z = ones(Nx, Nz) * 10e3
    Hr = zero(z)
    compute_radioactive_heat!(Hr, a, (; z=z))

    @test sum(H_r) ≈ 3.678794411714423e-7

    # Check that it works if we give a phase array (including with an empty field)
    Mat_tup = (
        SetMaterialParams(;
            Name="Mantle", Phase=1, RadioactiveHeat=ConstantRadioactiveHeat()
        ),
        SetMaterialParams(;
            Name="Crust", Phase=2, RadioactiveHeat=ExpDepthDependentRadioactiveHeat()
        ),
        SetMaterialParams(; Name="MantleLithosphere", Phase=3),
    )
    Mat_tup1 = (
        SetMaterialParams(;
            Name="Mantle", Phase=1, RadioactiveHeat=ConstantRadioactiveHeat()
        ),
        SetMaterialParams(;
            Name="Crust", Phase=2, RadioactiveHeat=ConstantRadioactiveHeat()
        ),
        SetMaterialParams(; Name="MantleLithosphere", Phase=3),
    )

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    Phases[:, :, 20:end] .= 2
    Phases[:, :, 60:end] .= 3

    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end

    Hr = zeros(size(Phases))
    z = ones(size(Phases)) * 10e3
    args = (z=z,)

    compute_radioactive_heat!(Hr, Mat_tup, Phases, args)
    @test minimum(Hr) ≈ 0.0
    @test maximum(Hr) ≈ 1e-6
    @test Hr[50, 50, 50] ≈ 3.678794411714423e-7

    args1 = (;)
    compute_radioactive_heat!(Hr, Mat_tup1, Phases, args1)
    @test Hr[50, 50, 50] ≈ 1e-6

    # num_alloc = @allocated compute_radioactive_heat!(Hr, Mat_tup, Phases, args)
    # @test num_alloc <= 32   # in the commandline this gives 0; while running the script not always

    compute_radioactive_heat!(Hr, Mat_tup, PhaseRatio, args)
    @test sum(Hr) ≈ 0.33715177646857664

    # -----------------------

    # Shear heating -------
    Χ = ConstantShearheating(1.0)
    @test isbits(Χ)

    # Define parameters as vectors
    τ    = (1, 2, 2, 3) .* 1e6
    ε    = (1.0, 0.1, 0.1, 1.0)
    ε_el = (0.01, 0.01, 0.01, 0.01)

    τ_2D = [1 2; 2 3] * 1e6
    ε_2D = [1 0.1; 0.1 1]
    ε_el_2D = [0.01 0.01; 0.01 0.01]

    # With elasticity
    H_s1 = compute_shearheating(Χ, τ, ε, ε_el)
    H_s2 = compute_shearheating(Χ, τ_2D, ε_2D, ε_el_2D)
    @test H_s1 ≈ 4.32e6
    @test H_s2 ≈ 4.32e6
 
    # No elasticity
    H_s3 = compute_shearheating(Χ, τ, ε)
    H_s4 = compute_shearheating(Χ, τ_2D, ε_2D)
    @test H_s3 ≈ 4.4e6
    @test H_s4 ≈ 4.4e6

    # symmetric tensors
    τ    = (1, 3, 2) .* 1e6
    ε    = (1.0, 1.0, 0.1)
    ε_el = (0.01, 0.01, 0.01)
    @test H_s1 == compute_shearheating(Χ, τ, ε, ε_el)
    @test H_s3 == compute_shearheating(Χ, τ, ε)
 
       
    # test in-place computation
    n = 12
    τ_xx = fill(1e6, n, n)
    τ_xy = fill(2e6, n, n)
    τ_yx = fill(3e6, n, n)
    τ_yy = fill(4e6, n, n)
    ε_xx = fill(1.0, n, n)
    ε_xy = fill(0.1, n, n)
    ε_yx = fill(0.1, n, n)
    ε_yy = fill(1.0, n, n)
    ε_el_xx = fill(0.01, n, n)
    ε_el_xy = fill(0.01, n, n)
    ε_el_yx = fill(0.01, n, n)
    ε_el_yy = fill(0.01, n, n)
    H_s = similar(τ_xx)
    τ = τ_xx, τ_xy, τ_yy, τ_yx
    ε = ε_xx, ε_xy, ε_yy, ε_yx
    ε_el = ε_el_xx, ε_el_xy, ε_el_yy, ε_el_yx
    compute_shearheating!(H_s, Χ, τ, ε, ε_el)
    @test all(x == 5.4e6 for x in H_s)

    compute_shearheating!(H_s, Χ, τ, ε)
    @test all(x == 5.5e6 for x in H_s)

    # Now in non-dimensional units
    τ = [1 2 3 4]
    ε = [1 0.1 0.1 1]
    ε_el = [0.01 0.01 0.01 0.01]

    τ_2D = [1 2; 3 4]
    ε_2D = [1 0.1; 0.1 1]
    ε_el_2D = [0.01 0.01; 0.01 0.01]
    Χ = nondimensionalize(Χ, CharUnits_GEO)
    info = param_info(Χ)

    H_s1 = compute_shearheating(Χ, τ, ε, ε_el)
    H_s2 = compute_shearheating(Χ, τ_2D, ε_2D, ε_el_2D)
    H_s3 = compute_shearheating(Χ, τ, ε)
    H_s4 = compute_shearheating(Χ, τ_2D, ε_2D)
    @test H_s1 ≈ 5.4
    @test H_s2 ≈ 5.4
    @test H_s3 ≈ 5.5
    @test H_s4 ≈ 5.5

    # test material structs    
    Mat_tup = (
        SetMaterialParams(;
            Name="Mantle", Phase=1, ShearHeat = ConstantShearheating(Χ=0.0NoUnits),
        ),
        SetMaterialParams(;
            Name="Crust", Phase=2, ShearHeat = ConstantShearheating(Χ=1.0NoUnits),
        ),
    )
    args = (τ=τ, ε=ε, ε_el=ε_el)
    @test compute_shearheating(Mat_tup[1], τ, ε, ε_el) == 0.0
    @test compute_shearheating(Mat_tup[2], τ, ε, ε_el) == 5.4
    @test compute_shearheating(Mat_tup, 1, τ, ε, ε_el) == 0.0
    @test compute_shearheating(Mat_tup, 2, τ, ε, ε_el) == 5.4

end

