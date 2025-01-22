using Test
using GeoParams
using StaticArrays
import ForwardDiff as FD
@testset "EnergyParameters.jl" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 10km)

    # Heat capacity ---------

    # Constant heat capacity
    Cp1 = ConstantHeatCapacity()
    info = param_info(Cp1)
    @test isbits(Cp1)
    @test Cp1.Cp.val == 1050.0
    @test GeoParams.get_Cp(Cp1) == 1050.0

    Cp1_nd = Cp1
    Cp1_nd = nondimensionalize(Cp1_nd, CharUnits_GEO)
    @test Cp1_nd.Cp.val ≈ 1.3368075000000002e22
    @test compute_heatcapacity(Cp1_nd) ≈ 1.3368075000000002e22  # compute
    @test compute_heatcapacity(Cp1_nd, (random_name = 1,)) ≈ 1.3368075000000002e22  # compute

    # Temperature-dependent heat capacity
    # dimensional
    T = 250.0:100:1250
    Cp2 = T_HeatCapacity_Whittington()
    Cp = similar(T)
    @test isbits(Cp2)
    args = (; T = T)
    compute_heatcapacity!(Cp, Cp2, args)
    @test sum(Cp) ≈ 11667.035717418683
    @test GeoParams.get_Tcutoff(Cp2) == 846.0

    FD.derivative(x -> compute_heatcapacity(Cp2, (; T = x)), T[1]) == 3.2721616015871584

    # nondimensional
    Cp2_nd = T_HeatCapacity_Whittington()
    Cp2_nd = nondimensionalize(Cp2_nd, CharUnits_GEO)
    T_nd = Float64.(T * K / CharUnits_GEO.Temperature)
    Cp_nd = similar(T)
    args = (; T = T_nd)
    compute_heatcapacity!(Cp_nd, Cp2_nd, args)
    @test sum(Cp_nd) ≈ 1.4853886523631602e23

    # Dimensionalize again and double-check the results
    @test sum(abs.(ustrip.(Cp_nd * CharUnits_GEO.heatcapacity) - Cp)) < 1.0e-11

    # heat capacity
    Cp3 = Vector_HeatCapacity(Cp = fill(1500.0, 100))
    index = fill(10, size(T))
    args = (; index = index)
    Cp = zero(T)
    @test compute_heatcapacity(Cp3, index = 10) == 1500.0

    # compute_heatcapacity!(Cp, Cp3, args)

    # Test with arrays
    T_array = T * ones(10)'
    Cp_array = similar(T_array)
    compute_heatcapacity!(Cp_array, Cp1, (;))
    @test Cp_array[1] ≈ 1050

    Cp_array = similar(T_array)
    compute_heatcapacity!(Cp_array, Cp2, (; T = T_array))
    @test sum(Cp_array[i, 1] for i in axes(Cp_array, 1)) ≈ 11667.035717418683

    T_array = T * ones(10)'
    Cp_array = zeros(size(T_array))
    compute_heatcapacity!(Cp_array, Cp2, (; T = T_array))
    @test sum(Cp_array[i, 1] for i in axes(Cp_array, 1)) ≈ 11667.035717418683

    # Check that it works if we give a phase array
    MatParam = Array{MaterialParams, 1}(undef, 3)
    MatParam[1] = SetMaterialParams(;
        Name = "Mantle", Phase = 1, HeatCapacity = ConstantHeatCapacity()
    )

    MatParam[2] = SetMaterialParams(;
        Name = "Crust", Phase = 2, HeatCapacity = T_HeatCapacity_Whittington()
    )

    MatParam[3] = SetMaterialParams(;
        Name = "Crust1", Phase = 3, HeatCapacity = Vector_HeatCapacity(Cp = fill(1500.0, 100))
    )

    Mat_tup = Tuple(MatParam)

    Mat_tup1 = (
        SetMaterialParams(; Name = "Mantle", Phase = 1, HeatCapacity = ConstantHeatCapacity()),
        SetMaterialParams(;
            Name = "Crust", Phase = 2, HeatCapacity = ConstantHeatCapacity(; Cp = 1100J / kg / K)
        ),
    )

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    @views Phases[:, :, 20:end] .= 2
    @views Phases[:, :, 50:end] .= 3

    Cp = zeros(size(Phases))
    T = fill(1500.0e0, size(Phases))
    P = zeros(size(Phases))

    args = (; T = T, index = fill(10, size(T)))
    compute_heatcapacity!(Cp, Mat_tup, Phases, args)    # computation routine w/out P (not used in most heat capacity formulations)
    @test sum(Cp[1, 1, k] for k in axes(Cp, 3)) ≈ 134023.72170619245

    # check with array of constant properties (and no required input args)
    args1 = (;)
    compute_heatcapacity!(Cp, Mat_tup1, Phases, args1)    # computation routine w/out P (not used in most heat capacity formulations)
    @test sum(Cp[1, 1, k] for k in axes(Cp, 3)) ≈ 52950.0

    num_alloc = @allocated compute_heatcapacity!(Cp, Mat_tup, Phases, args)
    @test sum(Cp[1, 1, k] for k in axes(Cp, 3)) ≈ 134023.72170619245

    @test num_alloc == 0

    # test if we provide phase ratios
    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end
    compute_heatcapacity!(Cp, Mat_tup, PhaseRatio, args)
    num_alloc = @allocated compute_heatcapacity!(Cp, Mat_tup, PhaseRatio, args)
    @test sum(Cp[1, 1, k] for k in axes(Cp, 3)) ≈ 134023.72170619245
    @test num_alloc == 0


    # Test latent heat based heat capacity
    CharUnits_GEO = GEO_units(; viscosity = 1.0e19, length = 10km)
    x_D = Latent_HeatCapacity(Q_L = 500.0e3 * J / kg)
    x_D1 = Latent_HeatCapacity(Cp = ConstantHeatCapacity())
    x_ND = nondimensionalize(x_D, CharUnits_GEO)
    x_ND1 = nondimensionalize(x_D1, CharUnits_GEO)

    @test isbits(x_D)
    @test isbits(x_D1)
    @test isbits(x_ND)
    @test isdimensional(x_D) == true
    @test isdimensional(x_D1) == true
    @test isdimensional(x_ND) == false
    @test isdimensional(x_ND1) == false

    dϕdT = 0.1
    dϕdT_ND = nondimensionalize(dϕdT / K, CharUnits_GEO)
    args = (; dϕdT = dϕdT, T = 300.0 + 273)
    args_ND = (; dϕdT = dϕdT_ND, T = 300.0 + 273)
    @test compute_heatcapacity(x_D, args) == 1050 + 500.0e3 * dϕdT

    @test compute_heatcapacity(x_D1, args) == 1050 + 400.0e3 * dϕdT

    x_ND = nondimensionalize(x_D, CharUnits_GEO)
    Cp_nd = compute_heatcapacity(x_ND, args_ND)
    @test compute_heatcapacity(x_D, args) ≈ dimensionalize(Cp_nd, J / kg / K, CharUnits_GEO).val

    x_ND1 = nondimensionalize(x_D1, CharUnits_GEO)
    Cp_nd1 = compute_heatcapacity(x_ND1, args_ND)
    @test compute_heatcapacity(x_D1, args) ≈ dimensionalize(Cp_nd1, J / kg / K, CharUnits_GEO).val

    #Temperature-dependent latent heat based heat capacity

    T = 300.0 + 273
    dϕdT = 0.1
    args = (; T = T, dϕdT = dϕdT)

    x_D = Latent_HeatCapacity(Cp = T_HeatCapacity_Whittington())

    @test isbits(x_D)
    @test isdimensional(x_D) == true

    @test compute_heatcapacity(x_D, args) == 41052.29268922852
    @test FD.derivative(x -> compute_heatcapacity(x_D, (; T = x, dϕdT = dϕdT)), T) == 0.6260890173995907
    @test FD.derivative(x -> compute_heatcapacity(x_D, (; T = T, dϕdT = x)), dϕdT) == 400000.0

    dϕdT_ND = nondimensionalize(dϕdT / K, CharUnits_GEO)
    args_ND = (; T = ustrip.(T * K / CharUnits_GEO.Temperature), dϕdT = dϕdT_ND)
    x_ND = nondimensionalize(x_D, CharUnits_GEO)

    @test isdimensional(x_ND) == false
    @test isbits(x_ND)

    Cp_nd = compute_heatcapacity(x_ND, args_ND)
    @test compute_heatcapacity(x_D, args) ≈ dimensionalize(Cp_nd, J / kg / K, CharUnits_GEO).val
    # -----------------------

    # Conductivity ----------

    # Constant

    # Constant conductivity
    cond = ConstantConductivity()
    @test isbits(cond)
    @test NumValue(cond.k) == 3.0
    @test cond.k.unit == u"W" / K / m
    @test GeoParams.get_k(cond) == 3.0

    cond = nondimensionalize(cond, CharUnits_GEO)
    @test NumValue(cond.k) ≈ 3.8194500000000007

    @test compute_conductivity(cond; T = 100.0) ≈ 3.8194500000000007 # compute

    # Temperature-dependent conductivity
    # dimensional
    T = collect(250.0e0:100:1250.0e0)
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
    args = (; T = ustrip.(T))
    compute_conductivity!(k1, cond2, args)
    @test sum(abs.(k - k1)) < 1.0e-13

    # Dimensionalize again and double-check the results
    @test sum(abs.(ustrip.(k_nd * CharUnits_GEO.conductivity) - k)) < 1.0e-11

    # Temperature-dependent parameterised conductivity
    # dimensional
    T = collect(250.0e0:100:1250.0e0)
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
    @test sum(abs.(ustrip.(k_nd * CharUnits_GEO.conductivity) - k)) < 1.0e-11

    # Check if we use arrays
    T_array = ustrip.(T) * ones(100)'
    k_array = copy(T_array)
    P_array = copy(T_array)

    args = (T = T_array, P = P_array)
    compute_conductivity!(k_array, cond, args)
    @test k_array[1] ≈ 3.8194500000000007

    compute_conductivity!(k_array, cond2, (; T = T_array))
    @test sum(k_array) ≈ 2755.3653387829254

    k_TP = Set_TP_Conductivity("LowerCrust")
    compute_conductivity!(k_array, k_TP, (P = P_array, T = T_array))
    @test sum(k_array) ≈ 2055.7129327367625

    # Check that it works if we give a phase array
    MatParam = Array{MaterialParams, 1}(undef, 3)
    MatParam[1] = SetMaterialParams(;
        Name = "Mantle", Phase = 1, Conductivity = ConstantConductivity()
    )

    MatParam[2] = SetMaterialParams(;
        Name = "Crust", Phase = 2, Conductivity = T_Conductivity_Whittington()
    )

    MatParam[3] = SetMaterialParams(;
        Name = "MantleLithosphere", Phase = 3, Conductivity = Set_TP_Conductivity("Mantle")
    )

    Mat_tup = Tuple(MatParam)

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    @views Phases[:, :, 20:end] .= 2
    @views Phases[:, :, 60:end] .= 3

    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end

    k = zeros(size(Phases))
    T = fill(1500.0e0, size(Phases))
    P = zeros(size(Phases))
    args = (P = P, T = T)

    compute_conductivity!(k, Mat_tup, Phases, args)
    @test sum(k) ≈ 1.9216938849389635e6
    num_alloc = @allocated compute_conductivity!(k, Mat_tup, Phases, args)
    @test num_alloc == 0

    compute_conductivity!(k, Mat_tup, PhaseRatio, args)
    @test sum(k) ≈ 1.9216938849389635e6
    num_alloc = @allocated compute_conductivity!(k, Mat_tup, PhaseRatio, args)
    @test num_alloc == 0

    ######

    # TP-dependent conductivity for different predefines cases
    T = Vector{Float64}(250:100:1250)
    P = 1.0e6 * ones(size(T)) / ustrip(uconvert(Pa, 1MPa))  # must be in MPa!
    List = "LowerCrust", "Mantle", "OceanicCrust", "UpperCrust"
    Sol_kT = 20.55712932736763, 28.700405819019323, 20.55712932736763, 19.940302462417037
    for i in eachindex(List)
        k_TP = Set_TP_Conductivity(List[i])
        k = compute_conductivity(k_TP, P, T)         # note that P must be in MPa
        @test sum(k) ≈ Sol_kT[i]

        k_TP_nd = deepcopy(k_TP)
        k_TP_nd = nondimensionalize(k_TP_nd, CharUnits_GEO)
        T_nd = Float64.(ustrip.(T / CharUnits_GEO.Temperature))
        P_nd = Float64.(ustrip(P / CharUnits_GEO.stress))
        k_nd = compute_conductivity(k_TP_nd, P_nd, T_nd)

        @test ustrip(sum(abs.(ustrip.(k_nd * CharUnits_GEO.conductivity) - k))) < 1.0e-11
    end

    T = [200.0 300.0; 400.0 500.0]
    k1 = compute_conductivity(cond2, T)
    @test FD.derivative(x -> compute_conductivity(cond2, (; T = x)), T[1]) == -0.007109905535
    # -----------------------

    # Latent heat -----------
    a = ConstantLatentHeat(Q_L = 400.0e3J / kg)
    Q_L = compute_latent_heat(a)
    @test isbits(a)
    @test Q_L == 400.0e3

    a = nondimensionalize(a, CharUnits_GEO)
    Q_L = compute_latent_heat(a)
    @test Q_L ≈ 4.0e21

    # Check that it works if we give a phase array (including with an empty field)
    Mat_tup = (
        SetMaterialParams(; Name = "Mantle", Phase = 1, LatentHeat = ConstantLatentHeat(Q_L = 400.0e3J / kg)),
        SetMaterialParams(;
            Name = "Crust", Phase = 2, LatentHeat = ConstantLatentHeat(; Q_L = 153.0e3J / kg)
        ),
        SetMaterialParams(; Name = "MantleLithosphere", Phase = 3),
    )

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    @views Phases[:, :, 20:end] .= 2
    @views Phases[:, :, 60:end] .= 3

    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end

    Hl = zeros(size(Phases))
    z = fill(10.0e3, size(Phases))
    args = (;)

    compute_latent_heat!(Hl, Mat_tup, Phases, args)
    @test minimum(Hl) ≈ 0.0
    @test maximum(Hl) ≈ 400.0e3
    @test Hl[50, 50, 50] ≈ 153.0e3

    compute_latent_heat!(Hl, Mat_tup, PhaseRatio, args)
    @test sum(Hl) ≈ 1.372e11

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
    z = 10.0e3
    H_r = compute_radioactive_heat(a, (z = z))
    @test isbits(a)
    @test H_r ≈ 3.678794411714423e-7

    Nx, Nz = 101, 101
    z = ones(Nx, Nz) * 10.0e3
    Hr = zero(z)
    compute_radioactive_heat!(Hr, a, (; z = z))

    @test sum(H_r) ≈ 3.678794411714423e-7

    # Check that it works if we give a phase array (including with an empty field)
    Mat_tup = (
        SetMaterialParams(;
            Name = "Mantle", Phase = 1, RadioactiveHeat = ConstantRadioactiveHeat()
        ),
        SetMaterialParams(;
            Name = "Crust", Phase = 2, RadioactiveHeat = ExpDepthDependentRadioactiveHeat()
        ),
        SetMaterialParams(; Name = "MantleLithosphere", Phase = 3),
    )
    Mat_tup1 = (
        SetMaterialParams(;
            Name = "Mantle", Phase = 1, RadioactiveHeat = ConstantRadioactiveHeat()
        ),
        SetMaterialParams(;
            Name = "Crust", Phase = 2, RadioactiveHeat = ConstantRadioactiveHeat()
        ),
        SetMaterialParams(; Name = "MantleLithosphere", Phase = 3),
    )

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    @views Phases[:, :, 20:end] .= 2
    @views Phases[:, :, 60:end] .= 3

    PhaseRatio = zeros(n, n, n, 3)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end

    Hr = zeros(size(Phases))
    z = fill(10.0e3, size(Phases))
    args = (z = z,)

    compute_radioactive_heat!(Hr, Mat_tup, Phases, args)
    @test minimum(Hr) ≈ 0.0
    @test maximum(Hr) ≈ 1.0e-6
    @test Hr[50, 50, 50] ≈ 3.678794411714423e-7

    args1 = (;)
    compute_radioactive_heat!(Hr, Mat_tup1, Phases, args1)
    @test Hr[50, 50, 50] ≈ 1.0e-6

    num_alloc = @allocated compute_radioactive_heat!(Hr, Mat_tup, Phases, args)
    @test num_alloc == 0

    compute_radioactive_heat!(Hr, Mat_tup, PhaseRatio, args)
    @test sum(Hr) ≈ 0.33715177646857664

    # -----------------------

    # Shear heating -------
    Χ = ConstantShearheating(1.0)
    @test isbits(Χ)

    # Define parameters as vectors
    τ = (1, 2, 2, 3) .* 1.0e6
    ε = (1.0, 0.1, 0.1, 1.0)
    ε_el = (0.01, 0.01, 0.01, 0.01)

    τ_2D = [1 2; 2 3] * 1.0e6
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
    τ = (1, 3, 2) .* 1.0e6
    ε = (1.0, 1.0, 0.1)
    ε_el = (0.01, 0.01, 0.01)
    @test H_s1 == compute_shearheating(Χ, τ, ε, ε_el)
    @test H_s3 == compute_shearheating(Χ, τ, ε)

    # test in-place computation
    n = 12
    τ_xx = fill(1.0e6, n, n)
    τ_xy = fill(2.0e6, n, n)
    τ_yx = fill(3.0e6, n, n)
    τ_yy = fill(4.0e6, n, n)
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
    @test H_s1 ≈ H_s2 ≈ 5.4
    @test H_s3 ≈ H_s4 ≈ 5.5

    # test material structs
    rheology = (
        SetMaterialParams(;
            Name = "Mantle", Phase = 1, ShearHeat = ConstantShearheating(Χ = 0.0NoUnits),
        ),
        SetMaterialParams(;
            Name = "Crust", Phase = 2, ShearHeat = ConstantShearheating(Χ = 1.0NoUnits),
        ),
    )
    @test iszero(compute_shearheating(rheology[1], τ, ε, ε_el))
    @test compute_shearheating(rheology[2], τ, ε, ε_el) == 5.4
    @test iszero(compute_shearheating(rheology, 1, τ, ε, ε_el))
    @test compute_shearheating(rheology, 2, τ, ε, ε_el) == 5.4

    # Test with phase ratios
    phase = SA[0.5, 0.5] # static array
    @test compute_shearheating(rheology, phase, τ, ε, ε_el) == 2.7
    phase = (0.5, 0.5) # tuple
    @test compute_shearheating(rheology, phase, τ, ε, ε_el) == 2.7

    # 3D shear heating tests
    τ_3D = [
        1.0e0 2.0e0 3.0e0
        2.0e0 5.0e0 6.0e0
        3.0e0 6.0e0 9.0e0
    ]
    ε_3D = [
        1.0e0  1.0e-1 1.0e-1
        1.0e-1 1.0e0  1.0e-1
        1.0e-1 1.0e-1 1.0e0
    ]
    ε_el_3D = fill(0.01, 3, 3)

    τ_3D_voigt = 1.0e0, 5.0e0, 9.0e0, 6.0e0, 3.0e0, 2.0e0
    ε_3D_voigt = 1.0e0, 1.0e0, 1.0e0, 0.1, 0.1, 0.1
    ε_el_3D_voigt = ntuple(_ -> 1.0e-2, Val(6))

    H_s5 = compute_shearheating(Χ, τ_3D, ε_3D, ε_el_3D)
    H_s6 = compute_shearheating(Χ, τ_3D_voigt, ε_3D_voigt, ε_el_3D_voigt)
    H_s7 = compute_shearheating(Χ, τ_3D, ε_3D)
    H_s8 = compute_shearheating(Χ, τ_3D_voigt, ε_3D_voigt)
    @test H_s5 ≈ H_s6 ≈ 16.83
    @test H_s7 ≈ H_s7 ≈ 17.2

    phase = SA[0.5, 0.5] # static array
    H_s9 = compute_shearheating(rheology, phase, τ_3D, ε_3D, ε_el_3D)
    H_s10 = compute_shearheating(rheology, 1, τ_3D, ε_3D, ε_el_3D)
    H_s11 = compute_shearheating(rheology[1], τ_3D, ε_3D, ε_el_3D)
    H_s12 = compute_shearheating(rheology, phase, τ_3D_voigt, ε_3D_voigt, ε_el_3D_voigt)
    H_s13 = compute_shearheating(rheology, 1, τ_3D_voigt, ε_3D_voigt, ε_el_3D_voigt)
    H_s14 = compute_shearheating(rheology[1], τ_3D_voigt, ε_3D_voigt, ε_el_3D_voigt)
    @test H_s9 ≈ H_s12 ≈ 8.415
    @test all(iszero, (H_s10, H_s11, H_s13, H_s14))
end
