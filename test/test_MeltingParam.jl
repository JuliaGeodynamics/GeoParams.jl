using Test
using LinearAlgebra
using GeoParams
using StaticArrays

@testset "MeltingParam.jl" begin

    #Make sure structure is isbits
    x = MeltingParam_Caricchi()
    @test isbits(x)

    # This tests the various melting parameterizations
    CharUnits_GEO = GEO_units(; viscosity=1e19, length=10km)

    T = collect(250:100:1250) * K .+ 273.15K
    T_nd = Float64.(T / CharUnits_GEO.Temperature)

    # Caricchi parameterization [in ND numbers, which is anyways the typical use case]
    p = MeltingParam_Caricchi()
    phi_dim = zeros(size(T))
    args = (; T=ustrip.(T))
    compute_meltfraction!(phi_dim, p, args)

    phi_dim1 = zeros(size(phi_dim))
    compute_meltfraction!(phi_dim1, p, args) # in-place routine

    p_nd = p
    p_nd = nondimensionalize(p_nd, CharUnits_GEO)
    phi_nd = zeros(size(T))
    args = (; T=T_nd)
    compute_meltfraction!(phi_nd, p_nd, args)

    # Do this computation manually, using the actual expression of Caricchi
    T_C = collect(250:100:1250)    # in celsius
    Phi_solid = 1.0 .- 1.0 ./ (1.0 .+ exp.((800.0 .- T_C) ./ 23.0))
    Phi_anal = 1.0 .- Phi_solid

    @test sum(phi_dim - Phi_anal) < 1e-12
    @test sum(phi_dim1 - Phi_anal) < 1e-12
    @test sum(phi_nd - Phi_anal) < 1e-12

    # test derivative vs T
    dϕdT_dim = zeros(size(T))
    args = (; T=ustrip.(T))
    compute_dϕdT!(dϕdT_dim, p, args)
    @test sum(dϕdT_dim) ≈ 0.008102237679214096

    #------------------------------
    # 5th order polynomial
    p = MeltingParam_5thOrder()
    compute_meltfraction!(phi_dim, p, args)
    @test sum(phi_dim) ≈ 4.708427909521561

    # experimental data to create the fit
    data = [
        1115 100
        1050 90
        969 85
        932 80
        907 70
        880 60
        850 54
        825 51.7
        800 52.9
        775 46.3
        750 44.9
        725 29.9
        700 14.9
        690 0
    ]

    data[:, 2] = data[:, 2] / 100
    Tdata = data[:, 1] .+ 273.15
    phi = zeros(size(Tdata))
    args = (; T=ustrip.(Tdata))
    compute_meltfraction!(phi, p, args)

    @test norm(data[:, 2] - phi) ≈ 0.07151515017819135

    # test derivative vs T
    dϕdT_dim = zeros(size(T))
    args = (; T=ustrip.(T))
    compute_dϕdT!(dϕdT_dim, p, args)
    @test sum(dϕdT_dim) ≈ 0.006484458453421382
    #------------------------------

    #------------------------------
    # 4th order polynomial
    p = MeltingParam_4thOrder()
    compute_meltfraction!(phi_dim, p, args)
    @test sum(phi_dim) ≈ 4.853749635538406

    # experimental data to create the fit
    data = [
        1000 100
        990 100
        975 93
        950 89.2
        925 76.3
        900 69.6
        875 59
        850 54
        825 51.7
        800 52.9
        775 46.3
        750 44.9
        725 29.9
        700 14.9
        690 0
    ]

    data[:, 2] = data[:, 2] / 100
    Tdata = data[:, 1] .+ 273.15
    phi = zeros(size(Tdata))
    args = (; T=ustrip.(Tdata))
    compute_meltfraction!(phi, p, args)

    @test norm(data[:, 2] - phi) ≈ 0.0678052542705406

    # test derivative vs T
    dϕdT_dim = zeros(size(T))
    args = (; T=ustrip.(T))
    compute_dϕdT!(dϕdT_dim, p, args)
    @test sum(dϕdT_dim) ≈ 0.00830985782591842
    #------------------------------

    #------------------------------
    # Quadratic parameterisation
    p = MeltingParam_Quadratic()
    compute_meltfraction!(phi_dim, p, args)
    @test sum(phi_dim) ≈ 5.0894901144641

    dϕdT_dim = zeros(size(T))
    compute_dϕdT!(dϕdT_dim, p, args)
    @test sum(dϕdT_dim) ≈ 0.009365244536940681
    #------------------------------

    #------------------------------
    # Assimilation parameterisation
    p = MeltingParam_Assimilation()
    compute_meltfraction!(phi_dim, p, args)
    @test sum(phi_dim) ≈ 4.995
    dϕdT_dim = zeros(size(T))
    compute_dϕdT!(dϕdT_dim, p, args)
    @test sum(abs.(dϕdT_dim)) ≈ 0.004605170185988078
    #------------------------------

    #------------------------------
    # Melnik parameterisation
    p = MeltingParam_Melnik()
    compute_meltfraction!(phi_dim, p, args)
    @test sum(phi_dim) ≈ 2.546
    dϕdT_dim = zeros(size(T))
    compute_dϕdT!(dϕdT_dim, p, args)
    @test sum(abs.(dϕdT_dim)) ≈ 0.009842552900718364
    #------------------------------

    # Test computation of melt parameterization for the whole computational domain, using arrays 
    MatParam = Vector{MaterialParams}(undef, 4)
    MatParam[1] = SetMaterialParams(;
        Name="Mantle", Phase=1, Melting=PerpleX_LaMEM_Diagram("test_data/Peridotite_dry.in")
    )

    MatParam[2] = SetMaterialParams(;
        Name="Crust", Phase=2, Melting=MeltingParam_Caricchi()
    )

    # No melting parameterization for this phase
    MatParam[3] = SetMaterialParams(;
        Name="UpperCrust", Phase=3, Melting=MeltingParam_5thOrder(), Density=PT_Density()
    )

    # No melting parameterization for this phase
    MatParam[4] = SetMaterialParams(; Name="LowerCrust", Phase=4, Density=PT_Density())

    Mat_tup = Tuple(MatParam)

    # test computing material properties
    n = 100
    Phases = ones(Int64, n, n, n)
    Phases[:, :, 20:end] .= 2
    Phases[:, :, 80:end] .= 3
    Phases[:, :, 90:end] .= 4

    ϕ = zeros(size(Phases))
    dϕdT = zeros(size(Phases))
    T = ones(size(Phases)) * 1500
    P = ones(size(Phases)) * 10
    args = (P=P, T=T)
    compute_meltfraction!(ϕ, Mat_tup, Phases, args) #allocations coming from computing meltfraction using PhaseDiagram_LookupTable
    @test sum(ϕ) / n^3 ≈ 0.7484802337240443

    compute_dϕdT!(dϕdT, Mat_tup, Phases, args)
    @test sum(dϕdT) / n^3 ≈ 4.695945929818635e-5

    # test computing material properties when we have PhaseRatios, instead of Phase numbers
    PhaseRatio = zeros(n, n, n, 4)
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i, iz)
        PhaseRatio[I] = 1.0
    end

    compute_meltfraction!(ϕ, Mat_tup, PhaseRatio, args)
    @test sum(ϕ) / n^3 ≈ 0.7484802337240443

    compute_dϕdT!(dϕdT, Mat_tup, PhaseRatio, args)
    @test sum(dϕdT) / n^3 ≈ 4.695945929818635e-5

    # Test smoothening of the melting curves:
    p = SmoothMelting(; p=MeltingParam_5thOrder())
    @test isbits(p)
    T = collect(250:100:1250) * K .+ 273.15K
    phi_dim = zeros(size(T))
    args = (; T=ustrip.(T))
    compute_meltfraction!(phi_dim, p, args)
    @test sum(phi_dim) ≈ 4.7084279086574226

    dϕdT = zeros(size(T))
    compute_dϕdT!(dϕdT, p, args)
    @test sum(dϕdT) ≈ 0.006484456307540648

    # try non-dimensionalisation
    p_nd = nondimensionalize(p, CharUnits_GEO)
    @test isbits(p)
    @test isdimensional(p_nd.p.a) == false
    @test p_nd.p.a ≈ 6968.639721576996

    p1 = dimensionalize(p_nd, CharUnits_GEO)
    @test isdimensional(p1.p.a) == true
    @test Value(p1.p.a) ≈ Value(p.p.a)

    #Mat_tup = ( SetMaterialParams( Name="Mantle", Phase=1, Melting=SmoothMelting(MeltingParam_4thOrder())),
    #            SetMaterialParams( Name="Crust", Phase=2, Melting=MeltingParam_5thOrder()),
    #            SetMaterialParams( Name="UpperCrust", Phase=3, Melting=SmoothMelting(MeltingParam_5thOrder()), Density=PT_Density()), 
    #            SetMaterialParams( Name="LowerCrust", Phase=4, Density=PT_Density())
    #                )
    Mat_tup = (
        SetMaterialParams(;
            Name="Mantle", Phase=1, Melting=SmoothMelting(MeltingParam_4thOrder())
        ),
        SetMaterialParams(; Name="Crust", Phase=2, Melting=MeltingParam_5thOrder()),
        SetMaterialParams(;
            Name="UpperCrust",
            Phase=3,
            Melting=SmoothMelting(MeltingParam_5thOrder()),
            Density=PT_Density(),
        ),
        SetMaterialParams(; Name="LowerCrust", Phase=4, Density=PT_Density()),
    )

    ϕ = zeros(size(Phases))
    dϕdT = zeros(size(Phases))
    T = ones(size(Phases)) * (800 + 273.15)
    P = ones(size(Phases)) * 10
    args = (P=P, T=T)
    compute_meltfraction!(ϕ, Mat_tup, Phases, args) #allocation free
    @test sum(ϕ) / n^3 ≈ 0.4409766063389861

    compute_dϕdT!(dϕdT, Mat_tup, Phases, args) #allocation free
    @test sum(dϕdT) / n^3 ≈ 0.0006838372430250584

    # test PhaseRatio and StaticArrays PhaseRatios as input
    args = (P=0.0, T=1000.0 + 273.15)
    PhaseRatio = (0.25, 0.25, 0.25, 0.25)
    @test 0.6991003705903673 ≈ compute_meltfraction_ratio(PhaseRatio, Mat_tup, args) 

    SvPhaseRatio = SA[0.25,0.25,0.25,0.25]
    @test 0.6991003705903673 ≈ compute_meltfraction_ratio(SvPhaseRatio, Mat_tup, args)

end
