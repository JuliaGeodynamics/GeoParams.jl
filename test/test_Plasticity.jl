using Test
using GeoParams
    
@testset "Plasticity.jl" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO   =   GEO_units(viscosity=1e19, length=10km);
            
    # DruckerPrager ---------

    # Constant heat capacity
    p        =  DruckerPrager()
    info     =  param_info(p)
    @test isbits(p)
    @test NumValue(p.ϕ) == 30

    p_nd = p;
    p_nd = nondimensionalize(p_nd,CharUnits_GEO)
    @test p_nd.C.val ≈ 1

    # Compute with dimensional units
    τII = 20e6;
    P = 1e6;
    args = (P=P, τII=τII)
    @test compute_yieldfunction(p,args) ≈ 1.0839745962155614e7  # compute
   
    
  #=

    # Test with arrays
    T_array     =  T*ones(10)'
    Cp_array    =  similar(T_array)
    compute_heatcapacity!(Cp_array, cp1, (;))
    @test Cp_array[1] ≈ 1.3368075000000002e22

    Cp_array    =  similar(T_array)
    compute_heatcapacity!(Cp_array, cp2, (;T=T_array))
    @test sum(Cp_array[:,1]) ≈ 11667.035717418683

    T_array     =  T*ones(10)'
    Cp_array    =  zeros(size(T_array))
    compute_heatcapacity!(Cp_array, cp2, (;T=T_array))
    @test sum(Cp_array[:,1]) ≈ 11667.035717418683


    # Check that it works if we give a phase array
    MatParam    =   Array{MaterialParams, 1}(undef, 2);
    MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                        HeatCapacity  = ConstantHeatCapacity());

    MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                        HeatCapacity  = T_HeatCapacity_Whittington());

    Mat_tup = Tuple(MatParam)

    Mat_tup1 =   (  SetMaterialParams(Name="Mantle", Phase=1,
                        HeatCapacity  = ConstantHeatCapacity()),
                    SetMaterialParams(Name="Crust", Phase=2,
                        HeatCapacity  = ConstantHeatCapacity(cp=1100J/kg/K))
                 );


    # test computing material properties
    n = 100;
    Phases              = ones(Int64,n,n,n);
    Phases[:,:,20:end] .= 2;

    Cp = zeros(size(Phases));
    T  =  ones(size(Phases))*1500;
    P  =  zeros(size(Phases));

    args=(;T=T);
    compute_heatcapacity!(Cp, Mat_tup, Phases, args)    # computation routine w/out P (not used in most heat capacity formulations)     
    @test sum(Cp[1,1,:]) ≈ 121399.0486067196

    # check with array of constant properties (and no required input args)
    args1=(;);
    compute_heatcapacity!(Cp, Mat_tup1, Phases, args1)    # computation routine w/out P (not used in most heat capacity formulations)     
    @test sum(Cp[1,1,:]) ≈ 109050.0

    num_alloc = @allocated compute_heatcapacity!(Cp, Mat_tup, Phases, args)
    @test sum(Cp[1,1,:]) ≈ 121399.0486067196
    @test num_alloc <= 32

    # test if we provide phase ratios
    PhaseRatio  = zeros(n,n,n,3);
    for i in CartesianIndices(Phases)
        iz = Phases[i]
        I = CartesianIndex(i,iz)
        PhaseRatio[I] = 1.0  
    end
    compute_heatcapacity!(Cp, Mat_tup, PhaseRatio, args)
    num_alloc = @allocated compute_heatcapacity!(Cp, Mat_tup, PhaseRatio, args)
    @test sum(Cp[1,1,:]) ≈ 121399.0486067196
    @test num_alloc <= 32
=#

    # -----------------------

   

end
