eousing Test
using GeoParams
@testset "SeismicVelocity.jl" begin
    # This tests the MaterialParameters structure
    CharUnits_GEO   =   GEO_units(viscosity=1e19, length=10km);
            
    # Constant seismic velocity capacity
    x = ConstantSeismicVelocity() 
    @test isbits(x)==true
    info     =  param_info(x)

    x_nd    = x
    x_nd     = nondimensionalize(x_nd,CharUnits_GEO)

    @test Value(x.Vp) ≈ 8.1km/s
    @test Value(x.Vs) ≈ 4.5km/s
    @test Value(x_nd.Vp) ≈ 8.1e11
    @test Value(x_nd.Vs) ≈ 4.5e11
    
    @test UnitValue(compute_pwave_velocity(x_nd, random_name=1)) ≈ 8.1e11


    # Check that it works if we give a phase array
    MatParam    =   Array{MaterialParams, 1}(undef, 2);
    MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                        SeismicVelocity  = ConstantSeismicVelocity());

    MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                        SeismicVelocity  = ConstantSeismicVelocity(Vp=10km/s));

    Mat_tup = Tuple(MatParam)

  
    # test computing material properties
    n = 100;
    Phases              = ones(Int64,n,n,n);
    Phases[:,:,20:end] .= 2;

    Vp = zeros(size(Phases));
    T  =  ones(size(Phases))*1500;
    P  =  zeros(size(Phases));

    args=(;T=T);
    compute_pwave_velocity!(Vp, Mat_tup, Phases, args)





    Vp_cor,Vs_cor = melt_correction(26.0,94.5,61.0,2802.0,3198.0,7.4,4.36,0.01,0.15) 
    @test  [Vp_cor,Vs_cor] ≈ [7.336238790906285, 4.314027804335563];

    Vs_anel = anelastic_correction(0,4.36734,5.0,1250.0)
    @test  Vs_anel ≈ 4.1182815519599325;

end
