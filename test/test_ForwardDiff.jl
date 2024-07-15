using Test, GeoParams
import ForwardDiff as AD

T = 1e3
P = 1e9

@testset "ForwardDiff compatibility with: density" begin
    ρ = ConstantDensity()
    @test AD.derivative(x -> compute_density(ρ, (; T=x, P=1e9)), T) == 0
    @test AD.derivative(x -> compute_density(ρ, (; T=1e3, P=x)), P) == 0

    ρ = PT_Density()
    @test AD.derivative(x -> compute_density(ρ, (; T=x, P=1e9)), T) == -0.08700000000000001
    @test AD.derivative(x -> compute_density(ρ, (; T=1e3, P=x)), P) == 2.9e-6

    ρ = Compressible_Density()
    @test AD.derivative(x -> compute_density(ρ, (; T=1e3, P=x)), P) == 7.883017302531231e-6

    ρ = MeltDependent_Density()
    @test AD.derivative(x -> compute_density(ρ, (; ϕ =x)), 0.2) == -700

    ρ = T_Density()
    @test AD.derivative(x -> compute_density(ρ, (; T=x, P=1e9)), T) == -0.08700000000000001
end

@testset "ForwardDiff compatibility with: heat capacity" begin
    Cp = T_HeatCapacity_Whittington()
    @test AD.derivative(x -> compute_heatcapacity(Cp, (; T=x)), T) == 0.145639823248696
    compute_heatcapacity(Cp, args)

    Cp = Latent_HeatCapacity(Q_L=500e3)
    @test AD.derivative(x -> compute_heatcapacity(Cp, (; T=x)), T) == 0 
end

@testset "ForwardDiff compatibility with: conductivity" begin
    K = ConstantConductivity()
    @test AD.derivative(x -> compute_conductivity(K, (; T=x)), T) == 0

    K = T_Conductivity_Whittington()
    @test AD.derivative(x -> compute_conductivity(K, (; T=x)), T) == -0.0001952210298486253

    K = T_Conductivity_Whittington_parameterised()
    @test AD.derivative(x -> compute_conductivity(K, (; T=x)), T) == -0.0006476655349999987
    
    K = TP_Conductivity()
    @test AD.derivative(x -> compute_conductivity(K, (; T=x)), T) == -0.000408645701590356
end

@testset "ForwardDiff compatibility with: Diffusion" begin
    import GeoParams.Diffusion
    # Define a linear viscous creep law ---------------------------------
    diffusion_law = Diffusion.dry_anorthite_Rybacki_2006
    p             = SetDiffusionCreep(diffusion_law; n = 1NoUnits)
    args          = (; T=T, P=P)
    TauII         = 1e6
    @test AD.derivative(x -> compute_εII(p, TauII, (; T=x, P=P)), T) == 5.756281285647566e-33
    @test AD.derivative(x -> compute_εII(p, TauII, (; T=T, P=x)), P) == -2.854354356519454e-40
    
    εII           = compute_εII(p, TauII, args)
    @test AD.derivative(x -> compute_τII(p, εII, (; T=x, P=P)), T) == -58211.55812135427
    @test AD.derivative(x -> compute_τII(p, εII, (; T=T, P=x)), P) == 0.0028865235432076496
end

@testset "ForwardDiff compatibility with: Dislocation" begin
    import GeoParams.Dislocation
    # Define a linear viscous creep law ---------------------------------
    diffusion_law = Dislocation.dry_olivine_Hirth_2003
    p             = SetDislocationCreep(diffusion_law; n = 1NoUnits)
    args          = (; T=T, P=P)
    TauII         = 1e6
    @test AD.derivative(x -> compute_εII(p, TauII, (; T=x, P=P)), T) == 4.152265494944879e-40
    @test AD.derivative(x -> compute_εII(p, TauII, (; T=T, P=x)), P) == -1.068597737669638e-47
    
    εII           = compute_εII(p, TauII, args)
    @test AD.derivative(x -> compute_τII(p, εII, (; T=x, P=P)), T) == -65427.866979373386
    @test AD.derivative(x -> compute_τII(p, εII, (; T=T, P=x)), P) == 0.0016838054002044623
end