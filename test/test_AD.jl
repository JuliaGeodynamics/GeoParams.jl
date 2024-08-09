using Test, GeoParams
import AbstractDifferentiation as AD, ReverseDiff, ForwardDiff

T        = 1e3
P        = 1e9
backends = AD.ReverseDiffBackend(), AD.ForwardDiffBackend()
pkg      = "ReverseDiff", "ForwardDiff"

for (name, backend) in zip(pkg, backends)
    @testset "name compatibility with: density" begin
        ρ = ConstantDensity()
        @test AD.derivative(backend, x -> compute_density(ρ, (; T=x, P=1e9)), T)[1] == 0
        @test AD.derivative(backend, x -> compute_density(ρ, (; T=1e3, P=x)), P)[1] == 0

        ρ = PT_Density()
        @test AD.derivative(backend, x -> compute_density(ρ, (; T=x, P=1e9)), T)[1] ≈ -0.087
        @test AD.derivative(backend, x -> compute_density(ρ, (; T=1e3, P=x)), P)[1] == 2.9e-6

        ρ = Compressible_Density()
        @test AD.derivative(backend, x -> compute_density(ρ, (; T=1e3, P=x)), P)[1] ≈ 7.88301730e-6

        ρ = MeltDependent_Density()
        @test AD.derivative(backend, x -> compute_density(ρ, (; ϕ =x)), 0.2)[1] == -700

        ρ = T_Density()
        @test AD.derivative(backend, x -> compute_density(ρ, (; T=x, P=1e9)), T)[1] ≈ -0.087
    end

    @testset "ForwardDiff compatibility with: heat capacity" begin
        Cp = T_HeatCapacity_Whittington()
        @test AD.derivative(backend, x -> compute_heatcapacity(Cp, (; T=x)), T)[1] ≈ 0.145639823

        Cp = Latent_HeatCapacity(Q_L=500e3)
        @test AD.derivative(backend, x -> compute_heatcapacity(Cp, (; T=x)), T)[1] == 0 
    end

    @testset "ForwardDiff compatibility with: conductivity" begin
        K = ConstantConductivity()
        @test AD.derivative(backend, x -> compute_conductivity(K, (; T=x)), T)[1] == 0

        K = T_Conductivity_Whittington()
        @test AD.derivative(backend, x -> compute_conductivity(K, (; T=x)), T)[1] ≈ -0.00019522103

        K = T_Conductivity_Whittington_parameterised()
        @test AD.derivative(backend, x -> compute_conductivity(K, (; T=x)), T)[1] ≈ -0.00064766553
        
        K = TP_Conductivity()
        @test AD.derivative(backend, x -> compute_conductivity(K, (; T=x)), T)[1] ≈ -0.00040864570
    end

    @testset "ForwardDiff compatibility with: Diffusion" begin
        import GeoParams.Diffusion
        # Define a linear viscous creep law ---------------------------------
        diffusion_law = Diffusion.dry_anorthite_Rybacki_2006
        p             = SetDiffusionCreep(diffusion_law; n = 1NoUnits)
        args          = (; T=T, P=P)
        TauII         = 1e6
        @test AD.derivative(backend, x -> compute_εII(p, TauII, (; T=x, P=P)), T)[1] ≈ 5.7562812856e-33
        @test AD.derivative(backend, x -> compute_εII(p, TauII, (; T=T, P=x)), P)[1] ≈ -2.8543543565e-40
        
        εII           = compute_εII(p, TauII, args)
        @test AD.derivative(backend, x -> compute_τII(p, εII, (; T=x, P=P)), T)[1] ≈ -58211.55812135427
        @test AD.derivative(backend, x -> compute_τII(p, εII, (; T=T, P=x)), P)[1] ≈ 0.0028865235432076496
    end

    @testset "ForwardDiff compatibility with: Dislocation" begin
        import GeoParams.Dislocation
        # Define a linear viscous creep law ---------------------------------
        diffusion_law = Dislocation.dry_olivine_Hirth_2003
        p             = SetDislocationCreep(diffusion_law; n = 1NoUnits)
        args          = (; T=T, P=P)
        TauII         = 1e6
        @test AD.derivative(backend, x -> compute_εII(p, TauII, (; T=x, P=P)), T)[1] ≈ 4.1522654949e-40
        @test AD.derivative(backend, x -> compute_εII(p, TauII, (; T=T, P=x)), P)[1] ≈ -1.0685977376e-47
        
        εII           = compute_εII(p, TauII, args)
        @test AD.derivative(backend, x -> compute_τII(p, εII, (; T=x, P=P)), T)[1] ≈ -65427.866979
        @test AD.derivative(backend, x -> compute_τII(p, εII, (; T=T, P=x)), P)[1] ≈ 0.00168380540
    end
end
