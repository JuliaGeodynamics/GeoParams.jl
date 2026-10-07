using Test, GeoParams
import DifferentiationInterface as DI
using DifferentiationInterface: AutoReverseDiff, AutoForwardDiff
import ReverseDiff, ForwardDiff

T = 1.0e3
P = 1.0e9
backends = AutoReverseDiff(), AutoForwardDiff()
pkg = "ReverseDiff", "ForwardDiff"

for (name, backend) in zip(pkg, backends)
    @testset "$name compatibility with: density" begin
        ρ = ConstantDensity()
        @test DI.derivative(x -> compute_density(ρ, (; T = x, P = 1.0e9)), backend, T) == 0
        @test DI.derivative(x -> compute_density(ρ, (; T = 1.0e3, P = x)), backend, P) == 0

        ρ = PT_Density()
        @test DI.derivative(x -> compute_density(ρ, (; T = x, P = 1.0e9)), backend, T) ≈ -0.087
        @test DI.derivative(x -> compute_density(ρ, (; T = 1.0e3, P = x)), backend, P) == 2.9e-6

        ρ = Compressible_Density()
        @test DI.derivative(x -> compute_density(ρ, (; T = 1.0e3, P = x)), backend, P) ≈ 7.8830173e-6

        ρ = MeltDependent_Density()
        @test DI.derivative(x -> compute_density(ρ, (; ϕ = x)), backend, 0.2) == -700

        ρ = T_Density()
        @test DI.derivative(x -> compute_density(ρ, (; T = x, P = 1.0e9)), backend, T) ≈ -0.087
    end

    @testset "$name compatibility with: heat capacity" begin
        Cp = T_HeatCapacity_Whittington()
        @test DI.derivative(x -> compute_heatcapacity(Cp, (; T = x)), backend, T) ≈ 0.145639823

        Cp = Latent_HeatCapacity(Q_L = 500.0e3)
        @test DI.derivative(x -> compute_heatcapacity(Cp, (; T = x)), backend, T) == 0
    end

    @testset "$name compatibility with: conductivity" begin
        k = ConstantConductivity()
        @test DI.derivative(x -> compute_conductivity(k, (; T = x)), backend, T) == 0

        # Regression: constructor must not attempt unit-dimension promotion across a,b,c,d.
        cname = Tuple("UpperCrust")
        K_nt = TP_Conductivity(cname, 0.64Watt / K / m, 807Watt / m, 77K, 0 / MPa)
        @test K_nt isa TP_Conductivity
        @test K_nt.Name == cname
        @test (typeof(K_nt.a.val) == typeof(K_nt.b.val)) && (typeof(K_nt.b.val) == typeof(K_nt.c.val)) && (typeof(K_nt.c.val) == typeof(K_nt.d.val))

        # This goes through the Parameters keyword constructor, which dispatches to String constructor.
        K_kw = TP_Conductivity(; Name = "UpperCrust", a = 0.64Watt / K / m, b = 807Watt / m, c = 77K, d = 0 / MPa)
        @test K_kw isa TP_Conductivity

        k = T_Conductivity_Whittington()
        @test DI.derivative(x -> compute_conductivity(k, (; T = x)), backend, T) ≈ -0.00019522103

        k = T_Conductivity_Whittington_parameterised()
        @test DI.derivative(x -> compute_conductivity(k, (; T = x)), backend, T) ≈ -0.00064766553

        k = TP_Conductivity()
        @test DI.derivative(x -> compute_conductivity(k, (; T = x)), backend, T) ≈ -0.0004086457
    end

    @testset "$name compatibility with: Diffusion" begin
        import GeoParams.Diffusion
        # Define a linear viscous creep law ---------------------------------
        diffusion_law = Diffusion.dry_anorthite_Rybacki_2006
        p = SetDiffusionCreep(diffusion_law; n = 1NoUnits)
        args = (; T = T, P = P)
        TauII = 1.0e6
        @test DI.derivative(x -> compute_εII(p, TauII, (; T = x, P = P)), backend, T) ≈ 4.6050250285180517e-26
        @test DI.derivative(x -> compute_εII(p, TauII, (; T = T, P = x)), backend, P) ≈ -2.283483485215563e-33

        εII = compute_εII(p, TauII, args)
        @test DI.derivative(x -> compute_τII(p, εII, (; T = x, P = P)), backend, T) ≈ -58211.55812135427
        @test DI.derivative(x -> compute_τII(p, εII, (; T = T, P = x)), backend, P) ≈ 0.0028865235432076496
    end

    @testset "$name compatibility with: Dislocation" begin
        import GeoParams.Dislocation
        # Define a linear viscous creep law ---------------------------------
        diffusion_law = Dislocation.dry_olivine_Hirth_2003
        p = SetDislocationCreep(diffusion_law; n = 1NoUnits)
        args = (; T = T, P = P)
        TauII = 1.0e6
        @test DI.derivative(x -> compute_εII(p, TauII, (; T = x, P = P)), backend, T) ≈ 4.1522654949e-40
        @test DI.derivative(x -> compute_εII(p, TauII, (; T = T, P = x)), backend, P) ≈ -1.0685977376e-47

        εII = compute_εII(p, TauII, args)
        @test DI.derivative(x -> compute_τII(p, εII, (; T = x, P = P)), backend, T) ≈ -65427.866979
        @test DI.derivative(x -> compute_τII(p, εII, (; T = T, P = x)), backend, P) ≈ 0.0016838054
    end

    @testset "$name compatibility with: CompositeRheology" begin
        # Define a range of rheological components
        v1 = SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006)
        v2 = SetDislocationCreep(Dislocation.dry_anorthite_Rybacki_2006)
        v3 = LinearViscous()
        el = ConstantElasticity()
        # composite rheology
        c1 = CompositeRheology(v1, v2, el)
        c2 = CompositeRheology(v3, el)       # composite rheology
        # arguments
        args = (T = 900.0, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)
        εII, τII = 1.0e-12, 2.0e6
        # test non-linear rheology
        @test DI.derivative(x -> compute_τII(c1, εII, (; T = x, P = P, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)), backend, T) ≈ -0.4417520137461898
        @test DI.derivative(x -> compute_τII(c1, εII, (; T = T, P = x, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)), backend, P) ≈ 1.594299818642621e-8
        @test DI.derivative(x -> compute_εII(c1, τII, (; T = x, P = P, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)), backend, T) ≈ 2.6560571677070596e-22
        @test DI.derivative(x -> compute_εII(c1, τII, (; T = T, P = x, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)), backend, P) ≈ -9.587013268449008e-30
        # test linear rheology
        @test iszero(DI.derivative(x -> compute_τII(c2, εII, (; T = x, P = P, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)), backend, T))
        @test iszero(DI.derivative(x -> compute_τII(c2, εII, (; T = T, P = x, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)), backend, P))
        @test iszero(DI.derivative(x -> compute_εII(c2, τII, (; T = x, P = P, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)), backend, T))
        @test iszero(DI.derivative(x -> compute_εII(c2, τII, (; T = T, P = x, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)), backend, P))
    end

    @testset "$name compatibility with: ChemicalDiffusion" begin

        Hf_Rt_para = Rutile.Rt_Hf_Cherniak2007_para_c
        Hf_Rt_para = SetChemicalDiffusion(Hf_Rt_para)
        @test DI.derivative(x -> compute_D(Hf_Rt_para, T = x, P = 0), backend, T) ≈ 2.7517698e-25 atol = 1.0e-28
        @test DI.derivative(x -> compute_D(Hf_Rt_para, T = 1273.15, P = x), backend, P) ≈ 0.0

        # test only if backend is ForwardDiff because it is not compatible with ReverseDiff
        if name == "ForwardDiff"
            melt_major = Melt.Melt_multicomponent_major_Guo2020_SiO2_basaltic
            melt_major = SetMulticompChemicalDiffusion(melt_major)
            @test DI.derivative(x -> compute_D(melt_major, T = x), backend, T)[1] ≈ -1.34176e-16 atol = 1.0e-20
        end

    end

end
