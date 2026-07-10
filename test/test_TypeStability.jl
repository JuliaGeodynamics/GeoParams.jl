# Type-stability tests for the scalar compute kernels.
#
# These kernels are contracted to be allocation-free and compatible with
# ForwardDiff/ReverseDiff; both break silently when a call infers to a
# non-concrete (`Union`/abstract) return type. `@inferred` asserts the return
# type is concrete, while JET's `@test_opt` catches runtime dispatch anywhere in
# the call tree. JET is scoped to `GeoParams` so instabilities originating in
# dependencies (Roots, Unitful, ...) do not mask our own.
#
# Everything here exercises the nondimensional `Float64` numeric path.

using Test, GeoParams
import GeoParams.Diffusion
import GeoParams.Dislocation
import GeoParams.Peierls
import GeoParams.NonLinearPeierls
import GeoParams.GBS

# JET tracks compiler internals and does not precompile on every Julia release.
# The `@inferred` layer runs everywhere; the JET `@test_opt` layer runs only
# where JET actually loads.
const JET_LOADABLE = try
    # Suppress the precompile-failure backtrace JET emits where it is
    # unsupported; the `@info` below records the fallback instead.
    redirect_stderr(devnull) do
        @eval import JET
    end
    true
catch
    @info "JET unavailable on this Julia version; type-stability tests run @inferred only."
    false
end

# Assert the return type is concretely inferred, and — where JET loads — that no
# runtime dispatch occurs inside GeoParams itself (dependency instabilities are
# out of our control, so they are filtered out via `target_modules`).
macro test_stable(ex)
    jet = JET_LOADABLE ? :(JET.@test_opt target_modules = (GeoParams,) $ex) : nothing
    return esc(
        quote
            @test (@inferred $ex; true)
            $jet
        end
    )
end

@testset "Type stability" begin

    @testset "Density" begin
        @test_stable compute_density(ConstantDensity(), (; T = 1.0e3, P = 1.0e9))
        @test_stable compute_density(PT_Density(), (; T = 1.0e3, P = 1.0e9))
        @test_stable compute_density(Compressible_Density(), (; T = 1.0e3, P = 1.0e9))
        @test_stable compute_density(MeltDependent_Density(), (; ϕ = 0.2))
        @test_stable compute_density(T_Density(), (; T = 1.0e3, P = 1.0e9))
    end

    @testset "Heat capacity" begin
        @test_stable compute_heatcapacity(T_HeatCapacity_Whittington(), (; T = 1.0e3))
        @test_stable compute_heatcapacity(Latent_HeatCapacity(Q_L = 500.0e3), (; T = 1.0e3))
    end

    @testset "Conductivity" begin
        @test_stable compute_conductivity(ConstantConductivity(), (; T = 1.0e3))
        @test_stable compute_conductivity(T_Conductivity_Whittington(), (; T = 1.0e3))
        @test_stable compute_conductivity(T_Conductivity_Whittington_parameterised(), (; T = 1.0e3))
        @test_stable compute_conductivity(TP_Conductivity(), (; T = 1.0e3))
    end

    @testset "Diffusion creep" begin
        p = SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006; n = 1NoUnits)
        args = (; T = 1.0e3, P = 1.0e9)
        τII = 1.0e6
        @test_stable compute_εII(p, τII, args)
        εII = compute_εII(p, τII, args)
        @test_stable compute_τII(p, εII, args)
    end

    @testset "Dislocation creep" begin
        p = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003; n = 1NoUnits)
        args = (; T = 1.0e3, P = 1.0e9)
        τII = 1.0e6
        @test_stable compute_εII(p, τII, args)
        εII = compute_εII(p, τII, args)
        @test_stable compute_τII(p, εII, args)
    end

    @testset "Peierls creep" begin
        p = SetPeierlsCreep(Peierls.dry_olivine_Goetze_1979)
        args = (; T = 873.15)
        τII = 1.0e9
        @test_stable compute_εII(p, τII, args)
        εII = compute_εII(p, τII, args)
        @test_stable compute_τII(p, εII, args)
    end

    @testset "Non-linear Peierls creep" begin
        # NonLinearPeierlsCreep implements only the stress -> strain direction.
        p = SetNonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010)
        @test_stable compute_εII(p, 1.0e9, (; T = 473.15))
    end

    @testset "Grain boundary sliding" begin
        p = SetGrainBoundarySliding(GBS.cold_dry_olivine_Hirth_2003)
        args = (; T = 923.15, d = 1.0e-3)
        τII = 1.0e6
        @test_stable compute_εII(p, τII, args)
        εII = compute_εII(p, τII, args)
        @test_stable compute_τII(p, εII, args)
    end

    @testset "Elasticity" begin
        a = ConstantElasticity()
        argsτ = (; τII_old = 15.0e6, dt = 1.0e6)
        argsp = (; P_old = 5.0e6, dt = 1.0e6)
        @test_stable compute_εII(a, 20.0e6, argsτ)
        @test_stable compute_τII(a, 1.0e-15, argsτ)
        @test_stable compute_εvol(a, 10.0e6, argsp)
    end

    @testset "Viscosity (linear rheology)" begin
        creep = LinearViscous()
        args = (; P = 1.0e9, T = 1.0e3, dt = 1.0e8)
        εII, τII = 1.0e-15, 1.0e6
        @test_stable compute_viscosity_εII(creep, εII, args)
        @test_stable compute_viscosity_τII(creep, τII, args)
        @test_stable compute_viscosity(creep, args)
    end

    @testset "CompositeRheology" begin
        c1 = CompositeRheology(
            SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006),
            SetDislocationCreep(Dislocation.dry_anorthite_Rybacki_2006),
            ConstantElasticity(),
        )
        c2 = CompositeRheology(LinearViscous(), ConstantElasticity())
        args = (; T = 900.0, P = 1.0e9, d = 100.0e-6, τII_old = 1.0e6, dt = 1.0e8)
        εII, τII = 1.0e-12, 2.0e6
        @test_stable compute_τII(c1, εII, args)
        @test_stable compute_εII(c1, τII, args)
        @test_stable compute_τII(c2, εII, args)
        @test_stable compute_εII(c2, τII, args)
    end

    @testset "Plasticity" begin
        p = DruckerPrager()
        @test_stable compute_yieldfunction(p, (; τII = 20.0e6, P = 1.0e6))
    end

    @testset "Melt fraction" begin
        args = (; T = 1000.0)
        for melt in (
                MeltingParam_Caricchi(),
                MeltingParam_Smooth3rdOrder(),
                MeltingParam_5thOrder(),
                MeltingParam_4thOrder(),
                MeltingParam_Quadratic(),
                MeltingParam_Assimilation(),
                SmoothMelting(),
            )
            @test_stable compute_meltfraction(melt, args)
        end
    end

    @testset "Energy" begin
        @test_stable compute_latent_heat(ConstantLatentHeat(Q_L = 400.0e3), (;))
        @test_stable compute_radioactive_heat(ConstantRadioactiveHeat(), (;))
        Χ = ConstantShearheating(1.0)
        τ = (1.0, 2.0, 2.0, 3.0) .* 1.0e6
        ε = (1.0, 0.1, 0.1, 1.0)
        ε_el = (0.01, 0.01, 0.01, 0.01)
        @test_stable compute_shearheating(Χ, τ, ε)
        @test_stable compute_shearheating(Χ, τ, ε, ε_el)
    end

    @testset "Gravity" begin
        @test_stable compute_gravity(ConstantGravity())
        @test_stable compute_gravity(DippingGravity(90, 0, 9.81))
    end

    @testset "Permeability" begin
        args = (; ϕ = 0.4)
        @test_stable compute_permeability(ConstantPermeability(), args)
        @test_stable compute_permeability(HazenPermeability(), args)
        @test_stable compute_permeability(PowerLawPermeability(), args)
    end

    @testset "Stress computations" begin
        c = CompositeRheology(LinearViscous(; η = 10.0), ConstantElasticity(; G = 1.0))
        ε = (1.0, -1.1, 2.3)
        τ_o = (0.0, 0.0, 0.0)
        args = (; dt = 0.1, τII_old = 0.0)
        @test_stable compute_τij(c, ε, args, τ_o)
        @test_stable compute_p_τij(c, ε, 5.0e6, args, τ_o)
    end

    @testset "Chemical diffusion" begin
        p = SetChemicalDiffusion(Rutile.Rt_Hf_Cherniak2007_para_c)
        @test_stable compute_D(p; T = 1273.15, P = 0.0)
    end

end
