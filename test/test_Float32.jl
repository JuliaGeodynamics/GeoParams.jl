using Test
using GeoParams
import ForwardDiff
import GeoParams: precision_of, convert_precision, argument_at
import GeoParams.Dislocation, GeoParams.Diffusion

# Laws whose scalar API is `f(law, x; T, P, ...)`; the tolerance is loose enough
# for the exponentials in the creep laws to differ between precisions.
const RTOL = Dict(Float64 => 1.0e-12, Float32 => 1.0e-4)

function invariant_laws()
    return (
        SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003),
        SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006),
        LinearViscous(η = 1.0e20Pa * s),
        ArrheniusType(),
    )
end

@testset "Float32 calculations" begin

    @testset "precision_of and convert_precision" begin
        @test precision_of(1.0f0) === Float32
        @test precision_of(1.0) === Float64
        @test precision_of(1) === Float64            # integers express no preference
        @test precision_of(1.0f6 * Pa) === Float32
        @test precision_of(Float32[1, 2]) === Float32
        @test precision_of((1.0f0, 2.0f0)) === Float32

        @test convert_precision(Float32, 1.0) === 1.0f0
        @test convert_precision(Float32, 1) === 1.0f0
        @test ustrip(convert_precision(Float32, 1.0Pa)) === 1.0f0
        @test unit(convert_precision(Float32, 1.0Pa)) === unit(1.0Pa)
        @test convert_precision(Float32, (1.0, 2.0)) === (1.0f0, 2.0f0)
        # arrays are left alone: converting one would allocate a copy per call
        v = [1.0, 2.0]
        @test convert_precision(Float32, v) === v

        @test argument_at(2.5, 7) === 2.5           # a shared value, not indexed
        @test argument_at([1.0, 2.0], 2) === 2.0
    end

    @testset "creep-law scalar paths" begin
        for law in invariant_laws(), T in (Float64, Float32)
            τ = T(1.0e6)
            args = (; T = 1200.0, P = 1.0e9)  # Float64 literals
            ε = compute_εII(law, τ, args)
            @test typeof(ε) === T
            @test typeof(compute_τII(law, ε, args)) === T
            @test typeof(dεII_dτII(law, τ, args)) === T
            @test typeof(dτII_dεII(law, ε, args)) === T

            # a Float32 evaluation must track the Float64 reference
            ref = compute_εII(law, Float64(τ), args)
            @test ε ≈ ref rtol = RTOL[T]
        end
    end

    @testset "elasticity follows the stress it is given" begin
        el = ConstantElasticity(G = 10.0e9Pa, Kb = 20.0e9Pa)
        for T in (Float64, Float32)
            τ, ε = T(1.0e6), T(1.0e-15)
            @test typeof(compute_εII(el, τ; τII_old = 0.0, dt = 1.0)) === T
            @test typeof(compute_τII(el, ε; τII_old = 0.0, dt = 1.0)) === T
            @test typeof(dεII_dτII(el, τ; τII_old = 0.0, dt = 1.0)) === T
            @test typeof(compute_εvol(el, τ; P_old = 0.0, dt = 1.0)) === T
            @test typeof(compute_p(el, ε; P_old = 0.0, dt = 1.0)) === T
            @test typeof(dτII_dεII(el, zero(T); τII_old = zero(T), dt = one(T))) === T
            @test typeof(dp_dεvol(el, zero(T); P_old = zero(T), dt = one(T))) === T
        end
    end

    @testset "in-place routines keep the destination eltype" begin
        law = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
        for T in (Float64, Float32)
            τ = fill(T(1.0e6), 4)
            ε = similar(τ)
            # scalar keyword arguments are shared by every element
            compute_εII!(ε, law, τ; T = 1200.0, P = 1.0e9)
            @test eltype(ε) === T
            @test all(x -> x ≈ compute_εII(law, τ[1]; T = 1200.0, P = 1.0e9), ε)

            # ... and array keyword arguments are read per element
            ε2 = similar(τ)
            compute_εII!(ε2, law, τ, (; T = fill(1200.0, 4), P = fill(1.0e9, 4)))
            @test ε2 == ε

            # the input array need not share the destination's eltype
            ε3 = similar(τ)
            compute_εII!(ε3, law, Float64.(τ); T = 1200.0, P = 1.0e9)
            @test eltype(ε3) === T
            @test ε3 ≈ ε rtol = RTOL[T]

            τ_out = similar(τ)
            compute_τII!(τ_out, law, ε; T = 1200.0, P = 1.0e9)
            @test eltype(τ_out) === T
            @test all(x -> x ≈ τ[1], τ_out)
        end
    end

    @testset "unsupported keyword arguments are ignored, not rejected" begin
        law = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
        ε = zeros(Float32, 2)
        compute_εII!(ε, law, fill(1.0f6, 2); T = 1200.0, P = 1.0e9, τII_old = 0.0)
        @test eltype(ε) === Float32
    end

    @testset "property laws follow their solver-state keyword" begin
        for T in (Float64, Float32)
            P, Temp, ϕ = T(1.0e9), T(1200.0), T(0.3)

            @test typeof(compute_density(PT_Density(), (; P, T = Temp))) === T
            @test typeof(compute_density(Compressible_Density(), (; P))) === T
            @test typeof(compute_conductivity(T_Conductivity_Whittington(), (; T = Temp))) === T
            @test typeof(compute_meltfraction(MeltingParam_Caricchi(), (; T = Temp))) === T
            @test typeof(compute_dϕdT(MeltingParam_Caricchi(), (; T = Temp))) === T
            @test typeof(compute_permeability(PowerLawPermeability(), (; ϕ))) === T
        end
    end

    @testset "Unitful input keeps its numeric type" begin
        law = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
        ε32 = compute_εII(law, 1.0f6Pa; T = 1200.0K, P = 1.0e9Pa)
        ε64 = compute_εII(law, 1.0e6Pa; T = 1200.0K, P = 1.0e9Pa)
        @test typeof(ustrip(ε32)) === Float32
        @test typeof(ustrip(ε64)) === Float64
        @test unit(ε32) == unit(ε64)
        @test ustrip(ε32) ≈ ustrip(ε64) rtol = RTOL[Float32]
    end

    @testset "stored parameters keep their published Float64 values" begin
        law = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
        @test typeof(law.n.val) === Float64
        @test typeof(compute_εII(law, 1.0f6; T = 1200.0, P = 1.0e9)) === Float32
    end

    @testset "default Float64 calls are unchanged" begin
        # Reference values, not bit patterns: the creep law's exp/pow chain and the
        # density's `@muladd` both land a few ULP apart across libm versions and
        # architectures, so an exact comparison passes only where it was recorded.
        # This tolerance is still ~4 orders tighter than any real regression.
        law = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
        @test compute_εII(law, 1.0e6; T = 1200.0, P = 1.0e9) ≈ 1.3638005232835048e-18 rtol = RTOL[Float64]
        @test compute_density(PT_Density(), (; P = 1.0e9, T = 1200.0)) ≈ 5719.36405 rtol = RTOL[Float64]
    end

    @testset "scratch storage follows the input precision" begin
        # partial-melt viscosity: the strain-rate guard used eps(Float64)
        melt = ViscosityPartialMelt_Costa_etal_2009(η = LinearMeltViscosity())
        @test typeof(compute_εII(melt, 1.0f6; T = 1000.0f0, ϕ = 0.5f0)) === Float32

        # elastic stress rotation allocated a Float64 rotation axis
        ω = ntuple(_ -> 1.0f-2, 3)
        τ = ntuple(i -> Float32(i), 6)
        @test eltype(rotate_elastic_stress(ω, τ, 1.0f0)) === Float32

        # TAS classification work buffers
        @test computeTASclassification(Float32[50.0, 4.0]) ==
            computeTASclassification([50.0, 4.0])
    end

    @testset "iterative and wide-range laws stay in the working precision" begin
        # Herschel-Bulkley solves a Newton iteration above the yield stress; its
        # convergence test has to be reachable in the working precision
        hb = HerschelBulkley()
        ε32 = compute_εII(hb, 2.0f8; T = 1273.0f0)
        @test typeof(ε32) === Float32
        @test Float64(ε32) ≈ compute_εII(hb, 2.0e8; T = 1273.0) rtol = 1.0e-4

        # Redlich-Kwong: Float64 coefficients evaluated at a Float32 state
        @test typeof(RedlichKwong_Density()(; P = 2.0f8, T = 1200.0f0)) === Float32
    end

    @testset "no allocations on the scalar and in-place paths" begin
        # Measure through a function barrier. `@allocated` written directly in a
        # testset body measures an unoptimized top-level thunk, in which the
        # keyword splat of the `compute_*(law, x, args)` forwarding methods is
        # materialized rather than elided — which is a property of the harness,
        # not of the code under test.
        allocs(f::F, args::Vararg{Any, N}) where {F, N} = (f(args...); @allocated f(args...))

        law = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
        τ, args = 1.0f6, (; T = 1200.0, P = 1.0e9)
        τv = fill(τ, 4)
        εv = similar(τv)
        ε = compute_εII(law, τ, args)
        @test allocs(compute_εII, law, τ, args) == 0
        @test allocs(compute_τII, law, ε, args) == 0
        @test allocs(compute_εII!, εv, law, τv) == 0
        @test allocs(compute_τII!, τv, law, εv) == 0
    end

    # Every concrete law of a family is swept, so a law added later is covered
    # without touching this file.
    @testset "every law follows its argument precision" begin
        MP = GeoParams.MaterialParameters
        args32 = (;
            T = 1.0f3, P = 1.0f8, τII = 1.0f6, τII_old = 0.0f0, εII = 1.0f-15,
            εvol = 1.0f-16, ϕ = 0.1f0, dt = 1.0f0, z = 1.0f3, d = 1.0f-3, f = 1.0f0,
        )
        cases = (
            (:density, MP.Density.AbstractDensity, l -> compute_density(l, args32)),
            (:conductivity, MP.Conductivity.AbstractConductivity, l -> compute_conductivity(l, args32)),
            (:heatcapacity, MP.HeatCapacity.AbstractHeatCapacity, l -> compute_heatcapacity(l, args32)),
            (:radioactive_heat, MP.RadioactiveHeat.AbstractRadioactiveHeat, l -> compute_radioactive_heat(l, args32)),
            (:latent_heat, MP.LatentHeat.AbstractLatentHeat, l -> compute_latent_heat(l, args32)),
            (:meltfraction, GeoParams.MeltingParam.AbstractMeltingParam, l -> compute_meltfraction(l, args32)),
            (:permeability, GeoParams.AbstractPermeability, l -> compute_permeability(l, args32)),
            (:creep_εII, GeoParams.AbstractCreepLaw, l -> compute_εII(l, args32.τII, args32)),
            (:creep_τII, GeoParams.AbstractCreepLaw, l -> compute_τII(l, args32.εII, args32)),
            (:yieldfunction, GeoParams.AbstractPlasticity, l -> compute_yieldfunction(l, args32)),
        )
        # Laws the sweep cannot drive, with the reason. Not precision failures.
        unsupported = Dict(
            (:density, :Vector_Density) => "reads ρ from a user vector; needs an `index`",
            (:heatcapacity, :Vector_HeatCapacity) => "reads Cp from a user vector; needs an `index`",
            (:meltfraction, :Vector_MeltingParam) => "reads ϕ from a user vector; needs an `index`",
            (:creep_τII, :NonLinearPeierlsCreep) => "no compute_τII method",
        )
        # Scanning GeoParams' own bindings keeps the sweep to this package's laws;
        # `subtypes` would also pick up laws defined by any other loaded package.
        concrete_laws(family) = sort!(
            Type[
                t for t in (getfield(GeoParams, n) for n in names(GeoParams; all = true) if isdefined(GeoParams, n))
                    if t isa Type && !isabstracttype(t) && t <: family
            ];
            by = nameof,
        )

        for (name, family, call) in cases
            @testset "$name/$(nameof(T))" for T in concrete_laws(family)
                haskey(unsupported, (name, nameof(T))) && continue
                @test call(T()) isa Float32
            end
        end
    end

    # Phase dispatch must not widen either: neither the matched branch nor the
    # no-match fallback may introduce a Float64.
    @testset "phase dispatch" begin
        args32 = (; T = 1.0f3, P = 1.0f8)
        phases = (
            SetMaterialParams(; Phase = 1, Density = PT_Density()),
            SetMaterialParams(; Phase = 2, Density = PT_Density()),
        )
        @test compute_density(phases, 1, args32) isa Float32
        @test_throws "phase not found in MaterialParams" compute_density(phases, 99, args32)
        @test compute_density(phases, (0.4f0, 0.6f0), args32) isa Float32
        @test_throws "phase not found in MaterialParams" GeoParams.nphase(v -> compute_density(v, args32), 99, phases)
        @test GeoParams.nphase_ratio(v -> compute_density(v, args32), (0.4f0, 0.6f0), phases) isa Float32

        # every family's phase-dispatch method agrees with a direct evaluation
        CR = GeoParams.MaterialParameters.ConstitutiveRelationships
        a = (; T = 1.0f3, P = 1.0f8, τII = 1.0f6, ϕ = 0.1f0)
        mp = SetMaterialParams(;
            Phase = 1,
            Conductivity = T_Conductivity_Whittington(),
            HeatCapacity = T_HeatCapacity_Whittington(),
            LatentHeat = ConstantLatentHeat(),
            RadioactiveHeat = ExpDepthDependentRadioactiveHeat(),
            Melting = MeltingParam_Caricchi(),
            Plasticity = DruckerPrager(Ψ = 10),
            SeismicVelocity = ConstantSeismicVelocity(),
        )
        τ3 = (1.0f6, 2.0f6, 3.0f6)
        for (f, x) in (
                (compute_conductivity, (a,)),
                (compute_heatcapacity, (a,)),
                (compute_latent_heat, (a,)),
                (compute_radioactive_heat, ((; z = 1.0f3),)),
                (compute_meltfraction, (a,)),
                (compute_dϕdT, (a,)),
                (compute_yieldfunction, (a,)),
                (CR.∂Q∂τ, (τ3,)),
                (CR.∂Q∂τII, (1.0f6,)),
                (CR.∂Q∂P, (1.0f8,)),
                (CR.plastic_strain_rate, (τ3, 1.0f-15)),
            )
            @test f((mp,), 1, x...) === f(mp, x...)
            @test eltype(f(mp, x...)) === Float32
        end
        @test compute_wave_velocity((mp,), 1, (; wave = :Vp)) == compute_wave_velocity(mp, (; wave = :Vp))
    end

    @testset "plastic flow direction and multiplier" begin
        CR = GeoParams.MaterialParameters.ConstitutiveRelationships
        τ3 = (1.0f0, 2.0f0, 3.0f0)
        τ6 = (1.0f0, 2.0f0, 3.0f0, 4.0f0, 5.0f0, 6.0f0)
        for p in (DruckerPrager(), DruckerPrager_regularised())
            @test CR.∂Q∂τxx(p, τ3) isa Float32
            @test CR.∂Q∂τxy(p, τ3) isa Float32
            @test CR.∂Q∂τzz(p, τ6) isa Float32
            @test CR.∂Q∂τyz(p, τ6) isa Float32
        end
        @test second_invariant(τ3) isa Float32
        @test second_invariant(τ6) isa Float32
        A4 = ntuple(i -> Float32(i), 4)
        @test second_invariant_staggered(A4, A4, 1.0f0) isa Float32
        @test second_invariant_staggered(A4, A4, A4, (1.0f0, 2.0f0, 3.0f0)) isa Float32
        pc = DruckerPragerCap()
        for g in (CR.∂Q∂τxx, CR.∂Q∂τyy, CR.∂Q∂τxy)
            @test g(pc, τ3; P = 1.0f6) isa Float32
        end
        # shear (Drucker-Prager) branch of the flow potential, Q = τII - sinΨ⋅P
        pc10 = DruckerPragerCap(Ψ = 10)
        @test CR.∂Q∂τII(pc10, 1.0f6; P = 1.0f8) === 0.5f0
        @test CR.∂Q∂P(pc10, 1.0f8; τII = 1.0f6) ≈ -sind(10.0f0)
        mp = SetMaterialParams(; Phase = 1, Plasticity = DruckerPrager())
        @test CR.plastic_strain(mp, τ3, 1.0f-15) isa Float32
        @test CR.plastic_strain((mp,), 1, τ3, 1.0f-15) isa Float32
        p = DruckerPrager()
        λ32 = CR.lambda(1.0f6, p, 1.0f20, 1.0f19; K = 1.0f10, dt = 1.0f10, h = 1.0f5, τij = τ3)
        @test λ32 isa Float32
        @test λ32 ≈ CR.lambda(1.0e6, p, 1.0e20, 1.0e19; K = 1.0e10, dt = 1.0e10, h = 1.0e5, τij = (1.0, 2.0, 3.0)) rtol = RTOL[Float32]
    end

    @testset "arg-independent laws follow the precision of their arguments" begin
        CR = GeoParams.MaterialParameters.ConstitutiveRelationships
        args = (; T = 1.0f3, P = 1.0f8, dt = 1.0f10)
        @test compute_viscosity(LinearViscous(), args) isa Float32
        @test compute_viscosity(ConstantElasticity(), args) isa Float32
        @test compute_viscosity_εII(LinearViscous(), 1.0f-15, args) isa Float32
        @test compute_viscosity_τII(LinearViscous(), 1.0f6, args) isa Float32
        @test compute_viscosity_εII(ConstantElasticity(), 1.0f-15, args) isa Float32
        @test compute_elastoviscosity(ConstantElasticity(), 1.0f20, 1.0f10) isa Float32
        for p in (DruckerPrager(), DruckerPrager_regularised())
            @test CR.∂Q∂P(p, 1.0f8) isa Float32
            @test CR.∂F∂P(p, 1.0f8) isa Float32
        end
    end

    @testset "melt, solubility and Herschel-Bulkley laws" begin
        @test compute_dϕdT(SmoothMelting(); T = 1.0f3) isa Float32
        cs = ViscosityPartialMelt_Costa_etal_2009()
        @test dτII_dεII(cs, 1.0f-15; ϕ = 0.1f0, T = 1.0f3) isa Float32
        @test dεII_dτII(cs, 1.0f6; ϕ = 0.1f0, T = 1.0f3) isa Float32
        gm = GiordanoMeltViscosity()
        @test dεII_dτII(gm, 1.0f6; T = 1.0f3) isa Float32
        for s in (Liu2005_Solubility(), Mafic_Solubility())
            @test all(x -> x isa Float32, compute_dissolved(s; P = 1.0f8, T = 1.1f3, X_co2 = 0.5f0))
            @test all(x -> x isa Float32, compute_dissolved(s; P = 1.0f8, T = 1.1f3))
        end
        # laboratory-scale viscosities: the Newton residual must stay inside Float32 range
        hb = HerschelBulkley()
        for τ in (1.0f0, 1.0f6, 1.0f8, 1.0f9)
            ε32 = compute_εII(hb, τ; T = 1.0f3)
            @test ε32 isa Float32
            @test ε32 ≈ compute_εII(hb, Float64(τ); T = 1.0e3) rtol = 1.0e-4
        end
        # derivatives stay finite and accurate where 1/η² underflows Float32
        for τ in (5.0e7, 2.0e8, 1.0e9)
            d32 = ForwardDiff.derivative(t -> compute_εII(hb, t; T = 1.0f3), Float32(τ))
            d64 = ForwardDiff.derivative(t -> compute_εII(hb, t; T = 1.0e3), τ)
            @test d32 isa Float32
            @test d32 ≈ d64 rtol = 1.0e-3
        end
        @test precision_of((; T = ForwardDiff.Dual(1.0f3, 1.0f0))) === Float32

        τ = Float32[1.0f6, 2.0f8]
        ε_ref = [compute_εII(hb, t; T = 1.0f3) for t in τ]
        ε = similar(τ)
        compute_εII!(ε, hb, τ; T = 1.0f3)
        @test ε == ε_ref
        compute_εII!(ε, hb, τ; T = fill(1.0f3, 2))
        @test ε == ε_ref
        η = compute_viscosity_τII(hb, 2.0f8, (; T = 1.0f3))
        @test η isa Float32
        @test η ≈ 2.0f8 / (2 * ε_ref[2]) rtol = RTOL[Float32]
        @test_throws "compute_hb_εII: iterations did not converge" compute_εII(hb, NaN32; T = 1.0f3)
    end

    # Float32 intermediates of the dislocation-creep derivatives leave the
    # Float32 range for answers that are inside it; the result is recomputed
    # in Float64 and narrowed back.
    @testset "wider recomputation on overflow and underflow" begin
        f(x) = x * x / x
        @test GeoParams.retry_wider(f, f(1.0f30), 1.0f30) === 1.0f30      # overflow
        @test GeoParams.retry_wider(f, f(1.0f-30), 1.0f-30) === 1.0f-30    # underflow
        @test GeoParams.retry_wider(f, f(Float16(300)), Float16(300)) === Float16(300)
        @test GeoParams.retry_wider(f, Inf, 1.0e300) === Inf               # Float64 is not widened

        law = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
        d32 = dτII_dεII(law, 1.0f-20; T = 400.0f0, P = 1.0f9)
        @test d32 isa Float32
        @test isfinite(d32)
        @test d32 ≈ dτII_dεII(law, 1.0e-20; T = 400.0, P = 1.0e9) rtol = RTOL[Float32]
        # the Float64 value, 7e-51, lies below the Float32 range
        @test dεII_dτII(law, 1.0f0; T = 800.0f0, P = 1.0f9) ===
            Float32(dεII_dτII(law, 1.0; T = 800.0, P = 1.0e9))
    end

    @testset "unitful input to property and creep laws" begin
        lv = LinearViscous(η = 1.0e20Pa * s)
        @test dεII_dτII(lv, 1.0e6Pa) == 5.0e-21 / (Pa * s)
        @test compute_τII(lv, 1.0e-15 / s) ≈ 2.0e5Pa
        @test dτII_dεII(lv, 1.0e-15 / s) == 2.0e20Pa * s

        cs = ViscosityPartialMelt_Costa_etal_2009()
        @test ustrip(compute_τII(cs, 1.0e-15 / s; ϕ = 0.1, T = 1000.0K)) ≈
            compute_τII(cs, 1.0e-15; ϕ = 0.1, T = 1000.0)
        @test ustrip(compute_heatcapacity(T_HeatCapacity_Whittington(); T = 1000.0K)) ≈
            compute_heatcapacity(T_HeatCapacity_Whittington(); T = 1000.0)
        # a porosity given in percent keeps its unit
        @test unit(PowerLawPermeability()(; ϕ = 10.0u"percent")) == u"m^2 * percent^3"
        @test unit(CarmanKozenyPermeability()(; ϕ = 10.0u"percent")) == u"m^2 * percent^3"

        @test precision_of((; T = 1.0f3K)) === Float32
    end

    @testset "typed unpacking of array-valued parameters" begin
        nt64 = (; x = GeoUnit([1.0, 2.0]))
        GeoParams.@unpack_val Float32 x = nt64
        @test x == Float32[1, 2]
        @test eltype(x) === Float32

        nt32 = (; x = GeoUnit(Float32[1, 2]))
        GeoParams.@unpack_val Float32 x = nt32
        @test x === nt32.x.val                  # already Float32: not copied
    end
end
