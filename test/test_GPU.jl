# Runs only with `--backend=CUDA|AMDGPU|Metal` (see runtests.jl). Every case evaluates a
# nondimensionalized Float32 law inside a one-thread kernel and compares the result with
# the same Float32 call on the host.
using Test, GeoParams
import Adapt, ForwardDiff

const backend = get(ENV, "GEOPARAMS_TEST_BACKEND", "CPU")
if backend == "CUDA"
    using CUDA
    const ArrayT = CuArray
    @eval launch(k, args...) = (@cuda threads = 1 k(args...); CUDA.synchronize())
elseif backend == "AMDGPU"
    using AMDGPU
    const ArrayT = ROCArray
    @eval launch(k, args...) = (@roc groupsize = 1 k(args...); AMDGPU.synchronize())
elseif backend == "Metal"
    using Metal
    const ArrayT = MtlArray
    @eval launch(k, args...) = (@metal threads = 1 k(args...); Metal.synchronize())
else
    error("test_GPU.jl needs GEOPARAMS_TEST_BACKEND = CUDA, AMDGPU or Metal, got $backend")
end

const CR = GeoParams.MaterialParameters.ConstitutiveRelationships
const CD = GEO_units()
nd(x) = nondimensionalize(x, CD)
f32(x) = convert_precision(Float32, x)

store!(out, r::Number) = (out[1] = r; nothing)
function store!(out, r::Tuple)
    for i in eachindex(r)
        out[i] = r[i]
    end
    return nothing
end
function kernel!(out, f::F, args) where {F}
    store!(out, f(args...))
    return nothing
end
flat(r::Number) = [Float64(r)]
flat(r::Tuple) = Float64[r...]

# `args` are host values; arrays inside them are moved to the device for the kernel
function test_gpu(f, args; rtol = 1.0e-4)
    ref = flat(f(args...))
    out = ArrayT(zeros(Float32, 8))
    launch(kernel!, out, f, Adapt.adapt(ArrayT, args))
    for (g, r) in zip(Array(out), ref)
        @test g ≈ r rtol = rtol atol = 1.0e-30
    end
    return
end

# same physical state for every law, nondimensionalized
const a = map(
    Float32, (;
        T = nd(1000.0K), P = nd(1.0e8Pa), τII = nd(1.0e6Pa), τII_old = nd(0.0Pa), εII = nd(1.0e-15 / s),
        εvol = nd(1.0e-16 / s), ϕ = 0.1, dt = nd(1.0e10s), z = nd(1.0e3m), d = nd(1.0e-3m), f = nd(1.0Pa),
    )
)

c_density(l, a) = compute_density(l, a)
c_conductivity(l, a) = compute_conductivity(l, a)
c_heatcapacity(l, a) = compute_heatcapacity(l, a)
c_radioactive(l, a) = compute_radioactive_heat(l, a)
c_latent(l, a) = compute_latent_heat(l, a)
c_meltfraction(l, a) = compute_meltfraction(l, a)
c_dϕdT(l, a) = compute_dϕdT(l, a)
c_permeability(l, a) = compute_permeability(l, a)
c_εII(l, a) = compute_εII(l, a.τII, a)
c_τII(l, a) = compute_τII(l, a.εII, a)
c_yield(l, a) = compute_yieldfunction(l, a)
c_visc_εII(l, a) = compute_viscosity_εII(l, a.εII, a)
c_visc_τII(l, a) = compute_viscosity_τII(l, a.τII, a)

@testset "GPU ($backend)" begin
    if backend == "Metal"
        @test Base.get_extension(GeoParams, :GeoParamsMetalExt) !== nothing
    end

    @testset "every law" begin
        MP = GeoParams.MaterialParameters
        cases = (
            (:density, MP.Density.AbstractDensity, c_density),
            (:conductivity, MP.Conductivity.AbstractConductivity, c_conductivity),
            (:heatcapacity, MP.HeatCapacity.AbstractHeatCapacity, c_heatcapacity),
            (:radioactive_heat, MP.RadioactiveHeat.AbstractRadioactiveHeat, c_radioactive),
            (:latent_heat, MP.LatentHeat.AbstractLatentHeat, c_latent),
            (:meltfraction, GeoParams.MeltingParam.AbstractMeltingParam, c_meltfraction),
            (:dϕdT, GeoParams.MeltingParam.AbstractMeltingParam, c_dϕdT),
            (:permeability, GeoParams.AbstractPermeability, c_permeability),
            (:creep_εII, GeoParams.AbstractCreepLaw, c_εII),
            (:creep_τII, GeoParams.AbstractCreepLaw, c_τII),
            (:yieldfunction, GeoParams.AbstractPlasticity, c_yield),
        )
        # read from an `index` into a user vector; tested below
        vector_laws = (:Vector_Density, :Vector_HeatCapacity, :Vector_MeltingParam)
        concrete_laws(family) = sort!(
            Type[
                t for t in (getfield(GeoParams, n) for n in names(GeoParams; all = true) if isdefined(GeoParams, n))
                    if t isa Type && !isabstracttype(t) && t <: family
            ];
            by = nameof,
        )
        for (name, family, f) in cases
            @testset "$name/$(nameof(T))" for T in concrete_laws(family)
                (nameof(T) in vector_laws || (name === :creep_τII && T === NonLinearPeierlsCreep)) && continue  # no compute_τII method
                test_gpu(f, (f32(nd(T())), a))
            end
        end
    end

    @testset "dislocation database" begin
        D = GeoParams.Dislocation
        laws = [
            getfield(D, n) for n in names(D; all = true)
                if isdefined(D, n) && getfield(D, n) isa Function && !startswith(string(n), "#") &&
                !(n in (:eval, :include)) && hasmethod(getfield(D, n), Tuple{})
        ]
        # At 1000 K, exp(-Q/RT) of the dry high-activation-energy laws lies below the
        # normal Float32 range, which Metal flushes to zero; their strain rates there
        # (~1e-29 /s) are physically zero anyway.
        aT = merge(a, (; T = Float32(nd(1200.0K))))
        @testset "$(nameof(law))" for law in laws
            l = f32(nd(SetDislocationCreep(law)))
            test_gpu(c_εII, (l, aT))
            test_gpu(c_τII, (l, aT))
        end
    end

    @testset "composite rheology and phases" begin
        vep = CompositeRheology(SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003), ConstantElasticity(), DruckerPrager())
        test_gpu(c_visc_εII, (f32(nd(vep)), a))
        test_gpu(c_visc_τII, (f32(nd(vep)), a))
        phases = f32(nd((SetMaterialParams(; Phase = 1, Density = PT_Density()), SetMaterialParams(; Phase = 2, Density = MeltDependent_Density()))))
        test_gpu((p, a) -> compute_density(p, 2, a), (phases, a))
        test_gpu((p, a) -> compute_density(p, (0.4f0, 0.6f0), a), (phases, a))
    end

    @testset "softening" begin
        test_gpu((s, x, m) -> s(x, m), (f32(LinearSoftening(0.0, 1.0, 0.0, 1.0)), 0.5f0, 1.0f0))
        test_gpu((s, x) -> s(x), (f32(NonLinearSoftening(; ξ₀ = 1.0, Δ = 0.5)), 0.8f0))
        test_gpu((s, x, m) -> s(x, m), (f32(DecaySoftening()), 1.0f-12, 1.0f0))
    end

    @testset "plastic flow" begin
        τ3, τ6 = (1.0f0, 2.0f0, 3.0f0), ntuple(Float32, 6)
        for p in (DruckerPrager(), DruckerPrager_regularised())
            test_gpu((p, τ) -> CR.∂Q∂τxx(p, τ), (f32(nd(p)), τ3))
            test_gpu((p, τ) -> CR.∂Q∂τyz(p, τ), (f32(nd(p)), τ6))
            test_gpu((p, P) -> CR.∂F∂P(p, P), (f32(nd(p)), a.P))
        end
        test_gpu((p, τ, P) -> CR.∂Q∂τxx(p, τ; P), (f32(nd(DruckerPragerCap())), τ3, a.P))
        mp = f32(nd(SetMaterialParams(; Phase = 1, Plasticity = DruckerPrager(), CreepLaws = LinearViscous())))
        test_gpu((m, τ) -> CR.∂Q∂τ(m, τ), (mp, τ3))
        test_gpu((m, τ, λ) -> CR.plastic_strain_rate(m, τ, λ), (mp, τ3, a.εII))
    end

    @testset "derivatives inside kernels" begin
        hirth = f32(nd(SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003)))
        test_gpu((l, a) -> ForwardDiff.derivative(t -> compute_εII(l, t, a), a.τII), (hirth, a))
        test_gpu((l, a) -> dεII_dτII(l, a.τII, a), (hirth, a))
        costa = f32(nd(ViscosityPartialMelt_Costa_etal_2009()))
        test_gpu((l, a) -> dεII_dτII(l, a.τII, a), (costa, a))
    end

    @testset "solubility and seismic velocity" begin
        for s in (Liu2005_Solubility(), Mafic_Solubility())
            test_gpu((s, P, T, X) -> compute_dissolved(s, P, T, X), (f32(s), 2.0f8, 1200.0f0, 0.3f0))
            test_gpu((s, P, T, X) -> ∂dissolved_∂P(s, P, T, X), (f32(s), 2.0f8, 1200.0f0, 0.3f0))
        end
        test_gpu(l -> compute_wave_velocity(l, (; wave = :Vs)), (f32(nd(ConstantSeismicVelocity())),))
    end

    @testset "laws holding arrays" begin
        td = joinpath(@__DIR__, "test_data")
        for pd in (PerpleX_LaMEM_Diagram(joinpath(td, "Peridotite_dry.in")), MAGEMin_Diagram(joinpath(td, "MAGEMin_Rhyolite.in")))
            test_gpu((pd, T, P) -> compute_density(pd, (; T, P)), (f32(pd), 1500.0f0, 1.0f9))
            test_gpu((ph, T, P) -> compute_density(ph, 1, (; T, P)), ((f32(SetMaterialParams(; Phase = 1, Density = pd)),), 1500.0f0, 1.0f9))
        end
        test_gpu((v, i) -> compute_density(v, (; index = i)), (Vector_Density(; rho = Float32[2900, 3000]), 2))
        test_gpu((v, i) -> compute_heatcapacity(v, (; index = i)), (Vector_HeatCapacity(; Cp = Float32[1000, 1100]), 2))
        test_gpu((v, i) -> compute_meltfraction(v, (; index = i)), (Vector_MeltingParam(; ϕ = Float32[0.1, 0.2]), 2))
    end
end
