using Test
using GeoParams
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
            # these two take τII_old/dt positionally, ahead of an args tuple
            @test typeof(dτII_dεII(el, zero(T), one(T), NamedTuple())) === T
            @test typeof(dp_dεvol(el, zero(T), one(T), NamedTuple())) === T
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
        law = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
        @test compute_εII(law, 1.0e6; T = 1200.0, P = 1.0e9) === 1.3638005232835048e-18
        @test compute_density(PT_Density(), (; P = 1.0e9, T = 1200.0)) === 5695.599999999999
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
        law = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
        τ, args = 1.0f6, (; T = 1200.0, P = 1.0e9)
        τv = fill(τ, 4)
        εv = similar(τv)
        ε = compute_εII(law, τ, args)
        compute_εII!(εv, law, τv)
        compute_τII!(τv, law, εv)
        @test @allocated(compute_εII(law, τ, args)) == 0
        @test @allocated(compute_τII(law, ε, args)) == 0
        @test @allocated(compute_εII!(εv, law, τv)) == 0
        @test @allocated(compute_τII!(τv, law, εv)) == 0
    end
end
