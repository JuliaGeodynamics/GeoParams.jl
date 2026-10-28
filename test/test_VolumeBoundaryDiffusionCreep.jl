using Test
using GeoParams
import GeoParams.VolumeBoundaryDiffusion

# Goldsby & Kohlstedt (2001), eq. (4): axial strain rate for differential stress σ.
function gk2001_eq4(σ, T, d; D0v = 9.1e-4, Qv = 59.4e3, D0b = 6.0e-4, Qb = 49.0e3, δ = 9.04e-10, Vm = 1.97e-5, R = 8.3145)
    Dv = D0v * exp(-Qv / (R * T))
    Db = D0b * exp(-Qb / (R * T))
    return 42 * σ * Vm / (R * T * d^2) * (Dv + π * δ / d * Db)
end

@testset "VolumeBoundaryDiffusionCreep" begin
    p = SetVolumeBoundaryDiffusionCreep(VolumeBoundaryDiffusion.ice_Goldsby_2001)
    @test isbits(p)
    @test p isa VolumeBoundaryDiffusionCreep
    @test param_info(p).Comment != ""
    @test occursin("10.1029/2000JB900336", param_info(p).BibTex_Reference)
    @test VolumeBoundaryDiffusion.ice_Goldsby_2001 in volumeboundarydiffusion_law_list()

    T, d, σ = 250.0, 1.0e-3, 0.1e6

    @testset "printed formula" begin
        # without tensor correction the law is eq. (4) verbatim
        @test compute_εII(remove_tensor_correction(p), σ; T, d) ≈ gk2001_eq4(σ, T, d) rtol = 1.0e-12
        # axial compression: τII = σ/√3, ε̇II = (√3/2) ε̇, hence η = RTd²/(126 Vm D)
        τII = σ / sqrt(3)
        @test compute_εII(p, τII; T, d) ≈ sqrt(3) / 2 * gk2001_eq4(σ, T, d) rtol = 1.0e-12
        D = 9.1e-4 * exp(-59.4e3 / (8.3145T)) + π * 9.04e-10 / d * 6.0e-4 * exp(-49.0e3 / (8.3145T))
        η = 8.3145T * d^2 / (126 * 1.97e-5 * D)
        @test compute_εII(p, τII; T, d) ≈ τII / (2η) rtol = 1.0e-12
        # the stored grain size is the default
        @test compute_εII(p, τII; T) == compute_εII(p, τII; T, d = 1.0e-3)
        @test compute_εII(p, τII; T, d = 1.0e-4) > compute_εII(p, τII; T)
        # Unitful path
        @test ustrip(s^-1, compute_εII(p, τII * Pa; T = T * K, d = d * m)) ≈ compute_εII(p, τII; T, d)
    end

    @testset "round trip and derivatives" begin
        args = (; T, d = 2.0e-4)
        τ = compute_τII(p, 1.0e-12, args)
        @test compute_εII(p, τ, args) ≈ 1.0e-12
        @test dεII_dτII(p, τ, args) ≈ 1.0e-12 / τ
        @test dτII_dεII(p, 1.0e-12, args) ≈ τ / 1.0e-12
        @test ustrip(Pa, compute_τII(p, 1.0e-12 / s; T = T * K, d = 2.0e-4m)) ≈ τ

        τv = [1.0e4, 1.0e5, 1.0e6]
        εv = similar(τv)
        compute_εII!(εv, p, τv; T = [240.0, 250.0, 260.0], d)
        @test εv ≈ [compute_εII(p, τv[i]; T = 230.0 + 10i, d) for i in 1:3]
        τ2 = similar(τv)
        compute_τII!(τ2, p, εv; T = [240.0, 250.0, 260.0], d)
        @test τ2 ≈ τv
    end

    @testset "nondimensionalization" begin
        CharDim = GEO_units()
        p_nd = SetVolumeBoundaryDiffusionCreep(VolumeBoundaryDiffusion.ice_Goldsby_2001, CharDim)
        @test all(k -> isfinite(NumValue(getfield(p_nd, k))), (:d, :D0v, :Qv, :D0b, :Qb, :δ, :Vm, :R))
        τn = nondimensionalize(1.0e5Pa, CharDim)
        Tn = nondimensionalize(T * K, CharDim)
        εn = compute_εII(p_nd, τn; T = Tn)
        @test ustrip(dimensionalize(εn, s^-1, CharDim)) ≈ compute_εII(p, 1.0e5; T) rtol = 1.0e-10
        @test compute_τII(p_nd, εn; T = Tn) ≈ τn
        @test Transform_VolumeBoundaryDiffusionCreep(VolumeBoundaryDiffusion.ice_Goldsby_2001, CharDim) == p_nd
    end

    @testset "Float32" begin
        p32 = convert_precision(Float32, p)
        for law in (p, p32), Tk in (100.0f0, 250.0f0)
            ε = compute_εII(law, 1.0f5; T = Tk)
            @test ε isa Float32
            @test ε ≈ compute_εII(p, 1.0e5; T = Float64(Tk)) rtol = 1.0e-5
            @test compute_τII(law, ε; T = Tk) ≈ 1.0f5 rtol = 1.0e-5
        end
        # At 50 K the Float32 coefficient underflows to zero; every routine recomputes in Float64
        args = (; T = 50.0f0)
        @test compute_εII(p, 1.0f30, args) ≈ Float32(compute_εII(p, 1.0e30; T = 50.0))
        @test dεII_dτII(p, 1.0f8, args) == Float32(dεII_dτII(p, 1.0e8; T = 50.0))
        @test compute_τII(p, 1.0f-30, args) == Float32(compute_τII(p, 1.0e-30; T = 50.0))
        @test dτII_dεII(p, 1.0f-30, args) == Float32(dτII_dεII(p, 1.0e-30; T = 50.0))
    end

    @testset "CompositeRheology" begin
        c = CompositeRheology(p, SetDislocationCreep(GeoParams.Dislocation.ice_dislocation_Goldsby_2001), ConstantElasticity())
        args = (; T = 265.0, d, τII_old = 1.0e4, dt = 1.0e9)   # the dislocation law holds above 258 K
        for ε in (1.0e-14, 1.0e-12, 1.0e-10)
            τ = compute_τII(c, ε, args)
            @test compute_εII(c, τ, args) ≈ ε rtol = 1.0e-6
        end
        # without elasticity the composite is the harmonic sum of the two viscous laws
        c2 = CompositeRheology(p, LinearViscous(; η = 1.0e14Pa * s))
        τ = compute_τII(c2, 1.0e-10, (; T, d))
        @test compute_εII(p, τ; T, d) + τ / 2.0e14 ≈ 1.0e-10 rtol = 1.0e-6
    end
end
