using Test
using GeoParams
import GeoParams: Dislocation, Diffusion, GBS, NonLinearPeierls

# Reference values from the script that made the creep-law figures of Becker & Faccenna, Tectonic
# Geodynamics (draft, 2026), chapter 5 (calc_creep_visc.awk; input "P T d ε̇ f" in SI units, output
# "viscosity stress", or "viscosity strain rate" for calc_strainrate = 1).
# They check the stress conventions documented in docs/src/man/rheology_becker_faccenna.md:
#  - Table 5.2 laws ("ol1"):  viscosity = τII / (2ε̇II) of GeoParams, printed stress = 2τII
#  - all other laws (stress scale B = 2^((1+3/n)/2) A^(-1/n)): ε̇ = A σ^n / 2^((n+3)/2), i.e. at the same
#    differential stress σd = √3 τII the GeoParams ε̇II is √3 2^((n+1)/2) times larger.

const pts = ((1.0e9, 1523.15, 3.0e-3, 1.0e-15), (5.0e8, 873.15, 1.0e-3, 1.0e-13), (2.0e9, 1173.15, 1.0e-2, 1.0e-12))
k_B(n) = sqrt(3) * 2^((n + 1) / 2)

@testset "Becker & Faccenna textbook conventions" begin
    # Table 5.2 dry dislocation and diffusion creep ("ol1"): equal viscosity, stress = 2τII
    dis = SetDislocationCreep(Dislocation.dry_olivine_Hirth_2003)
    dif = SetDiffusionCreep(Diffusion.dry_olivine_Hirth_2003; V = 6.0e-6m^3 / mol)
    awk_dis = ((1.20712962e20, 4.82851849e5), (3.17274757e22, 1.26909903e10), (5.09550867e19, 2.03820347e8))
    awk_dif = ((6.97861306e19, 2.79144522e5), (9.10735562e27, 3.64294225e15), (3.78070507e25, 1.51228203e14))
    for ((P, T, d, ε̇), (η_d, s_d), (η_f, s_f)) in zip(pts, awk_dis, awk_dif)
        τ = compute_τII(dis, ε̇; T, P)
        @test τ / (2ε̇) ≈ η_d rtol = 2.0e-3          # Table 5.2 rounds A to 7.4e-15 (7.37e-15)
        @test 2τ ≈ s_d rtol = 2.0e-3
        τ = compute_τII(dif, ε̇; T, P, d)
        @test τ / (2ε̇) ≈ η_f rtol = 1.0e-3
        @test 2τ ≈ s_f rtol = 1.0e-3
    end

    # B-convention power laws: printed stress = (√3 2^((n+1)/2))^(1/n) σd
    gbs_cold = SetGrainBoundarySliding(GBS.cold_dry_olivine_Hirth_2003)
    gbs_hot = SetGrainBoundarySliding(GBS.hot_dry_olivine_Hirth_2003)
    tokle = SetDislocationCreep(Dislocation.wet_quartzite_Tokle_2019)
    rybacki = SetDislocationCreep(Dislocation.wet_anorthite_Rybacki_2006)
    fw = (3.0e8, 1.0e8, 1.0e9)                       # water fugacity [Pa]
    awk_gbs = (9.69898236e6, 1.51606444e10, 3.92359629e9)
    awk_tokle = (4.08269944e5, 2.79939699e7, 6.2227951e6)
    awk_rybacki = (6.11961603e4, 3.10978457e8, 3.03783155e7)
    for (i, (P, T, d, ε̇)) in enumerate(pts)
        gbs = T < 1523.15 ? gbs_cold : gbs_hot
        @test k_B(3.5)^(1 / 3.5) * sqrt(3) * compute_τII(gbs, ε̇; T, P, d) ≈ awk_gbs[i] rtol = 1.0e-3
        @test k_B(3.0)^(1 / 3) * sqrt(3) * compute_τII(tokle, ε̇; T, P, f = fw[i]) ≈ awk_tokle[i] rtol = 1.0e-3
        @test k_B(3.0)^(1 / 3) * sqrt(3) * compute_τII(rybacki, ε̇; T, P, f = fw[i]) ≈ awk_rybacki[i] rtol = 1.0e-3
    end

    # Peierls creep of Mei et al. (2010), strain rate at given differential stress
    mei = SetNonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010)
    for (σ, T, ε̇_awk) in ((1.0e9, 1073.15, 1.70014569e-11), (5.0e8, 873.15, 1.66446313e-16), (2.0e9, 1273.15, 3.24178376e-7))
        @test compute_εII(mei, σ / sqrt(3); T) / k_B(2) ≈ ε̇_awk rtol = 1.0e-3
    end
end
