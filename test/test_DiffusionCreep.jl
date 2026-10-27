using Test
using GeoParams
using LaTeXStrings
import GeoParams.Diffusion

@testset "DiffusionCreepLaws" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1Pa * s, length = 1m)

    # Define a linear viscous creep law ---------------------------------
    x1 = DiffusionCreep()
    @test isbits(x1)
    @test param_info(x1).Equation === L"$\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}$"
    str = sprint(show, x1)
    @test occursin("DiffusionCreep:", str)
    @test Value(x1.p) == -3.0
    @test Value(x1.A) == 1.5e6Pa^-1 * s^-1 * m^3

    # perform a computation with the dislocation creep laws
    # Calculate EpsII, using a set of pre-defined values
    CharDim = GEO_units(;
        length = 1000km, viscosity = 1.0e19Pa * s, stress = 100MPa, temperature = 1000C
    )
    EpsII = GeoUnit(1.0s^-1.0)
    TauII = GeoUnit(0.3MPa)
    P = GeoUnit(1.0e9Pa)
    T = GeoUnit(1400C)
    f = GeoUnit(1000NoUnits)
    f_nd = nondimensionalize(f, CharDim)
    d = GeoUnit(10mm)
    d_nd = nondimensionalize(d, CharDim)

    # compute a pure diffusion creep rheology
    diffusion_law = Diffusion.dry_anorthite_Rybacki_2006
    p = SetDiffusionCreep(diffusion_law; n = 1NoUnits)
    @test isbits(p)
    @test p.n.val == 1.0

    # compute a pure diffusion creep rheology
    diffusion_law = Diffusion.dry_anorthite_Rybacki_2006
    p = SetDiffusionCreep(diffusion_law)
    T = 650 + 273.15
    args = (; T = T)
    TauII = 1.0e6
    ε = compute_εII(p, TauII, args)
    @test ε ≈ 1.4177666788096435e-25

    # same but while removing the tensor correction
    ε_notensor = compute_εII(remove_tensor_correction(p), TauII, args)
    @test ε_notensor ≈ 9.451777858730959e-26

    # test with arrays
    τII_array = ones(10) * 1.0e6
    ε_array = similar(τII_array)
    T_array = ones(size(τII_array)) * (650.0 + 273.15)
    args_array = (; T = T_array)
    compute_εII!(ε_array, p, τII_array, args_array)
    @test ε_array[1] ≈ ε

    # compute when args are scalars
    compute_εII!(ε_array, p, τII_array, args)
    @test ε_array[1] ≈ ε

    #---------------------------
    # This is data from a matlab script implementation of the rheology (which was again benchmarked vs. LaMEM)
    for itest in 1:2
        eII = 1.0e-22
        PPa = 0.0
        gsiz = 100
        TK = 650 + 273.15
        logA = [12.1 12.7]    # Logarithm of pre-exponential factor
        npow = [1 3]          # Power law exponent
        Qact = [460 641]      # Activation energy (KJ)
        m_gr0 = [3 0]          # Grain size Exponent (will convert to negative)
        r_fug = [0 0]          # Exponent of Fugacity
        Vact = [24 24]        # Activation Volume cm-3
        fugH = [1 1]          # Fugacity of water MPa
        R = 8.3145         # Gas Constant
        MPa2Pa = 1.0e6            # MPa  -> Pa
        cm32m3 = 1.0e-6           # cm3  -> m3
        J2kJ = 1.0e-3           # Joul -> kJoule
        A0 = 10.0 .^ (logA) # in MPa^-n s^-1 micrometer^m
        m_gr = -m_gr0
        PMPa = PPa / MPa2Pa

        i_flow = itest
        FG_e =
            1 / (
            2^((npow[i_flow] - 1) ./ npow[i_flow]) *
                3^((npow[i_flow] + 1) ./ (2 * npow[i_flow]))
        )
        FG_s = 1 / (3^((npow[i_flow] + 1) ./ 2))

        mu1 =
            FG_e .* eII .^ (1 / npow[i_flow] - 1) *
            A0[i_flow]^(-1.0 / npow[i_flow]) *
            gsiz^(-m_gr[i_flow] / npow[i_flow]) *
            fugH[i_flow]^(-r_fug[i_flow] / npow[i_flow]) *
            exp(
            (Qact[i_flow] + PMPa * MPa2Pa .* Vact[i_flow] * cm32m3 * J2kJ) /
                (R * J2kJ * TK * npow[i_flow]),
        )
        mu = mu1 .* MPa2Pa #In Pa.s
        Tau = 2 * mu * eII     # stress

        #---------------------------

        # Do the same but using GeoParams:
        pp = if itest == 1
            SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006; d = gsiz * 1.0e-6m)
        elseif itest == 2
            SetDislocationCreep(GeoParams.Dislocation.dry_anorthite_Rybacki_2006)
        end

        # using SI units
        τII = compute_τII(pp, eII / s, (; T = TK * K))
        η = τII / (2 * eII / s)
        @test Tau ≈ ustrip(τII)
        @test mu ≈ ustrip(η)

        εII = compute_εII(pp, τII, (; T = TK * K))
        @test eII ≈ ustrip(εII)

        # using Floats
        τII = compute_τII(pp, eII, (; T = TK))
        η = τII / (2 * eII)
        @test Tau ≈ τII
        @test mu ≈ η

        εII = compute_εII(pp, τII, (; T = TK))
        @test eII ≈ ustrip(εII)

        # using arrays for some of the variables
        TK_vec = ones(10) .* TK
        eII_vec = ones(size(TK_vec)) * eII
        τII_vec = zero(eII_vec)
        args = (; T = TK_vec)
        # gsiz_vec = ones(size(TK_vec)) * gsiz * 1.0e-6
        args1 = (; T = TK_vec)

        compute_τII!(τII_vec, pp, eII_vec, args)
        η_vec = τII_vec ./ (2 * eII_vec)
        @test Tau ≈ τII_vec[1]
        @test mu ≈ η_vec[1]

        εII_vec = zero(τII_vec)
        compute_εII!(εII_vec, pp, τII_vec, args)
    end

    # # test overriding the default values
    # a = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)"; V=1e-6m^3 / mol)
    # @test Value(a.V) == 1e-6m^3 / mol

    # Do some basic checks on all creeplaws in the DB
    CharDim = GEO_units()
    creeplaw_list = diffusion_law_list()

    for fun in creeplaw_list
        p = SetDiffusionCreep(fun)                    # original creep law
        p_nd = nondimensionalize(p, CharDim)          # non-dimensionalized
        @test p_nd == SetDiffusionCreep(fun, CharDim) # check that the non-dimensionalized version is the same as the original
        p_dim = dimensionalize(p, CharDim)            # dimensionalized

        # Check that values are the same after non-dimensionalisation & dimensionalisation
        for field in fieldnames(typeof(p_dim))
            val_original = getfield(p, field)
            val_final = getfield(p_dim, field)
            if isa(val_original, GeoUnit)
                @test Value(val_original) == Value(val_final)
            end
        end

        # Perform computations with the rheology
        args = (T = 900.0, d = 100.0e-6, τII_old = 1.0e6)
        ε = 1.0e-15
        τ = compute_τII(p, ε, args)
        ε_test = compute_εII(p, τ, args)
        @test ε ≈ ε_test
    end

    # param_info resolves the builder from the stored Name (find_creep_law)
    info_named = param_info(SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006))
    @test info_named isa MaterialParamsInfo
    @test !isempty(info_named.BibTex_Reference)
    @test isempty(param_info(DiffusionCreep()).BibTex_Reference)               # empty-name early return
    @test isempty(param_info(DiffusionCreep(; Name = "nope_xyz")).BibTex_Reference)  # unknown -> fallback

    # grain size: stored value by default, a runtime `d` overrides it
    p = SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006)       # stored d = 5 mm, p = -3
    p1 = SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006; d = 1mm)
    @test Value(p1.d) ≈ 1.0e-3m                                       # stored in m
    args = (; T = 1100.0, P = 5.0e8)
    τ0, ε0 = 1.0e6, 1.0e-15
    ε_stored = compute_εII(p, τ0, args)
    @test ε_stored ≈ compute_εII(p, τ0, (; args..., d = 5.0e-3))       # omitting d == stored d
    @test compute_εII(p, τ0, (; args..., d = 1.0e-3)) ≈ 125 * ε_stored  # runtime d overrides
    @test compute_εII(p, τ0, (; args..., d = 1.0e-3)) ≈ ustrip(compute_εII(p1, τ0 * Pa; T = 1100.0K, P = 5.0e8Pa))
    τ_stored = compute_τII(p, ε0, args)
    @test τ_stored ≈ compute_τII(p, ε0, (; args..., d = 5.0e-3))
    @test compute_τII(p, ε0, (; args..., d = 1.0e-3)) ≈ τ_stored / 125
    @test dεII_dτII(p, τ0, (; args..., d = 1.0e-3)) ≈ 125 * dεII_dτII(p, τ0, args)
    @test dτII_dεII(p, ε0, (; args..., d = 1.0e-3)) ≈ dτII_dεII(p, ε0, args) / 125

    # unit-ful methods, d = 1mm
    uargs = (; T = 1100.0K, P = 5.0e8Pa)
    @test ustrip(compute_εII(p, τ0 * Pa, uargs)) ≈ ε_stored
    @test ustrip(compute_εII(p, τ0 * Pa; uargs..., d = 1mm)) ≈ 125 * ε_stored
    @test ustrip(compute_εII(p, τ0 * Pa, (; uargs..., d = 1mm))) ≈ 125 * ε_stored
    @test ustrip(Pa, compute_τII(p, ε0 / s; uargs..., d = 1mm)) ≈ τ_stored / 125
    @test ustrip(dεII_dτII(p, τ0 * Pa; uargs..., d = 1mm)) ≈ 125 * dεII_dτII(p, τ0, args)
    @test ustrip(dτII_dεII(p, ε0 / s; uargs..., d = 1mm)) ≈ dτII_dεII(p, ε0, args) / 125

    # arrays: runtime d as array, scalar, or the stored value
    τv = fill(τ0, 3)
    εv = zeros(3)
    compute_εII!(εv, p, τv, args)
    @test all(εv .≈ ε_stored)
    compute_εII!(εv, p, τv, (; args..., d = [5.0e-3, 1.0e-3, 5.0e-3]))
    @test εv ≈ [1, 125, 1] .* ε_stored
    τv2 = zeros(3)
    compute_τII!(τv2, p, fill(ε0, 3); args..., d = 1.0e-3)
    @test all(τv2 .≈ τ_stored / 125)

    # allocation-free and type stable, with and without runtime d
    f_stored(p, τ, T) = compute_εII(p, τ; T = T)
    f_rt(p, τ, T, d) = compute_εII(p, τ; T = T, d = d)
    f_stored(p, τ0, 1100.0)
    f_rt(p, τ0, 1100.0, 1.0e-3)
    @test (@allocated f_stored(p, τ0, 1100.0)) == 0
    @test (@allocated f_rt(p, τ0, 1100.0, 1.0e-3)) == 0
    @test (@inferred f_rt(p, τ0, 1100.0, 1.0e-3)) isa Float64
    @test f_rt(p, Float32(τ0), 1100.0f0, 1.0f-3) isa Float32

    # composite rheology forwards a runtime d
    c = CompositeRheology(p)
    @test compute_εII(c, τ0, (; args..., d = 1.0e-3)) ≈ 125 * compute_εII(c, τ0, args)
end

# Corrected database entries, each against its flow law written in the paper's units:
# ε = A σⁿ d^p fH2O^r exp(-(Q + PV)/RT) with σ (differential stress) and fH2O in MPa, d in μm,
# V in cm³/mol and P in MPa. GeoParams takes τII, P, fH2O and d in SI; σ and ε go through FT and FE.
@testset "Diffusion creep database vs. paper units" begin
    σ, T, P, fH2O, d = 10.0, 1473.0, 300.0, 200.0, 10.0
    for (law, A, n, p, r, Q, V) in (
            (Diffusion.dry_olivine_Hirth_2003, 1.5e9, 1.0, -3.0, 0.0, 375.0e3, 10.0),   # HK03 Table 1
            (Diffusion.wet_olivine_Hirth_2003, 2.5e7, 1.0, -3.0, 1.0, 375.0e3, 20.0),
            (Diffusion.wet_olivine_Mei_2000a, 10^4.7, 1.0, -3.0, 1.0, 295.0e3, 20.0),
            (Diffusion.wet_anorthite_rybacki_2006, 10^-0.7, 1.0, -3.0, 1.0, 159.0e3, 38.0), # Rybacki 2006 Table 5
            (Diffusion.wet_clinopyroxene_HierMajumder2005, 10^7.9 / 64.9e3, 1.0, -3.0, 1.4, 340.0e3, 14.0),
            (Diffusion.dry_olivine_Faul_2007, 2.0e10, 1.37, -3.0, 0.0, 484.0e3, 0.0),     # Faul & Jackson 2007
            (Diffusion.wet_diopside_Dimanov_2005, 10^8.1, 1.0, -3.0, 0.0, 337.0e3, 0.0),
        )
        c = SetDiffusionCreep(law)
        ε = compute_εII(c, σ * 1.0e6 / c.FT; T, P = P * 1.0e6, f = fH2O * 1.0e6, d = d * 1.0e-6) * c.FE
        @test ε ≈ A * σ^n * d^p * fH2O^r * exp(-(Q + V * P) / (NumValue(c.R) * T)) rtol = 1.0e-6
    end

    # HK03 Table 1, constant C_OH: f is C_OH in H/10^6 Si (dimensionless). Footnote c:
    # 1.0e6 (0.3) (10000)^-3 (1000) exp(-(335000 + 1e9 * 4e-6)/(8.314 * 1673)) = 7.8e-15 s^-1
    c = SetDiffusionCreep(Diffusion.wet_olivine1_Hirth_2003)
    ε = compute_εII(c, 0.3e6 / c.FT; T = 1673.0, P = 1.0e9, f = 1000.0, d = 1.0e-2) * c.FE
    @test ε ≈ 1.0e6 * 0.3 * 1.0e4^-3 * 1000 * exp(-(335.0e3 + 1.0e9 * 4.0e-6) / (NumValue(c.R) * 1673.0)) rtol = 1.0e-6
    @test ε ≈ 7.8e-15 rtol = 1.0e-2
end
