using Test
using GeoParams
using LaTeXStrings
import GeoParams.Peierls

@testset "PeierlsCreep" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity = 1Pa * s, length = 1m)

    # Define a linear viscous creep law ---------------------------------
    x1 = PeierlsCreep()
    @test isbits(x1)
    @test param_info(x1).Equation === L"$\tau_{ij} = 2 \eta  \dot{\varepsilon}_{ij}$"
    str = sprint(show, x1)
    @test occursin("PeierlsCreep:", str)
    @test Value(x1.n) == 1.0
    @test Value(x1.q) == 2.0
    @test Value(x1.A) == 5.7e11s^(-1.0)
    @test repr("text/plain", x1) isa String

    # perform a computation with the dislocation creep laws
    # Calculate EpsII, using a set of pre-defined values
    CharDim = GEO_units(;
        length = 1000km, viscosity = 1.0e19Pa * s, stress = 100MPa, temperature = 1000C
    )
    EpsII = GeoUnit(1.0s^-1.0)
    TauII = GeoUnit(1.0e6MPa)
    T = GeoUnit(600C)

    # compute a pure diffusion creep rheology
    p = SetPeierlsCreep(Peierls.dry_olivine_Goetze_1979)
    T = 600 + 273.15
    args = (; T = T)
    TauII = 1.0e9
    ε = compute_εII(p, TauII, args)
    @test ε ≈ 2.3181145468270706e-9

    # same but while removing the tensor correction
    ε_notensor = compute_εII(remove_tensor_correction(p), TauII, args)
    @test ε_notensor ≈ 6.197048346429745e-14

    # test with arrays
    τII_array = ones(10) * 1.0e9
    ε_array = similar(τII_array)
    T_array = ones(size(τII_array)) * (600.0 + 273.15)
    args_array = (; T = T_array)

    compute_εII!(ε_array, p, τII_array, args_array)
    @test ε_array[1] ≈ ε

    # compute when args are scalars
    compute_εII!(ε_array, p, τII_array, args)
    @test ε_array[1] ≈ ε
    # ===

    # test overriding the default values
    # a = SetPeierlsCreep("Dry Olivine | Goetze and Evans (1979)"; E=535.0kJ / mol)
    # @test Value(a.E) == 535.0kJ / mol

    # Do some basic checks on all creeplaws in the DB
    CharDim = GEO_units()
    creeplaw_list = peierls_law_list()       # all creeplaws in database
    for fun in creeplaw_list
        p = SetPeierlsCreep(fun)                # original creep law
        p_nd = nondimensionalize(p, CharDim)        # non-dimensionalized
        @test p_nd == SetPeierlsCreep(fun, CharDim) # check that the non-dimensionalized version is the same as the original
        p_dim = dimensionalize(p, CharDim)           # dimensionalized

        # Check that values are the same after non-dimensionalisation & dimensionalisation
        for field in fieldnames(typeof(p_dim))
            val_original = getfield(p, field)
            val_final = getfield(p_dim, field)
            if isa(val_original, GeoUnit)
                @test Value(val_original) == Value(val_final)
            end
        end

        # Perform computations with the rheology
        args = (T = 500.0, τII_old = 2.2e9)
        ε = 1.0e-7
        τ = compute_τII(p, ε, args)
        ε_test = compute_εII(p, τ, args)
        @test ε ≈ ε_test
    end

    # param_info resolves the builder from the stored Name (find_creep_law)
    info_named = param_info(SetPeierlsCreep(Peierls.dry_olivine_Goetze_1979))
    @test info_named isa MaterialParamsInfo
    @test !isempty(info_named.BibTex_Reference)
    @test isempty(param_info(PeierlsCreep()).BibTex_Reference)                # empty-name early return
    @test isempty(param_info(PeierlsCreep(; Name = "nope_xyz")).BibTex_Reference)  # unknown -> fallback

    @testset "PeierlsCreep stress paths" begin
        a = SetPeierlsCreep(Peierls.dry_olivine_Goetze_1979)
        # dimensional (Quantity) compute_τII — regression value
        τq = compute_τII(a, 1.0e-15 / s; T = 1000K)
        @test ustrip(uconvert(u"Pa", τq)) ≈ 1.1562073728905828e8 rtol = 1.0e-6
        # array compute_τII! matches the scalar result
        ε = fill(1.0e-15, 4)
        τ = similar(ε)
        compute_τII!(τ, a, ε; T = fill(1000.0, 4))
        @test τ[1] ≈ 1.1562073728905828e8 rtol = 1.0e-6
        # AD derivative — regression value
        @test dτII_dεII(a, 1.0e-15; T = 1000.0) ≈ 3.898110451060781e22 rtol = 1.0e-6
    end

    # Database laws against the papers: ε̇ = A exp(-E/RT (1 - (σ/σP)^p)^q), with σ the differential stress
    # (σ = √3 τII) and ε̇ the axial strain rate (ε̇II = ε̇ / (2/√3)). σ = 1 GPa, T = 873 K.
    σ, T = 1.0e9, 873.0
    for (law, A, E, σP, pp, q) in (
            (Peierls.dry_olivine_Goetze_1979, 5.7e11, 536.0e3, 8.5e9, 1.0, 2.0),     # Goetze & Evans (1979), eq. (1a)
            (Peierls.dry_olivine_Demouchy_2013, 1.0e6, 450.0e3, 15.0e9, 0.5, 2.0),   # Demouchy et al. (2013), eq. (4)
            (Peierls.dry_olivine_Idrissi_2016, 1.0e6, 566.0e3, 3.8e9, 0.5, 2.0),     # Idrissi et al. (2016), eq. (1)
        )
        a = SetPeierlsCreep(law)
        @test compute_εII(a, σ / sqrt(3); T) ≈ A * exp(-E / (NumValue(a.R) * T) * (1 - (σ / σP)^pp)^q) / (2 / sqrt(3)) rtol = 1.0e-6
    end
end
