using Test
using GeoParams
import GeoParams.Diffusion

@testset "DiffusionCreepLaws" begin

    # This tests the MaterialParameters structure
    CharUnits_GEO = GEO_units(; viscosity=1Pa * s, length=1m)

    # Define a linear viscous creep law ---------------------------------
    x1 = DiffusionCreep()
    @test Value(x1.n) == 1.0
    @test Value(x1.p) == -3.0
    @test Value(x1.A) == 1.5MPa^-1 * s^-1 * m^3

    # perform a computation with the dislocation creep laws 
    # Calculate EpsII, using a set of pre-defined values
    CharDim = GEO_units(;
        length=1000km, viscosity=1e19Pa * s, stress=100MPa, temperature=1000C
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
    p = SetDiffusionCreep(diffusion_law)

    T = 650 + 273.15

    args = (; T=T)
    TauII = 1e6
    ε = compute_εII(p, TauII, args)
    @test ε ≈ 1.7722083485120549e-32

    # same but while removing the tensor correction
    ε_notensor = compute_εII(remove_tensor_correction(p), TauII, args)
    @test ε_notensor ≈ 1.1814722323413693e-32

    # test with arrays
    τII_array = ones(10) * 1e6
    ε_array = similar(τII_array)
    T_array = ones(size(τII_array)) * (650.0 + 273.15)

    args_array = (; T=T_array)

    compute_εII!(ε_array, p, τII_array, args_array)
    @test ε_array[1] ≈ ε

    # compute when args are scalars
    compute_εII!(ε_array, p, τII_array, args)
    @test ε_array[1] ≈ ε

    #---------------------------
    # This is data from a matlab script implementation of the rheology (which was again benchmarked vs. LaMEM)
    for itest in 1:2
        eII = 1e-22
        PPa = 0.0
        gsiz = 100
        TK = 650 + 273.15

        logA = [12.1 12.7] #Logarithm of pre-exponential factor
        npow = [1 3] #Power law exponent
        Qact = [460 641] #Activation energy (KJ)
        m_gr0 = [3 0] #Grain size Exponent (will convert to negative)
        r_fug = [0 0] #Exponent of Fugacity
        Vact = [24 24] #Activation Volume cm-3
        fugH = [1 1] #Fugacity of water MPa 

        R = 8.3145 #Gas Constant
        MPa2Pa = 1e6   #MPa  -> Pa
        cm32m3 = 1e-6  #cm3  -> m3
        J2kJ = 1e-3  #Joul -> kJoule

        A0 = 10.0 .^ (logA)  # in MPa^-n s^-1 micrometer^m
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
            SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006)
        elseif itest == 2
            SetDislocationCreep(GeoParams.Dislocation.dry_anorthite_Rybacki_2006)
        end

        # using SI units
        τII = compute_τII(pp, eII / s, (; T=TK * K, d=gsiz * 1e-6m))
        η = τII / (2 * eII / s)
        @test Tau ≈ ustrip(τII)
        @test mu ≈ ustrip(η)

        εII = compute_εII(pp, τII, (; T=TK * K, d=gsiz * 1e-6m))
        @test eII ≈ ustrip(εII)

        # using Floats
        τII = compute_τII(pp, eII, (; T=TK, d=gsiz * 1e-6))
        η = τII / (2 * eII)
        @test Tau ≈ τII
        @test mu ≈ η

        εII = compute_εII(pp, τII, (; T=TK, d=gsiz * 1e-6))
        @test eII ≈ ustrip(εII)

        # using arrays for some of the variables
        TK_vec = ones(10) .* TK
        eII_vec = ones(size(TK_vec)) * eII
        τII_vec = zero(eII_vec)
        args = (; T=TK_vec, d=gsiz * 1e-6)
        gsiz_vec = ones(size(TK_vec)) * gsiz * 1e-6
        args1 = (; T=TK_vec, d=gsiz_vec)

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
        args = (T=900.0, d=100e-6, τII_old=1e6)
        ε = 1e-15
        τ = compute_τII(p, ε, args)
        ε_test = compute_εII(p, τ, args)
        @test ε ≈ ε_test
    end
end