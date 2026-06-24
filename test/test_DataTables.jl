using Test
using GeoParams
import GeoParams.Dislocation, GeoParams.Diffusion, GeoParams.GBS, GeoParams.Peierls, GeoParams.NonLinearPeierls
import GeoParams: Rutile, Garnet, Olivine, Melt

# Enumerates the (un-exported) law-builder functions in a creep `Data` module that,
# when called with no arguments, return a `(law::T, MaterialParamsInfo)` tuple.
function _creep_laws(mod, ::Type{T}) where {T}
    laws = Any[]
    for name in names(mod; all = true, imported = false)
        s = string(name)
        (startswith(s, "#") || name in (:eval, :include)) && continue
        f = getfield(mod, name)
        isa(f, Function) || continue
        try
            r = f()
            if isa(r, Tuple) && length(r) == 2 && r[1] isa T
                push!(laws, f)
            end
        catch
        end
    end
    return laws
end

@testset "Data tables sweep" begin
    CharDim = GEO_units(; viscosity = 1.0e20, length = 1000km)
    args = (; T = 1373.0, P = 1.0e9, f = 1.0, d = 1.0e-3)

    # ----- viscous creep-law databases -----
    creep_cases = (
        (Dislocation, DislocationCreep, SetDislocationCreep, true),
        (Diffusion, DiffusionCreep, SetDiffusionCreep, true),
        (GBS, GrainBoundarySliding, SetGrainBoundarySliding, true),
        (Peierls, PeierlsCreep, SetPeierlsCreep, false),  # inverse compute_τII is domain-limited; covered in test_PeierlsCreep
        (NonLinearPeierls, NonLinearPeierlsCreep, SetNonLinearPeierlsCreep, false), # no closed-form compute_τII
    )

    for (mod, T, setfun, has_τII) in creep_cases
        laws = _creep_laws(mod, T)
        @testset "$(nameof(mod)) ($(length(laws)) laws)" begin
            @test !isempty(laws)
            for f in laws
                # builds the law (covers the data-file body) + Transform to SI units
                p = setfun(f)
                @test p isa T
                # non-dimensional transform branch
                @test setfun(f, CharDim) isa T
                # show + remove_tensor_correction
                @test sprint(show, p) isa String
                @test remove_tensor_correction(p) isa T
                # computational routines (plain-float path)
                ε = compute_εII(p, 1.0e6, args)
                @test ε isa Number && !isnan(ε)
                @test dεII_dτII(p, 1.0e6, args) isa Number
                if has_τII
                    τ = compute_τII(p, 1.0e-15, args)
                    @test τ isa Number && !isnan(τ)
                    @test dτII_dεII(p, 1.0e-15, args) isa Number
                end
            end
        end
    end

    # ----- chemical-diffusion databases -----
    @testset "ChemicalDiffusion data" begin
        for mod in (Rutile, Garnet, Olivine, Melt)
            fns = mod.chemical_diffusion_list()
            ndiff = 0
            nmulti = 0
            for f in fns
                isa(f, Function) || continue
                local r
                try
                    r = f()
                catch
                    continue
                end
                (isa(r, Tuple) && length(r) == 2) || continue
                data = r[1]
                if data isa DiffusionData
                    ndiff += 1
                    p = SetChemicalDiffusion(f)
                    @test p isa DiffusionData
                    @test sprint(show, p) isa String
                    # plain-float and dimensional (units) paths of compute_D
                    @test compute_D(p; T = 1273.0, P = 1.0e9, X = 0.1, fO2 = 1.0e-7) isa Number
                    @test compute_D(p; T = 1273.0K, P = 1.0e9Pa, X = 0.1NoUnits, fO2 = 1.0e-7NoUnits) isa Quantity
                    # in-place path
                    D_arr = zeros(3)
                    compute_D!(D_arr, p; T = fill(1273.0, 3), P = fill(1.0e9, 3), X = fill(0.1, 3), fO2 = fill(1.0e-7, 3))
                    @test all(isfinite, D_arr)
                elseif data isa MeltMulticompDiffusionData
                    nmulti += 1
                    p = SetMulticompChemicalDiffusion(f)
                    @test p isa MeltMulticompDiffusionData
                    @test sprint(show, p) isa String
                    D = compute_D(p; T = 1500.0)
                    @test size(D, 1) ≥ 1
                end
            end
            @test ndiff + nmulti > 0
        end

        # missing from test_ChemicalDiffusion.jl: Rt_Zr_Sasaki1985_para_c
        Zr_para = SetChemicalDiffusion(Rutile.Rt_Zr_Sasaki1985_para_c)
        D = ustrip(compute_D(Zr_para, T = 1200C))
        @test D isa Float64 && D > 0
    end
end
