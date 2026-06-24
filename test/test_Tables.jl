using Test
using GeoParams
using Unidecode
import GeoParams: Dislocation, Diffusion, GBS, Peierls, NonLinearPeierls
import GeoParams.Tables: detachFloatfromExponent, extract_parameters_from_phases, Dict2LatexTable, extract_parameters_from_phases_md, Dict2MarkdownTable, ParameterTable, create_latex_symbol, get_material_reference_info
import GeoParams.Tables: get_rheology_info, format_number, format_markdown_number

@testset "Tables.jl" begin
    MatParam = (
        SetMaterialParams(;
            Name = "Viscous Matrix",
            Phase = 1,
            Density = ConstantDensity(),
            CreepLaws = SetDislocationCreep(Dislocation.quartz_diorite_HansenCarter_1982),
        ),
        SetMaterialParams(;
            Name = "Viscous Sinker",
            Phase = 2,
            Density = PT_Density(),
            CreepLaws = LinearViscous(; η = 1.0e21Pa * s),
        ),
        SetMaterialParams(;
            Name = "Viscous Bottom",
            Phase = 3,
            Density = PT_Density(),
            CreepLaws = SetDislocationCreep(Dislocation.diabase_Caristan_1982),
        ),
    )

    # test detachFloatfromExponent()
    n1 = "1"
    n2 = "10.1234"
    n3 = "1.12345e6"
    n4 = "12.34567e-8"
    dig, num, ex = detachFloatfromExponent(n1)
    @test dig == 0
    @test num == "1"
    @test ex == "1"
    dig, num, ex = detachFloatfromExponent(n2)
    @test dig == 4
    @test num == "10.1234"
    @test ex == "1"
    dig, num, ex = detachFloatfromExponent(n3)
    @test dig == 5
    @test num == "1.12345"
    @test ex == "6"
    dig, num, ex = detachFloatfromExponent(n4)
    @test dig == 5
    @test num == "12.34567"
    @test ex == "-8"

    # test extract_parameters_from_phases()
    dict, ref = extract_parameters_from_phases(MatParam)
    dval = MatParam[1].Density[1].ρ.val
    @test dict["ρ Density 1"][1] == "$dval"
    @test dict["ρ Density 1"][2] == "\\" * "$(unidecode("ρ"))"
    @test dict["ρ Density 1"][3] == ""
    @test dict["ρ Density 1"][4] == "1"
    @test dict["ρ Density 1"][5] == ""
    @test dict["ρ Density 1"][6] == ""
    @test dict["Name 1"][1] == "$(unsafe_string(MatParam[1].Name))"
    @test dict["Name 1"][2] == "$(length(MatParam))"
    @test dict["Name 1"][3] == "1"
    @test dict["Name 1"][4] == ""
    @test dict["Name 1"][5] == ""
    @test dict["Name 1"][6] == ""

    # test extract_parameters_from_phases_md()
    dictMd = extract_parameters_from_phases_md(MatParam)
    dvalMd = MatParam[1].Density[1].ρ.val
    @test dictMd["ρ Density 1"][1] == "$dvalMd"
    @test dictMd["ρ Density 1"][2] == "ρ"
    @test dictMd["ρ Density 1"][3] == ""
    @test dictMd["ρ Density 1"][4] == "1"
    @test dictMd["ρ Density 1"][5] == ""
    @test dictMd["ρ Density 1"][6] == ""
    @test dictMd["Name 1"][1] == "$(unsafe_string(MatParam[1].Name))"
    @test dictMd["Name 1"][2] == "$(length(MatParam))"
    @test dictMd["Name 1"][3] == "1"
    @test dictMd["Name 1"][4] == ""
    @test dictMd["Name 1"][5] == ""
    @test dictMd["Name 1"][6] == ""

    # test Dict2LatexTable()
    Dict2LatexTable(dict, ref)
    @test "ParameterTable.tex" in readdir()
    @test "References.bib" in readdir()

    rm("ParameterTable.tex"; force = true)
    rm("References.bib"; force = true)

    # test Dict2MarkdownTable()
    Dict2MarkdownTable(dictMd)
    @test "ParameterTable.md" in readdir()

    rm("ParameterTable.md"; force = true)

    # test ParameterTable()
    # for Latex
    ParameterTable(MatParam)
    @test "ParameterTable.tex" in readdir()
    @test "References.bib" in readdir()
    rm("ParameterTable.tex"; force = true)
    rm("References.bib"; force = true)

    ParameterTable(MatParam; format = "TEX")
    @test "ParameterTable.tex" in readdir()
    @test "References.bib" in readdir()
    rm("ParameterTable.tex"; force = true)
    rm("References.bib"; force = true)

    ParameterTable(MatParam; format = "LaTeX", filename = "TestTable")
    @test "TestTable.tex" in readdir()
    rm("TestTable.tex"; force = true)
    rm("References.bib"; force = true)

    # for Markdown
    ParameterTable(MatParam; format = "Markdown")
    @test "ParameterTable.md" in readdir()
    rm("ParameterTable.md"; force = true)

    ParameterTable(MatParam; format = "MD")
    @test "ParameterTable.md" in readdir()
    rm("ParameterTable.md"; force = true)

    filename = "TestTable"
    ParameterTable(MatParam; format = "MaRkDoWn", filename = filename)
    @test "TestTable.md" in readdir()
    rm("TestTable.md"; force = true)

    # test create_latex_symbol()
    @test create_latex_symbol("ρ") == "\\rho"
    @test create_latex_symbol("rho") == "\\rho"
    @test create_latex_symbol("η0") == "\\eta_0"
    @test create_latex_symbol("T0") == "T_0"
    @test create_latex_symbol("P0") == "P_0"
    @test create_latex_symbol("A_diff") == "A_{\\text{diff}}"
    @test create_latex_symbol("n") == "n"
    @test create_latex_symbol("E") == "E"

    # test get_material_reference_info()
    disl_creep = SetDislocationCreep(Dislocation.quartz_diorite_HansenCarter_1982)
    ref_info = get_material_reference_info(disl_creep)
    @test ref_info !== nothing
    @test hasfield(typeof(ref_info), :BibTex_Reference)

    # test phase with CompositeRheology field
    v1 = SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006)
    c1 = CompositeRheology(
        v1,
        SetDislocationCreep(Dislocation.diabase_Caristan_1982),
        LinearViscous(; η = 1.0e21Pa * s),
        v1,
    )
    MatParam = (
        SetMaterialParams(;
            Name = "Viscous Matrix",
            Phase = 1,
            Density = ConstantDensity(),
            CreepLaws = SetDislocationCreep(Dislocation.quartz_diorite_HansenCarter_1982),
        ),
        SetMaterialParams(;
            Name = "Viscous Sinker", Phase = 2, Density = PT_Density(), CompositeRheology = c1
        ),
        SetMaterialParams(;
            Name = "Viscous Bottom",
            Phase = 3,
            Density = PT_Density(),
            CreepLaws = SetDislocationCreep(Dislocation.diabase_Caristan_1982),
        ),
    )

    # test extract_parameters_from_phases() for CompositeRheology
    dict, ref = extract_parameters_from_phases(MatParam)
    dval = MatParam[2].CompositeRheology[1][3].η.val
    @test dict["η CompositeRheology LinVisc 2.3"][1] == "$dval"
    @test dict["η CompositeRheology LinVisc 2.3"][2] == "\\" * "$(unidecode("η"))"
    @test dict["η CompositeRheology LinVisc 2.3"][3] ==
        "CompoRheo(DiffCreep,DislCreep,LinVisc,)"
    @test dict["η CompositeRheology LinVisc 2.3"][4] == "2"
    @test dict["η CompositeRheology LinVisc 2.3"][5] == "1"
    @test dict["η CompositeRheology LinVisc 2.3"][6] == "LinVisc"

    # test extract_parameters_from_phases_md() for CompositeRheology
    dictMd = extract_parameters_from_phases_md(MatParam)
    dvalMd = MatParam[2].CompositeRheology[1][3].η.val
    @test dictMd["η CompositeRheology LinVisc 2.3"][1] == "$dvalMd"
    @test dictMd["η CompositeRheology LinVisc 2.3"][2] == "η"
    @test dictMd["η CompositeRheology LinVisc 2.3"][3] ==
        "CompoRheo(DiffCreep,DislCreep,LinVisc,)"
    @test dictMd["η CompositeRheology LinVisc 2.3"][4] == "2"
    @test dictMd["η CompositeRheology LinVisc 2.3"][5] == "1"
    @test dictMd["η CompositeRheology LinVisc 2.3"][6] == "LinVisc"

    # ---- get_material_reference_info: GBS / Peierls / NonLinearPeierls branches ----
    v_gbs = SetGrainBoundarySliding(GBS.cold_dry_olivine_Hirth_2003)
    @test get_material_reference_info(v_gbs) isa MaterialParamsInfo

    v_peierls = SetPeierlsCreep(Peierls.dry_olivine_Demouchy_2013)
    @test get_material_reference_info(v_peierls) isa MaterialParamsInfo

    v_nlp = SetNonLinearPeierlsCreep(NonLinearPeierls.dry_olivine_Mei_2010)
    @test get_material_reference_info(v_nlp) isa MaterialParamsInfo

    # LinearViscous has no database entry -> returns nothing (fallback path)
    @test get_material_reference_info(LinearViscous()) === nothing

    # ---- get_rheology_info: cover all named branches + fallback ----
    @test get_rheology_info(SetDislocationCreep(Dislocation.quartz_diorite_HansenCarter_1982)) == ("DislCreep", "DislCreep")
    @test get_rheology_info(SetDiffusionCreep(Diffusion.dry_anorthite_Rybacki_2006)) == ("DiffCreep", "DiffCreep")
    @test get_rheology_info(LinearViscous()) == ("LinVisc", "LinVisc")
    @test get_rheology_info(PowerlawViscous()) == ("PowerVisc", "PowerVisc")
    @test get_rheology_info(v_gbs) == ("GBS", "GBS")
    @test get_rheology_info(v_peierls) == ("PeierlsCreep", "PeierlsCreep")
    @test get_rheology_info(v_nlp) == ("NonLinPeierls", "NonLinPeierls")
    @test get_rheology_info(LinearMeltViscosity()) == ("MeltVisc", "MeltVisc")
    @test get_rheology_info(ViscosityPartialMelt_Costa_etal_2009()) == ("PartialMelt", "PartialMelt")
    @test get_rheology_info(GiordanoMeltViscosity()) == ("GiordanoMelt", "GiordanoMelt")
    # fallback: unknown type -> type name
    r1, r2 = get_rheology_info(ConstantDensity())
    @test r1 isa String && !isempty(r1)

    # ---- ParameterTable with more field types (Elasticity, HeatCapacity, Conductivity) ----
    MatParam3 = (
        SetMaterialParams(;
            Name = "Elastic Phase",
            Phase = 1,
            Density = ConstantDensity(),
            Elasticity = ConstantElasticity(),
            HeatCapacity = ConstantHeatCapacity(),
            Conductivity = ConstantConductivity(),
        ),
    )
    ParameterTable(MatParam3)
    @test "ParameterTable.tex" in readdir()
    rm("ParameterTable.tex"; force = true)
    rm("References.bib"; force = true)

    ParameterTable(MatParam3; format = "Markdown")
    @test "ParameterTable.md" in readdir()
    rm("ParameterTable.md"; force = true)

    # ---- format_number / format_markdown_number: all four branches ----
    # dig ≤ rdigits, expo ≠ "1"  -> scientific, no rounding
    @test format_number("1.23", "6", 2, 4) == "1.23 \\times 10^{6}"
    @test format_markdown_number("1.23", "6", 2, 4) == "1.23 × 10^6"
    # dig ≤ rdigits, expo == "1"  -> plain
    @test format_number("1.23", "1", 2, 4) == "1.23"
    @test format_markdown_number("1.23", "1", 2, 4) == "1.23"
    # dig > rdigits, expo ≠ "1"   -> rounded scientific
    @test format_number("1.23456", "6", 5, 2) == "1.23 \\times 10^{6}"
    @test format_markdown_number("1.23456", "6", 5, 2) == "1.23 × 10^6"
    # dig > rdigits, expo == "1"  -> rounded plain
    @test format_number("1.23456", "1", 5, 2) == "1.23"
    @test format_markdown_number("1.23456", "1", 5, 2) == "1.23"
end
