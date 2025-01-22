using Test
using GeoParams
using Unidecode
import GeoParams: Dislocation, Diffusion

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

    # test Phase2Dict()
    dict, ref = Phase2Dict(MatParam)
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

    # test Phase2DictMd()
    dictMd = Phase2DictMd(MatParam)
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

    ParameterTable(MatParam; format = "MD")
    @test "ParameterTable.md" in readdir()
    rm("ParameterTable.md"; force = true)

    filename = "TestTable"
    ParameterTable(MatParam; format = "MaRkDoWn", filename = "TestTable")
    @test "TestTable.md" in readdir()
    rm("TestTable.md"; force = true)

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

    # test Phase2Dict() for CompositeRheology
    dict, ref = Phase2Dict(MatParam)
    dval = MatParam[2].CompositeRheology[1][3].η.val
    @test dict["η CompositeRheology LinVisc 2.3"][1] == "$dval"
    @test dict["η CompositeRheology LinVisc 2.3"][2] == "\\" * "$(unidecode("η"))"
    @test dict["η CompositeRheology LinVisc 2.3"][3] ==
        "CompoRheo(DiffCreep,DislCreep,LinVisc,)"
    @test dict["η CompositeRheology LinVisc 2.3"][4] == "2"
    @test dict["η CompositeRheology LinVisc 2.3"][5] == "1"
    @test dict["η CompositeRheology LinVisc 2.3"][6] == "LinVisc"

    # test Phase2DictMd() for CompositeRheology
    dictMd = Phase2DictMd(MatParam)
    dvalMd = MatParam[2].CompositeRheology[1][3].η.val
    @test dictMd["η CompositeRheology LinVisc 2.3"][1] == "$dvalMd"
    @test dictMd["η CompositeRheology LinVisc 2.3"][2] == "η"
    @test dictMd["η CompositeRheology LinVisc 2.3"][3] ==
        "CompoRheo(DiffCreep,DislCreep,LinVisc,)"
    @test dictMd["η CompositeRheology LinVisc 2.3"][4] == "2"
    @test dictMd["η CompositeRheology LinVisc 2.3"][5] == "1"
    @test dictMd["η CompositeRheology LinVisc 2.3"][6] == "LinVisc"
end
