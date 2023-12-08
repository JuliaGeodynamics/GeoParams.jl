using Test
using GeoParams
using Unidecode

### 1. Values für E_e und ..._melt parameter in Table einfügen
### 2. Noch die Strings im Dict anpassen z.B. so einer: Compo(Disl, Diff,...)
### 3. Test cases dafür anpassen
MatParam = (
               #Name="UpperCrust"
               SetMaterialParams(;
                   Phase   = 1,
                   Density  = PT_Density(ρ0=2700kg/m^3, β=1e10/Pa),
                   HeatCapacity = ConstantHeatCapacity(cp=1050J/kg/K),
                   Conductivity = ConstantConductivity(k=3.0Watt/K/m),
                   LatentHeat = ConstantLatentHeat(Q_L=350e3J/kg),
                   CompositeRheology = CompositeRheology((LinearViscous(; η=1e21 * Pa * s), SetConstantElasticity(; G=25e9Pa, ν=0.5), DruckerPrager_regularised(; C=10MPa, ϕ=30.0, η_vp=1.0e14Pas, Ψ=0), )),
                   Melting = MeltingParam_Caricchi(),
                   Elasticity = SetConstantElasticity(; G=25e9Pa, ν=0.5),
                   ),

               #Name="Magma"
               SetMaterialParams(;
                   Phase   = 2,
                   Density  = PT_Density(ρ0=2600kg/m^3, β=1e8/Pa),
                   HeatCapacity = ConstantHeatCapacity(cp=1050J/kg/K),
                   Conductivity = ConstantConductivity(k=1.5Watt/K/m),
                   LatentHeat = ConstantLatentHeat(Q_L=350e3J/kg),
                   CompositeRheology = CompositeRheology(Parallel((LinearViscous(), SetConstantElasticity(; G=10e9Pa, ν=0.3))), SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)"), SetDiffusionCreep("Dry Diopside | Dimanov & Dresen (2005)")),
                   Melting = MeltingParam_Caricchi(),
                   Elasticity = SetConstantElasticity(; G=10e9Pa, ν=0.3),
                   ),

               #Name="Thermal Anomaly"
               SetMaterialParams(;
                   Phase   = 3,
                   Density  = PT_Density(ρ0=2600kg/m^3, β=1e8/Pa),
                   HeatCapacity = ConstantHeatCapacity(cp=1050J/kg/K),
                   Conductivity = ConstantConductivity(k=1.5Watt/K/m),
                   LatentHeat = ConstantLatentHeat(Q_L=350e3J/kg),
                   CompositeRheology = Parallel(CompositeRheology((LinearViscous(; η=1e16 * Pa * s),SetConstantElasticity(; G=10e9Pa, ν=0.3))), LinearViscous(; η=1e15*Pa*s)),
                   Melting = MeltingParam_Caricchi(),
                   Elasticity = SetConstantElasticity(; G=10e9Pa, ν=0.3),
                   ),

               #Name="Sticky Air"
               SetMaterialParams(;
                   Phase   = 4,
                   Density   = ConstantDensity(ρ=100kg/m^3,),
                   HeatCapacity = ConstantHeatCapacity(cp=1000J/kg/K),
                   Conductivity = ConstantConductivity(k=15Watt/K/m),
                   LatentHeat = ConstantLatentHeat(Q_L=0.0J/kg),
                   CompositeRheology = CompositeRheology(Parallel((LinearViscous(; η=1e16*Pa*s),SetConstantElasticity(; ν=0.5, Kb=0.101MPa)),LinearViscous(;η=1e17*Pa*s))),
                   Elasticity = SetConstantElasticity(; ν=0.5, Kb=0.101MPa),
                   ),
                   )

@testset "Tables.jl" begin
    MatParam = (
        SetMaterialParams(;
            Name="Viscous Matrix",
            Phase=1,
            Density=ConstantDensity(),
            CreepLaws=SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)"),
        ),
        SetMaterialParams(;
            Name="Viscous Sinker",
            Phase=2,
            Density=PT_Density(),
            CreepLaws=LinearViscous(; η=1e21Pa * s),
        ),
        SetMaterialParams(;
            Name="Viscous Bottom",
            Phase=3,
            Density=PT_Density(),
            CreepLaws=SetDislocationCreep("Diabase | Caristan (1982)"),
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
    @test dict["Name 1"][1] == "$(join(MatParam[1].Name))"
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
    @test dictMd["Name 1"][1] == "$(join(MatParam[1].Name))"
    @test dictMd["Name 1"][2] == "$(length(MatParam))"
    @test dictMd["Name 1"][3] == "1"
    @test dictMd["Name 1"][4] == ""
    @test dictMd["Name 1"][5] == ""
    @test dictMd["Name 1"][6] == ""

    # test Dict2LatexTable()
    Dict2LatexTable(dict, ref)
    @test "ParameterTable.tex" in readdir()
    @test "References.bib" in readdir()

    rm("ParameterTable.tex"; force=true)
    rm("References.bib"; force=true)

    # test Dict2MarkdownTable()
    Dict2MarkdownTable(dictMd)
    @test "ParameterTable.md" in readdir()

    rm("ParameterTable.md"; force=true)

    # test ParameterTable()
    # for Latex
    ParameterTable(MatParam)
    @test "ParameterTable.tex" in readdir()
    @test "References.bib" in readdir()

    ParameterTable(MatParam; format="TEX")
    @test "ParameterTable.tex" in readdir()
    @test "References.bib" in readdir()
    rm("ParameterTable.tex"; force=true)
    rm("References.bib"; force=true)

    ParameterTable(MatParam; format="LaTeX", filename="TestTable")
    @test "TestTable.tex" in readdir()
    rm("TestTable.tex"; force=true)
    rm("References.bib"; force=true)

    # for Markdown
    ParameterTable(MatParam; format="Markdown")
    @test "ParameterTable.md" in readdir()

    ParameterTable(MatParam; format="MD")
    @test "ParameterTable.md" in readdir()
    rm("ParameterTable.md"; force=true)

    filename = "TestTable"
    ParameterTable(MatParam; format="MaRkDoWn", filename="TestTable")
    @test "TestTable.md" in readdir()
    rm("TestTable.md"; force=true)

    # test phase with CompositeRheology field
    v1 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
    c1 = CompositeRheology(
        v1,
        SetDislocationCreep("Diabase | Caristan (1982)"),
        LinearViscous(; η=1e21Pa * s),
        v1,
    )
    MatParam = (
        SetMaterialParams(;
            Name="Viscous Matrix",
            Phase=1,
            Density=ConstantDensity(),
            CreepLaws=SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)"),
        ),
        SetMaterialParams(;
            Name="Viscous Sinker", Phase=2, Density=PT_Density(), CompositeRheology=c1
        ),
        SetMaterialParams(;
            Name="Viscous Bottom",
            Phase=3,
            Density=PT_Density(),
            CreepLaws=SetDislocationCreep("Diabase | Caristan (1982)"),
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
