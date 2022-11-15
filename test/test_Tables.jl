using Test
using GeoParams
using Unidecode

@testset "Tables.jl" begin
MatParam = (SetMaterialParams(Name="Viscous Matrix", Phase=1, Density=ConstantDensity(),CreepLaws = SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)")),
            SetMaterialParams(Name="Viscous Sinker", Phase=2, Density= PT_Density(),CreepLaws = LinearViscous(η=1e21Pa*s)),
            SetMaterialParams(Name="Viscous Bottom", Phase=3, Density= PT_Density(),CreepLaws = SetDislocationCreep("Diabase | Caristan (1982)")))

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

#test Phase2Dict()
dict, ref = Phase2Dict(MatParam)
dval = MatParam[1].Density[1].ρ.val
@test dict["ρ Density 1"][1] == "$dval"
@test dict["ρ Density 1"][2] == "\\" * "$(unidecode("ρ"))"
@test dict["ρ Density 1"][3] == ""
@test dict["ρ Density 1"][4] == "1"
@test dict["Name 1"][1] == "$(join(MatParam[1].Name))"
@test dict["Name 1"][2] == "$(length(MatParam))"
@test dict["Name 1"][3] == "1"
@test dict["Name 1"][4] == ""

#test Phase2DictMd()
dictMd = Phase2DictMd(MatParam)
dvalMd = MatParam[1].Density[1].ρ.val
@test dictMd["ρ Density 1"][1] == "$dvalMd"
@test dictMd["ρ Density 1"][2] == "ρ"
@test dictMd["ρ Density 1"][3] == ""
@test dictMd["ρ Density 1"][4] == "1"
@test dictMd["Name 1"][1] == "$(join(MatParam[1].Name))"
@test dictMd["Name 1"][2] == "$(length(MatParam))"
@test dictMd["Name 1"][3] == "1"
@test dictMd["Name 1"][4] == ""

#test Dict2LatexTable()


end