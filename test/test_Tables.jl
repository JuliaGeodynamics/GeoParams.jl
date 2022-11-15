using Test
using GeoParams

@testset "Tables.jl" begin
MatParam = (SetMaterialParams(Name="Viscous Matrix", Phase=1, Density=ConstantDensity(),CreepLaws = SetDislocationCreep("Quartz Diorite | Hansen & Carter (1982)")),
            SetMaterialParams(Name="Viscous Sinker", Phase=2, Density= PT_Density(),CreepLaws = LinearViscous(η=1e21Pa*s)),
            SetMaterialParams(Name="Viscous Bottom", Phase=3, Density= PT_Density(),CreepLaws = SetDislocationCreep("Diabase | Caristan (1982)")))
n1 = "1"
n2 = "10.1234"
n3 = "1.12345e6"
n4 = "12.34567e-8"
dig, num, ex = detachFloatfromExponent(n1)
@test dig == 0
@test num == "1"
@test ex == "1"


end