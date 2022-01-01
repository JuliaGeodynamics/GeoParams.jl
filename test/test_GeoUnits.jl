# Tests the GeoUnits
using Test
using GeoParams

@testset "Units" begin

# test GeoUnits
a = GeoUnit(100km);
@test a.val==100
@test Unit(a)==km
@test isdimensional(a)==true

b = GeoUnit(100);
@test b.val==100
@test Unit(b)==NoUnits
@test isdimensional(b)==false

a = GeoUnit(100km:10km:200km);
@test a.val==100:10:200
@test Unit(a)==km
@test isdimensional(a)==true

b = GeoUnit(100:10:200);
@test b.val==100:10:200
@test Unit(b)==NoUnits
@test isdimensional(b)==false

a = GeoUnit([100km 200km; 300km 400km]);
@test a.val==[100 200; 300 400]
@test Unit(a)==km
@test isdimensional(a)==true

b = GeoUnit([100.0 200; 300 400]);
@test b.val==[100.0 200.0; 300.0 400.0]
@test Unit(b)==NoUnits
@test isdimensional(b)==false

# test creating structures
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
@test CharUnits_GEO.length  ==  1000km
@test CharUnits_GEO.Pa      ==  10000000Pa
@test CharUnits_GEO.Mass    ==  1.0e37kg
@test CharUnits_GEO.Time    ==  1.0e12s
@test CharUnits_GEO.Length  ==  1000000m


CharUnits_SI   =   SI_units();
@test CharUnits_SI.length   ==  1000m
@test CharUnits_SI.Pa       ==  10Pa

CharUnits_NO   =   NO_units();
@test CharUnits_NO.length   ==  1
@test CharUnits_NO.Pa       ==  1

# test nondimensionization of various parameters
@test Nondimensionalize(10cm/yr,CharUnits_GEO) ≈ 0.0031688087814028945   rtol=1e-10
@test Nondimensionalize(10cm/yr,CharUnits_SI) ≈ 3.168808781402895e7   rtol=1e-10

A = 10MPa*s^(-1)   # should give 1e12 in ND units
@test Nondimensionalize(A,CharUnits_GEO)  ≈ 1e12 

B = GeoUnit(10MPa*s^(-1))   
@test Nondimensionalize(B,CharUnits_GEO)  ≈ 1e12 

A = 10Pa^1.2*s^(-1)   
@test Nondimensionalize(A,CharUnits_GEO)  ≈ 39810.71705534975 

CharUnits   =   GEO_units(viscosity=1e23, length=100km, stress=0.00371MPa);
@test Nondimensionalize(A,CharUnits)  ≈ 1.40403327e16 

A = (1.58*10^(-25))*Pa^(-4.2)*s^(-1)        # calcite
@test Nondimensionalize(A,CharUnits_GEO)  ≈ 3.968780561785161e16

R=8.314u"J/mol/K"
@test Nondimensionalize(R,CharUnits_SI)  ≈ 8.314e-7

# test Dimensionalize in case we provide a number and units
v_ND      =   Nondimensionalize(3cm/yr, CharUnits_GEO); 
@test Dimensionalize(v_ND, cm/yr, CharUnits_GEO) == 3.0cm/yr

# Test the GeoUnit struct
x=GeoUnit(8.1cm/yr)
@test  x.val==8.1
x_ND  = Nondimensionalize(x,CharUnits_GEO)
@test  x_ND ≈ 0.002566735112936345 rtol=1e-8        # Nondimensionalize a single value
x_dim = Dimensionalize(x,CharUnits_GEO)
@test  x_dim.val == 8.1                              # Dimensionalize again

y   =   x+2cm/yr
@test  y.val==10.1                            # errors

xx  =GeoUnit([8.1cm/yr; 10cm/yr]);
@test  xx.val==[8.1; 10]
yy=xx/(1cm/yr);                                         # transfer to no-unt


z=GeoUnit([100.0km 1000km 11km; 10km 2km 1km]);       # array
@test z/1km*1.0 == [100.0 1000.0 11.0; 10.0 2.0 1.0]    # The division by 1km transfer it to a GeoUnit structure with no units; the multiplying with a float creates a float array

zz = GeoUnit((1:10)km)


# Test non-dimensionalisation if z is an array
z_ND = Nondimensionalize(z,CharUnits_GEO)
@test z_ND.val==[0.1   1.0    0.011; 0.01  0.002  0.001]
@test isdimensional(z_ND)==false          

z_D = Dimensionalize(z_ND,CharUnits_GEO)
@test z_D.val==[100 1000 11; 10 2 1]          # transform back
@test isdimensional(z_D)==true           

# test extracting a value from a GeoUnit array
@test NumValue(z_D[2,2]) == 2.0

# test setting a new value 
z_D[2,1]=3
@test z_D[2,1].val == 3.0

# Conversion to different units
@test convert(GeoUnit,Float64(10.1)).val==10.1
@test convert(GeoUnit,10.1:.1:20).val == 10.1:.1:20
@test Unit(convert(GeoUnit, 10km/s))==km/s
@test convert(Float64,  GeoUnit(10.2)) == 10.2
@test convert(Float64,  GeoUnit([10.2 11.2])) == [10.2 11.2]


# test various calculations (using arrays with and without units)
#T       =   273K:10K:500K        # using units
T       =   400.0K
T_nd    =   10:10:200            # no units  
α       =   GeoUnit(3e-5/K)
T₀      =   GeoUnit(293K)
ρ₀      =   GeoUnit(3300kg/m^3)

ρ = ρ₀*(1.0 .- α*(T-T₀))  # The second expression is turned into a GeoUnit, so ρ should be GeoUnits



# FEW DEBUGGING AND TESTS, mostly to determine the speed
using BenchmarkTools
function f(x,y)
    r=x*y
end
# with which we can access values
A = GeoUnit(10.1);
B = GeoUnit(12.1);

struct ConDensity
    ρ::Float64
end

struct ConDensity1
    test::String
    ρ::GeoUnit
end

struct ConDensity2 <: AbstractMaterialParam
    test::String
    ρ::GeoUnit 
end

mutable struct ConDensity4 <: AbstractMaterialParam
    test::String
    ρ::GeoUnit 
end

rho=ConDensity(3300.1)
rho1=ConDensity1("t",GeoUnit(3300.1kg/m^3))
rho2=ConDensity2("t",3300.1kg/m^3)
rho4=ConDensity4("t",3300.1kg/m^3)

b = 20.1
c = 10.1;
r = 0.0
#@show @btime f($b,$c)

#@btime f($rho.ρ,$b)

#64f(rho2.ρ, c)

end