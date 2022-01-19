# Tests the GeoUnits
using Test
using GeoParams
#using Parameters

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
@test nondimensionalize(10cm/yr,CharUnits_GEO) ≈ 0.0031688087814028945   rtol=1e-10
@test nondimensionalize(10cm/yr,CharUnits_SI) ≈ 3.168808781402895e7   rtol=1e-10

@test nondimensionalize([10 1; 20 30.0],CharUnits_GEO) == [10 1; 20 30.0]

A = 10MPa*s^(-1)   # should give 1e12 in ND units
@test nondimensionalize(A,CharUnits_GEO)  ≈ 1e12 

B = GeoUnit(10MPa*s^(-1))   
@test nondimensionalize(B,CharUnits_GEO)  ≈ 1e12 

A = 10Pa^1.2*s^(-1)   
@test nondimensionalize(A,CharUnits_GEO)  ≈ 39810.71705534975 

CharUnits   =   GEO_units(viscosity=1e23, length=100km, stress=0.00371MPa);
@test nondimensionalize(A,CharUnits)  ≈ 1.40403327e16 

A = (1.58*10^(-25))*Pa^(-4.2)*s^(-1)        # calcite
@test nondimensionalize(A,CharUnits_GEO)  ≈ 3.968780561785161e16

R=8.314u"J/mol/K"
@test nondimensionalize(R,CharUnits_SI)  ≈ 8.314e-7

# test Dimensionalize in case we provide a number and units
v_ND      =   nondimensionalize(3cm/yr, CharUnits_GEO); 
@test Dimensionalize(v_ND, cm/yr, CharUnits_GEO) == 3.0cm/yr

# Test the GeoUnit struct
x=GeoUnit(8.1cm/yr)
@test  x.val==8.1
x_ND  = nondimensionalize(x,CharUnits_GEO)
@test  x_ND ≈ 0.002566735112936345 rtol=1e-8        # nondimensionalize a single value
x_dim = Dimensionalize(x,CharUnits_GEO)
@test  x_dim.val == 8.1                              # Dimensionalize again

y   =   x+2cm/yr
@test  y.val==10.1                            # errors

xx  =GeoUnit([8.1cm/yr; 10cm/yr]);
@test  xx.val==[8.1; 10]
yy=xx/(1cm/yr);                                         # transfer to no-unt


z=GeoUnit([100.0km 1000km 11km; 10km 2km 1km]);       # array
@test z/1km == [100.0 1000.0 11.0; 10.0 2.0 1.0]    # The division by 1km transfer it to a GeoUnit structure with no units; the multiplying with a float creates a float array

zz = GeoUnit((1:10)km)
@test zz/1km == 1:10


# Test non-dimensionalisation if z is an array
z_ND = nondimensionalize(z,CharUnits_GEO)
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
T_vec   =   273K:10K:500K        # using units
T       =   400.0K               # Unitful quanity
T_nd    =   10:10:200            # no units  
α       =   GeoUnit(3e-5/K)
T₀      =   GeoUnit(293K)
ρ₀      =   GeoUnit(3300kg/m^3)

T_geo   =   GeoUnit(400K)
T_f = 400

# Calculations with two GeoUnits
t = T_geo+T₀
@test isa(t, GeoUnit)    # addition
@test Unit(t) == K
@test t ≈ 693.0

t = T_geo-T₀
@test isa(t, GeoUnit)    # subtraction
@test t ≈ 107
@test Unit(t) == K

t = T_geo*T₀
@test isa(t, GeoUnit)    # multiplication
@test t ≈ 117200.0
@test Unit(t) == K^2

t = T_geo/T₀
@test isa(t, GeoUnit)    # division
@test t ≈ 1.3651877133105803
@test Unit(t) == NoUnits

# Case in which one of them is a Unitful quantity 
t = T+T₀
@test isa(t, Quantity)    # addition
@test t==693K

t = T-T₀
@test isa(t, Quantity)    # subtraction
@test t ≈ 107K

t = T*T₀
@test isa(t, Quantity)    # multiplication
@test t ≈ 117200.0K^2

t = T/T₀
@test isa(t, Float64)    # division
@test t ≈ 1.3651877133105803

# case in which one of them is a Float
t = T_f+T₀
@test isa(t, Float64)    # addition
@test t==693

t = T_f-T₀
@test isa(t, Float64)    # subtraction
@test t ≈ 107

t = T_f*T₀
@test isa(t, Float64)    # multiplication
@test t ≈ 117200.0

t = T_f/T₀
@test isa(t, Float64)    # division
@test t ≈ 1.3651877133105803


# case in which one of them is an array with Floats
t = T_nd+T₀
@test isa(t, AbstractArray)    # addition
@test t==303.0:10.0:493.0

t = T_nd-T₀
@test isa(t, AbstractArray)    # subtraction
@test t == -283.0:10.0:-93.0

t = T_nd*T₀
@test isa(t, AbstractArray)    # multiplication
@test t == 2930.0:2930.0:58600.0

t = T_nd/T₀
@test isa(t, AbstractArray)    # division
@test t ≈ 0.034129692832764506:0.034129692832764506:0.6825938566552902

# case in which one of them is an array with Quantities
t = T_vec+T₀
@test isa(t, AbstractArray)    # addition
@test t==(566.0:10.0:786.0)K

t = T_vec-T₀
@test isa(t, AbstractArray)    # subtraction
@test t == (-20.0:10.0:200.0)K

t = T_vec*T₀
@test isa(t, AbstractArray)    # multiplication
@test t == (79989.0:2930.0:144449.0)*K^2

t = T_vec/T₀
@test isa(t, AbstractArray)    # division
@test t ≈ 0.931740614334471:0.034129692832764506:1.68259385665529


ρ = ρ₀*(1.0 .- α*(T-T₀))  # The second expression is turned into a GeoUnit, so ρ should be GeoUnits
@test ρ ≈ 3289.4069999999997

ρ = ρ₀*(1.0 .- α*(T_nd-T₀))  # The second expression is turned into a GeoUnit, so ρ should be GeoUnits
@test ρ ≈ 3328.017:-0.99:3309.207


# test conversion:
b = convert(GeoUnit{Float64}, 3300kg/m^3)
@test b.val == 3300.0

c = convert(GeoUnit{Float32}, 3300.0kg/m^3)
@test c.val == 3300.0f0

# test broadcasting with float and unit arrays
ar   = [1.0 2;  4 5]            # Float array
ar1  = [1.0km 2km;  4km 5km]    # Unit array
test = GeoUnit(10km)

@test ar .+ test == [11.0 12.0; 14.0 15.0]
@test test .- ar == [9.0 8.0; 6.0 5.0]
@test test .* ar1 == [10.0km^2 20.0km^2; 40.0km^2 50.0km^2]
@test ar ./ ar1 == [1/km 1/km ; 1/km 1/km ]


# The way we define different structures heavily affects the number of allocs 
# that are done while computing with parameters defined within the struct 
struct ConDensity   # 1 allocation
    ρ::Float64
end

# Define struct with types and dimensions
struct ConDensity1{T} <: AbstractMaterialParam  # 1 allocation
    test::String
    ρ::GeoUnit{T}
    v::GeoUnit{T}
end

struct ConDensity2 <: AbstractMaterialParam # 3001 allocation
    test::String
    ρ::GeoUnit 
    v::GeoUnit
end

# No type info
mutable struct ConDensity3 <: AbstractMaterialParam  # 3001 allocation
    test::String
    ρ::GeoUnit 
end

# add type info in name, but no default values
mutable struct ConDensity4{T <: AbstractFloat} <: AbstractMaterialParam  # 1 allocation
    test::String
    ρ::GeoUnit{T}
    v::GeoUnit{T}
end

# use keywords, but no info about type
Base.@kwdef struct ConDensity5 <: AbstractMaterialParam # 3001 allocation
    test::String =""
    ρ::GeoUnit = GeoUnit(3300.1kg/m^3)
    v::GeoUnit = GeoUnit(100/s)
end

# This is slow
Base.@kwdef struct ConDensity6 <: AbstractMaterialParam # 3001 allocation
    test::String = ""
    ρ::GeoUnit = GeoUnit(3300.1kg/m^3)
    v::GeoUnit = GeoUnit(100/s)
end

# Here, we indicate the type:
Base.@kwdef struct ConDensity7{T} <: AbstractMaterialParam # 1 allocation
    test::String =""
    ρ::GeoUnit{T} = GeoUnit(3300.1kg/m^3)
    v::GeoUnit{T}   = GeoUnit(100/s)
end

# This seems to work, but has 3001 allocations
Base.@kwdef struct ConDensity8 <: AbstractMaterialParam # 3001 allocation
    test::String =""
    ρ::GeoUnit = 3300.1kg/m^3
    v::GeoUnit   = 100/s
end

# This works with 1 allocation, and is the preferred way to use it within GeoParams
Base.@kwdef struct ConDensity9{T <: AbstractFloat} <: AbstractMaterialParam # 1 allocation
    test::String =""
    ρ::GeoUnit{T} = 3300.1kg/m^3
    v::GeoUnit{T} = 100.0m/s
end
ConDensity9(a...) = ConDensity9{Float64}(a...)  # in case we do not give a type

rho  = ConDensity(3300.1)
rho1 = ConDensity1("t",GeoUnit(3300.1kg/m^3), GeoUnit(100/s))
rho2 = ConDensity2("t",3300.1kg/m^3, GeoUnit(100km/s))
rho3 = ConDensity3("t",3300.1kg/m^3)
rho4 = ConDensity4("t",GeoUnit(3300.1kg/m^3), GeoUnit(100/s))
rho5 = ConDensity5(test="t")
rho6 = ConDensity6("t",GeoUnit(3300.1kg/m^3), GeoUnit(100/s))
rho7 = ConDensity7()
rho8 = ConDensity8()
rho9 = ConDensity9(ρ=2800kg/m^3)
rho9_ND = nondimensionalize(rho9, GEO_units())


# test automatic nondimensionalization of a MaterialsParam struct:
CD = GEO_units()
rho2_ND = nondimensionalize(rho2, CD)
@test rho2_ND.ρ ≈ 3.3000999999999995e-18
@test rho2_ND.v ≈ 9.999999999999999e11

b = 20.1
c = 10.1;
r = 0.0

# Simple function to test speed
function f!(r,x,y)
    for i=1:1000
        r += x.ρ * y          # compute
    end
    r 
end

function g!(r,x,y)
    for i=1:1000
        r += x.Density[1].ρ * y          # compute
    end
    r 
end

#Phase1 = SetMaterialParams(Name="test1", Phase=22, Density  = ConDensity9())

#=
# testing speed (# of allocs)
r = 0.0
@btime f!($r, $rho,  $c)    # 1 allocations
@btime f!($r, $rho1, $c)    # 1 allocation
@btime f!($r, $rho2, $c)    # 3001 allocations
@btime f!($r, $rho3, $c)    # 3001 allocations
@btime f!($r, $rho4, $c)    # 1 allocation
@btime f!($r, $rho5, $c)    # 3001 allocation
@btime f!($r, $rho6, $c)    # 3001 allocation (so also with keywords, it is crucial to indicate the type)
@btime f!($r, $rho7, $c)    # 1 allocation (shows that we need to encode the units)
@btime f!($r, $rho8, $c)    # 3001 allocation 
@btime f!($r, $rho9, $c)    # 1 allocation (shows that we need to encode the units)
@btime f!($r, $rho9_ND, $c) 
=#

end