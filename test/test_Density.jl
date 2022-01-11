using Test
using GeoParams
@testset "Density.jl" begin

# This tests the MaterialParameters structure
CharUnits_GEO   =   GEO_units(viscosity=1e19, length=1000km);
                
# Define a linear viscous creep law
x1      =   ConstantDensity(ρ=2900kg/m^3)
@test x1.ρ.val == 2900

x1 = Nondimensionalize(x1,CharUnits_GEO)
@test x1.ρ.val ≈ 2.9e-16

x2      =   PT_Density()
@test x2.α.val==3e-5
@test x2.ρ0.val==2900

x2 = Nondimensionalize(x2,CharUnits_GEO)
@test x2.T0.val≈0.21454659702313156

# Compute with density
@test compute_density(1.0,1.0, x2) ≈ 2.8419999999999996e-16
@test compute_density(1.0,1.0, x1) ≈ 2.9e-16

# Read Phase diagram interpolation object
fname   =   "./test_data/Peridotite.in"
PD_data =   PerpleX_LaMEM_Diagram(fname);
@test PD_data.meltFrac(1500,1e7) ≈ 0.24368492372485706
@test PD_data.Rho(1500,1e7) ≈ 3042.836820256982
@test PD_data.meltRho(1500,1e7) ≈ 2662.227167592414
@test PD_data.rockRho(1500,1e7) ≈ 3165.467673917775


@test compute_density(1e7, 1500, PD_data) ≈ 3042.836820256982

# Do the same but non-dimensionalize the result
CharDim  =  GEO_units();
PD_data1 =  PerpleX_LaMEM_Diagram(fname, CharDim=CharDim);

rho_ND   =  PD_data1.Rho(Nondimensionalize(1500K,CharDim),Nondimensionalize(1e8*Pa,CharDim)) 
Vp_ND    =  PD_data1.Vp(Nondimensionalize(1500K,CharDim),Nondimensionalize(1e8*Pa,CharDim)) 
Vs_ND    =  PD_data1.Vs(Nondimensionalize(1500K,CharDim),Nondimensionalize(1e8*Pa,CharDim)) 

# redimensionalize and check with value from original structure that did not use non-dimensionalization 
@test   ustrip(Dimensionalize(rho_ND,kg/m^3,CharDim)) ≈ PD_data.Rho(1500,1e8) 
@test   ustrip(Dimensionalize(Vp_ND, km/s,  CharDim)) ≈ PD_data.Vp(1500,1e8) 
@test   ustrip(Dimensionalize(Vs_ND, km/s,  CharDim)) ≈ PD_data.Vs(1500,1e8) 


# Test computation of density for the whole computational domain, using arrays 
MatParam    =   Array{MaterialParams, 1}(undef, 3);
MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=0,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PerpleX_LaMEM_Diagram("test_data/Peridotite.in"));

MatParam[2] =   SetMaterialParams(Name="Crust", Phase=1,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = ConstantDensity(ρ=2900kg/m^3));

MatParam[3] =   SetMaterialParams(Name="UpperCrust", Phase=2,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PT_Density());

# test computing material properties
Phases              = ones(Int64,400,400)*0;
Phases[:,20:end] .= 1
Phases[:,300:end] .= 2

#Phases .= 2;

rho     = zeros(size(Phases))
T       =  ones(size(Phases))
P       =  ones(size(Phases))*10

compute_density!(rho, Phases, P,T, MatParam)   # 82 allocations
@test sum(rho)/400^2 ≈ 2920.6148898225

# test computing material properties when we have PhaseRatios, instead of Phase numbers
PhaseRatio  = zeros(size(Phases)...,length(MatParam));
for i in CartesianIndices(Phases)
    iz = Phases[i]
    I = CartesianIndex(i,iz+1)
    PhaseRatio[I] = 1.0  
end

compute_density!(rho, PhaseRatio, P,T, MatParam)
@test sum(rho)/400^2 ≈ 2920.6148898225

# test speed of assigning density 
rho1 = ConstantDensity();

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

function h!(r::Float64,x,y::Float64)
    ρ = x.Density.ρ
    for i=1:1000
        r += ρ * y          # compute
    end
    r 
end

function h1!(r::Float64,x,y::Float64)
    ρ = x.Density.ρ0
    for i=1:1000
        r += ρ * y          # compute
    end
    r 
end

c = 10.1;
r = 0.0


# This is the way this was originally implemented 
Base.@kwdef mutable struct TestMat1
    Density = nothing 
end

Base.@kwdef struct TestMat2{T}
    Density =  Vector{AbstractDensity{T}}(undef, 1); 
end

Base.@kwdef struct TestMat3{T}
    Density::AbstractDensity{T}
end

# By explicitly defining Density as a tuple, we gain speed
Base.@kwdef struct TestMat6{V <: Tuple}
    Density::V = ()
end

Base.@kwdef struct TestMat9{V1<:Tuple, V2<:Tuple}
    Density::V1 = ()
    Conductity::V2 = ()
end

# This is how it is implemented (we have to use different tuples for type stability) 
Base.@kwdef struct TestMat9cc{Float64, V1<:Tuple, V2<:Tuple}  
    Density::V1 = ()
    Conductity::V2 = ()
end

test1 = TestMat1((ConstantDensity(),))
test2 = TestMat2{Float64}((ConstantDensity(),))
test3 = TestMat3{Float64}(ConstantDensity{Float64}())
test6 = TestMat6((ConstantDensity{Float64}(),))
test9 = TestMat9(Density = (ConstantDensity{Float64}(), ) )
#test9b = TestMat9(Density = (ConstantDensity{Float64}(), ConstantDensity{Float64}()) )

#test9 = TestMat9cc(Density = (ConstantDensity{Float64}()) )
#test9b = TestMat9k(Density = (ConstantDensity{Float64}(), ConstantDensity{Float64}()) )

#test9 = TestMat9cc(Density = [ConstantDensity{Float64}()] )


#=
using BenchmarkTools
@btime f!($r, $rho1, $c) # 1 allocation (as expected)
@btime f!($r, $(test1.Density[1]), $c)      # 1 allocation (but non-typical usage)
@btime g!($r, $test1, $c)       # typical usage: 6001 allocations
@btime g!($r, $test2, $c)       # typical usage: 6001 allocations
@btime g!($r, $test6, $c)       # typical usage: 1 allocations
@btime g!($r, $test9, $c)       # typical usage: 1 allocations [as implemented]
=#


## testing summing up the fractions in a pointwise manner

# Mimic the way we store material properties
N = 5
MatProp = Vector{TestMat9}(undef,N)
for i=1:N
    MatProp[i] = test9
end
Frac = ones(N)/N

#using StructArrays
#MatStructArr = StructArray(MatProp)


# test different ways to compute the density without allocating
den  = ConstantDensity()
den1 = PT_Density()


# Generate a 2D array with density types:
bb=Array{AbstractDensity{Float64},2}(undef,5,2)
[bb[i] = den for i in eachindex(bb)]
bb[end] = den1
bb[3]   = No_Density()

P,T     = 10.0,11.0
rho_bb  = zeros(size(bb))


# This statement has 0 allocations and compute the result for all densities simultaneously:
#
# Yet, it has the disadvantage 
#rho_bb .= compute_density.(P, T, bb)

#using BenchmarkTools
#@btime $rho_bb .= compute_density.($P, $T, $bb)
#  31.946 ns (0 allocations: 0 bytes)

# Next, lets try a vector of Vectors
den  = ConstantDensity()
den1 = PT_Density()
aa=Vector{Vector{AbstractDensity{Float64}}}(undef,5)
[aa[i] = [den,den1] for i=1:5]
aa[end] = [den1, No_Density()]

P,T     = 10.0,11.0

rho_aa=Vector{Vector{Float64}}(undef,5)
[rho_aa[i] = [0.,0.]  for i=1:5]

# For some reason this sometimes allocates and sometimes not, which puzzles me:
#compute_density!(rho_aa, P, T, aa)
#@btime compute_density!($rho_aa, $P, $T, $aa)
# 71.674 ns (0 allocations: 0 bytes)
# yet, if I restart and recompile this it gives:
#138.251 ns (5 allocations: 160 bytes)


#= no allocations:
julia> @code_warntype compute_density!(rho_aa, P, T, aa)
Variables
  #self#::Core.Const(GeoParams.MaterialParameters.Density.compute_density!)
  ρ::Vector{Vector{Float64}}
  P::Float64
  T::Float64
  s::Vector{Vector{AbstractDensity{Float64}}}
Body::Nothing
1 ─ %1 = Base.broadcasted(GeoParams.MaterialParameters.Density.compute_density!, ρ, P, T, s)::Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(compute_density!), Tuple{Vector{Vector{Float64}}, Float64, Float64, Vector{Vector{AbstractDensity{Float64}}}}}
│        Base.materialize!(ρ, %1)
└──      return nothing
=#

#= with 5 allocations: 
julia> @code_warntype compute_density!(rho_aa, P, T, aa)
Variables
  #self#::Core.Const(GeoParams.MaterialParameters.Density.compute_density!)
  ρ::Vector{Vector{Float64}}
  P::Float64
  T::Float64
  s::Vector{Vector{AbstractDensity{Float64}}}
Body::Nothing
1 ─ %1 = Base.broadcasted(GeoParams.MaterialParameters.Density.compute_density!, ρ, P, T, s)::Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(compute_density!), Tuple{Vector{Vector{Float64}}, Float64, Float64, Vector{Vector{AbstractDensity{Float64}}}}}
│        Base.materialize!(ρ, %1)
└──      return nothing
=#

# Next, lets try a vector of Tuples
den  = ConstantDensity()
den1 = PT_Density()
P,T     = 10.0,11.0

cc=Vector{NTuple{2,AbstractDensity{Float64}}}(undef,5)
#cc = @SVector [SA[den,den1] for i=1:5]
[cc[i] = (den,den1) for i=1:5]
cc[end] = (den1, No_Density())


rho_cc=Vector{NTuple{2,Float64}}(undef,5)
[rho_cc[i] = (0.,0.)  for i=1:5]


#=
cc=Vector{NTuple{1,AbstractDensity{Float64}}}(undef,5)
[cc[i] = (den,) for i=1:5]
cc[end] = (den1,)
rho_cc=Vector{NTuple{1,Float64}}(undef,5)
[rho_cc[i] = (0.,)  for i=1:5]
=#


# Static arrays:
using StaticArrays
cc = @SVector [SA[den,den1] for i=1:5]
rho_cc = @SVector [SA[0.,0.] for i=1:5]

compute_density!(rho_cc, P, T, cc)



#---------------------------------------------------------------------------------------------------------------#
#Structures involving AbstractDensity

#Using a vector
rho = zeros(100)
P,T = 9.,10.
s = [PT_Density() for i in 1:100]

#@btime compute_density!($rho, $P, $T, $s) -> No allocations (132.506 ns (0 allocations: 0 bytes))

#The same with a tuple 
s_tup = ntuple(x->PT_Density(), Val(100))
  
#@btime compute_density!($rho, $P, $T, $s_tup) -> Still no allocations (195.721 ns (0 allocations: 0 bytes))


#Now using vector of vectors and vector of tuples 
den = ConstantDensity()
den1 = PT_Density()
den2 = No_Density()
P,T = 1.,2.

#Vector
rho_ss = Vector{Vector{Float64}}(undef,15)
[rho_ss[i] = [0.,0.,0.]  for i=1:15]

ss = Vector{Vector{AbstractDensity{Float64}}}(undef,15)
[ss[i] = [den,den1,den] for i in 1:14]
ss[end] = [den2,den2,den2]

#@btime compute_density!($rho_ss, $P, $T, $ss) -> This sometimes allocates, (690.066 ns (16 allocations: 768 bytes))

#Tuple
rho_ss=Vector{Vector{Float64}}(undef,15)
[rho_ss[i] = [0.,0.,0.]  for i=1:15]

rho_tup = Vector{NTuple{3,Float64}}(undef,15)
[rho_tup[i] = (0.,0.,0.) for i = 1:15]

ss_vectup = Vector{NTuple{3,AbstractDensity{Float64}}}(undef,15)
[ss_vectup[i] = (den,den1,den) for i=1:14]
ss_vectup[end] = (den1, den2, den2)
ss_tup = Tuple(ss_vectup)

#@btime compute_density!($rho_ss, $P, $T, $ss_vectup)-> 902.500 ns (34 allocations: 1.03 KiB)
#@btime compute_density!($rho_tup, $P, $T, $ss_tup) -> 35.650 ns (0 allocations: 0 bytes) works!!


#---------------------------------------------------------------------------------------------------------#

#Structures involving MaterialParams
P,T = 1.,1.
rho = zeros(3)

MatParam = Vector{MaterialParams}(undef, 3)
MatParam[1] = SetMaterialParams(Name="Crust", Phase=1,
                            CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pas)),
                            Density   = ConstantDensity(ρ=2900kg/m^3))
MatParam[2] = SetMaterialParams(Name="Lower Crust", Phase=2,
                            CreepLaws= (PowerlawViscous(n=5.), LinearViscous(η=1e21Pas)),
                            Density  = PT_Density(ρ0=3000kg/m^3))
MatParam[3] = SetMaterialParams(Name="Lower Crust", Phase=3,
                            CreepLaws= (PowerlawViscous(n=1.), LinearViscous(η=1e21Pas)),
                            Density  = No_Density())

Mat_tup = Tuple(MatParam)

#@btime compute_density!($rho,$P,$T,$MatParam)
#@btime compute_density!($rho,$P,$T,$Mat_tup)

#now using Vector of Tuples 
rho_mm = Vector{Vector{Float64}}(undef, 10)
[rho_mm[i] = [.0,.0,.0] for i in 1:10]

rho_tup = Vector{NTuple{3,Float64}}(undef, 10)
[rho_tup[i] = (0.,0.,0.) for i in 1:10]

mm = Vector{NTuple{3,MaterialParams}}(undef, 10)
[mm[i] = Mat_tup for i = 1:10]
mm_tup = Tuple(mm)

#@btime compute_density!($rho_mm, $P, $T, $mm) -> 564.571 ns (21 allocations: 368 bytes)
#@btime compute_density!($rho_tup, $P, $T, $mm) -> 934.615 ns (25 allocations: 1.05 KiB)
#@btime compute_density!($rho_tup, $P, $T, $mm_tup) -> 25.502 ns (0 allocations: 0 bytes) !!







#=
#using AbstractTensors
den = PT_Density()
den1 = ConstantDensity()
den2 =  No_Density()
P,T = 1.,1.
rho = zeros(100)

tup = ntuple(x->PT_Density(), Val(100))
s = Values(tup)
#@btime compute_density!($rho,$P,$T,$s) -> 117.538 ns (0 allocations: 0 bytes) 

#now with Values of Values
rho_ss = Vector{Vector{Float64}}(undef,15)
[rho_ss[i] = [0.,0.,0.]  for i=1:15]

#= Can construct as:
        s_val = Values(den, den1, den)
        ss_val = Values(ntuple(x->s_val, Val(15)))
    
    Or even do:
        s_end = Values(den2, den2, den2)
        ss = Vector{Values{3,AbstractDensity{Float64}}}(undef,15)
        [ss[i] = s_val for i=1:14]
        ss[end] = s_end
        ss_tup = Tuple(ss)
        ss_val = Values(ss_tup)
=#

#@btime compute_density!($rho_ss,$P,$T,$ss_val) -> 1.390 μs (42 allocations: 4.59 KiB)



#using StaticArrays
den = PT_Density()
den1 = No_Density()
P,T = 1.,1.
rho = zeros(100)
s = SVector{100}([den for i in 1:100])
#@btime compute_density!($rho,$P,$T,$s) ->  116.975 ns (0 allocations: 0 bytes)

#now using SVector of SVector
rho_ss = Vector{Vector{Float64}}(undef,15)
[rho_ss[i] = [0.,0.,0.]  for i=1:15]
ss_vec = SVector{15}([SVector{3}([ss[i][j] for j in 1:3]) for i in 1:15]) #@SVector[@SVector[...] for i in 1:15]

#@btime compute_density!($rho_ss,$P,$T,$ss_vec) -> 1.310 μs (43 allocations: 3.44 KiB)



=#

end