module Density

# This implements different methods to compute density
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using ..PhaseDiagrams
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
import Base.show

abstract type AbstractDensity{T} <: AbstractMaterialParam end

export  compute_density,        # calculation routines
        compute_density!,       # in place calculation
        ConstantDensity,        # constant
        PT_Density,             # P & T dependent density
<<<<<<< HEAD
        No_Density,
        AbstractDensity


# Define default struct in case nothing is defined
struct No_Density{_T} <: AbstractDensity{_T} end
No_Density() = No_Density{Float64}();
compute_density!(rho::Number,P::Number,T::Number, s::No_Density{_T}) where {_T} = zero(_T)
compute_density(P::Number,T::Number, s::No_Density{_T}) where {_T} = zero(_T)

=======
        No_Density,             # nothing defined
        AbstractDensity 
>>>>>>> origin/boris_refactorGeoUnit


# Define default struct in case nothing is defined
struct No_Density{_T} <: AbstractDensity{_T} end
No_Density() = No_Density{Float64}();
compute_density!(rho::Number,P::Number,T::Number, s::No_Density{_T}) where {_T} = zero(_T)
compute_density(P::Number,T::Number, s::No_Density{_T}) where {_T} = zero(_T)


# Constant Density -------------------------------------------------------
"""
    ConstantDensity(ρ=2900kg/m^3)
    
Set a constant density:
```math  
    \\rho  = cst
```
where ``\\rho`` is the density [``kg/m^3``].
"""
@with_kw_noshow struct ConstantDensity{_T} <: AbstractDensity{_T} 
    equation::LaTeXString   =   L"\rho = cst"     
    ρ::GeoUnit{_T}          =   2900.0kg/m^3                # density
end
ConstantDensity(a...) = ConstantDensity{Float64}(a...) 

# Calculation routines
function compute_density(P::Quantity,T::Quantity, s::ConstantDensity{_T}) where _T
    @unpack_units ρ   = s
    return ρ
end

function compute_density(P::Number,T::Number, s::ConstantDensity{_T}) where _T
    @unpack_val ρ   = s
    return ρ
end

function compute_density(P::AbstractArray,T::AbstractArray, s::ConstantDensity{_T}) where _T
    @unpack_val ρ   = s
    return ρ*ones(size(T))
end

<<<<<<< HEAD
#rho can't take new value
=======
>>>>>>> origin/boris_refactorGeoUnit
function compute_density!(rho::Number,P::Number,T::Number, s::ConstantDensity{_T}) where _T
    @unpack_val ρ   = s
    #rho = ρ
    return ρ
end

function compute_density!(rho::AbstractArray,P::Number,T::Number, s::ConstantDensity{_T}) where _T
    @unpack_val ρ   = s
    rho[:] .= ρ
    return nothing
end

function compute_density!(rho::AbstractArray,P::AbstractArray,T::AbstractArray, s::ConstantDensity{_T}) where _T
    @unpack_val ρ   = s
    rho .= ρ
    return nothing
end

# Print info 
function show(io::IO, g::ConstantDensity)  
    print(io, "Constant density: ρ=$(g.ρ.val)")  
end
#-------------------------------------------------------------------------

# Pressure & Temperature dependent density -------------------------------
"""
    PT_Density(ρ0=2900kg/m^3, α=3e-5/K, β=1e-9/Pa, T0=0C, P=0MPa)
    
Set a pressure and temperature-dependent density:
```math  
    \\rho  = \\rho_0 (1.0 - \\alpha (T-T_0) + \\beta  (P-P_0) )  
```
where ``\\rho_0`` is the density [``kg/m^3``] at reference temperature ``T_0`` and pressure ``P_0``,
``\\alpha`` is the temperature dependence of density and ``\\beta`` the pressure dependence.

"""
@with_kw_noshow struct PT_Density{_T} <: AbstractDensity{_T}
    equation::LaTeXString    =   L"\rho = \rho_0(1.0-\alpha (T-T_0) + \beta (P-P_0)"     
    ρ0::GeoUnit{_T}          =   2900.0kg/m^3                # density
    α::GeoUnit{_T}           =   3e-5/K                      # T-dependence of density
    β::GeoUnit{_T}           =   1e-9/Pa                     # P-dependence of density
    T0::GeoUnit{_T}          =   0.0C                        # Reference temperature
    P0::GeoUnit{_T}          =   0.0MPa                      # Reference pressure
end
PT_Density(a...) = PT_Density{Float64}(a...) 

# Calculation routine in case units are provided
function compute_density(P::Quantity,T::Quantity, s::PT_Density{_T}) where _T
    @unpack_units ρ0,α,β,P0, T0   = s
    
    ρ = ρ0*(1.0 - α*(T - T0) + β*(P - P0) )

    return ρ
end

# All other numbers 
function compute_density(P::Number,T::Number, s::PT_Density{_T}) where _T
    @unpack_val ρ0,α,β,P0, T0   = s
    
    ρ = ρ0*(1.0 - α*(T - T0) + β*(P - P0) )

    return ρ
end

# in-place calculation (ρ is immutable, its value won't change)
function compute_density!(ρ::Number, P::Number,T::Number, s::PT_Density{_T}) where _T
    @unpack ρ0,α,β,P0, T0   = s
    
    #ρ[:] = ρ0*(1.0 - α*(T-T0) + β*(P-P0) )

    return ρ0*(1.0 - α*(T-T0) + β*(P-P0) )
end


function compute_density!(ρ::AbstractArray, P::Number,T::Number, s::PT_Density{_T}) where _T
    @unpack ρ0,α,β,P0, T0   = s
    
    ρ .= ρ0*(1.0 - α*(T-T0) + β*(P-P0) )

    return nothing
end

function compute_density!(ρ::AbstractArray, P::Number,T::Number, s::PT_Density{_T}) where _T
    @unpack ρ0,α,β,P0, T0   = s
    
    ρ .= ρ0*(1.0 - α*(T-T0) + β*(P-P0) )

    return nothing
end

function compute_density!(ρ::AbstractArray,P::AbstractArray,T::AbstractArray, s::PT_Density{_T})  where _T
    @unpack_val ρ0,α,β,P0, T0   = s     # only values are required as we have floats as input

    # use the dot to avoid allocations
    @.  ρ  = ρ0*(1.0 - α*(T - T0) +  β*(P - P0))
    
    return nothing
end

# Print info 
function show(io::IO, g::PT_Density)  
    print(io, "P/T-dependent density: ρ0=$(g.ρ0.val), α=$(g.α.val), β=$(g.β.val), T0=$(g.T0.val), P0=$(g.P0.val)")  
end



#-----------------------------------------------------------------------------------------------------------------#

#AbstractDensity

#Vector
function compute_density!(ρ::Vector{_T}, P::Number,T::Number, s::Vector{AbstractDensity{_T}}) where {_T}
    ρ .= compute_density.(P,T,s)
end


#Tuple (this method allocates)
function compute_density!(ρ::Vector{_T}, P::Number, T::Number, s::NTuple{N, AbstractDensity{_T}}) where {N, _T}
    ρ .= compute_density.(P,T,s) 
    #ρ .=  map(x->compute_density(P,T,x), s)
end


#Tuple of Tuples
function compute_density(P::Number, T::Number, s::NTuple{N, AbstractDensity{_T}}) where {N, _T}
    #ρ .= compute_density.(P,T,s)
    map(x->compute_density(P,T,x), s)
end

function compute_density(ρ::NTuple{N,_T}, P::Number, T::Number, s::NTuple{N, AbstractDensity{_T}}) where {N, _T}
    compute_density!.(ρ,P,T,s)
end

function compute_density!(rho::Vector{NTuple{M,_T}}, P::Number, T::Number, s::NTuple{N, NTuple{M, AbstractDensity{_T}}}) where {N,M,_T}
    #rho .= compute_density!.(rho,P,T,s)
    rho .= map(x->compute_density(P,T,x),s)
end


#------------------------------------------------------------------------------------------------------------------#
#using AbstractMaterialParamsStruct


function compute_density(P::Number, T::Number, s::AbstractMaterialParamsStruct)
    return compute_density(P,T,s.Density[1])
end

function compute_density!(rho::Number, P::Number, T::Number, s::AbstractMaterialParamsStruct)
    return compute_density!(rho,P,T,s.Density[1])
end

#with Tuple
function compute_density!(rho::Vector{_T}, P::Number, T::Number, MatParam::NTuple{N,AbstractMaterialParamsStruct}) where {N,_T}
    rho .= map(x->compute_density(P,T,x), MatParam)
end


function compute_density(P::Number, T::Number, MatParam::NTuple{N,AbstractMaterialParamsStruct}) where N
    #compute_density.(P,T,MatParam)
    map(x->compute_density(P,T,x), MatParam)
end

#now with Tuple of Tuples
function compute_density!(rho::Vector{NTuple{M,_T}}, P::Number, T::Number, MatParam::NTuple{N, NTuple{M, AbstractMaterialParamsStruct}}) where {N,M,_T}
    rho .= map(x->compute_density(P,T,x), MatParam)
end


#-------------------------------------------------------------------------
# Phase diagrams
"""
    compute_density(P,T, s::PhaseDiagram_LookupTable)

Interpolates density as a function of `T,P`   
"""
function compute_density(P,T, s::PhaseDiagram_LookupTable)
    return s.Rho.(T,P)
end

"""
    compute_density!(rho::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::PhaseDiagram_LookupTable)

In-place computation of density as a function of `T,P`, in case we are using a lookup table.    
"""
function compute_density!(rho::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::PhaseDiagram_LookupTable)
    rho[:] = s.Rho.(T,P)
    return nothing
end

#-------------------------------------------------------------------------

"""
    compute_density!(rho::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractAbstractMaterialParamsStructStruct})

In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.

# Example
```julia
julia> MatParam    =   Array{MaterialParams, 1}(undef, 2);
julia> MatParam[1] =   SetMaterialParams(Name="Mantle", Phase=1,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PerpleX_LaMEM_Diagram("test_data/Peridotite.in"));
julia> MatParam[2] =   SetMaterialParams(Name="Crust", Phase=2,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pas)),
                        Density   = ConstantDensity(ρ=2900kg/m^3));
julia> Phases = ones(Int64,400,400);
julia> Phases[:,20:end] .= 2
julia> rho     = zeros(size(Phases))
julia> T       =  ones(size(Phases))
julia> P       =  ones(size(Phases))*10
julia> compute_density!(rho, Phases, P,T, MatParam)
julia> rho
400×400 Matrix{Float64}:
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46  …  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46  …  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
    ⋮                                            ⋮     ⋱                     ⋮                            
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46  …  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 3334.46  3334.46  3334.46  3334.46  3334.46  3334.46     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
```
The routine is made to minimize allocations:
```julia
julia> using BenchmarkTools
julia> @btime compute_density!(\$rho, \$Phases, \$P, \$T, \$MatParam)
 1.300 ms (44 allocations: 3.77 MiB)
```
"""
function compute_density!(rho::AbstractArray{<:AbstractFloat, N}, Phases::AbstractArray{<:Integer, N}, P::AbstractArray{<:AbstractFloat, N},T::AbstractArray{<:AbstractFloat, N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct, 1}) where {N,_T}

    for i = 1:length(MatParam)

        if !isempty(MatParam[i].Density)
            # Create views into arrays (so we don't have to allocate)
            ind = Phases .== MatParam[i].Phase
            rho_local   =   view(rho, ind )
            P_local     =   view(P  , ind )
            T_local     =   view(T  , ind )

           # density::Union{AbstractDensity{_T},PhaseDiagram_LookupTable}     =   MatParam[i].Density[1]
            density    =   MatParam[i].Density[1]
           
            #@time compute_density!(rho_local, P_local, T_local, dens ) 
            compute_density!(rho_local, P_local, T_local, density ) 
           
        end
        
    end

end



"""
    compute_density!(rho::AbstractArray{<:AbstractFloat,N}, PhaseRatios::AbstractArray{<:AbstractFloat, M}, P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)

"""
function compute_density!(rho::AbstractArray{<:AbstractFloat, N}, PhaseRatios::AbstractArray{<:AbstractFloat, M}, P::AbstractArray{<:AbstractFloat, N},T::AbstractArray{<:AbstractFloat, N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct, 1}) where {N,M}
    
    if M!=(N+1)
        error("The PhaseRatios array should have one dimension more than the other arrays")
    end
    
    rho1        =  zeros(size(rho))
    rho         .=  0.0;
    for i = 1:length(MatParam)

        Fraction    = selectdim(PhaseRatios,M,i);
        if (!isempty(MatParam[i].Density))

            ind         =   Fraction .> 0.0
            rho_local   =   view(rho1, ind )
            P_local     =   view(P   , ind )
            T_local     =   view(T   , ind )

            density = MatParam[i].Density[1]
            compute_density!(rho_local, P_local, T_local, density ) 

            rho[ind] .+= rho_local.*Fraction[ind]
        end
        
    end

end

<<<<<<< HEAD
end
=======

#function compute_density(P::Number,T::Number, s::Vector{AbstractDensity{_T}}) where {_T}
#    return compute_density.(P,T,s)
#end

#function compute_density!(ρ::Vector{_T}, P::Number,T::Number, s::Vector{AbstractDensity{_T}}) where {_T}
#    ρ .= compute_density.(P,T,s)
#end


@noinline function compute_density!(ρ::AbstractVector{_T},P::Number,T::Number, s::AbstractVector{AbstractDensity{_T}}) where {_T}
    @show s
    ρ .= compute_density.(P,T,s)
    #map!(x -> compute_density(P,T,x), ρ, s)

    return   
end

#multiple dis
function compute_density!(ρ::AbstractVector{AbstractVector{_T}},P::Number,T::Number, s::AbstractVector{AbstractVector{AbstractDensity{_T}}}) where {_T}

    #ρ = compute_density.(P,T,s)

#    map(x -> compute_density!(ρ, P,T,x), ρ,   s)

    return 
end

function compute_density!(P::Number,T::Number, s::NTuple{N,AbstractDensity{_T}}) where {N,_T}
   # @show N, s Base.Iterators.Flatten(s)
  #  ss = [0.0, 0.0];
  #  for (i, ss) in enumerate(s)
     #  [compute_density(P,T,ss) for ss in s]
        #ρ[i] = compute_density(P,T,ss)
        k = compute_density.(P,T,s)
        #@show ρ

      #  ρ = copy.(r)

   # end
 #   @show ss
    return sum(k)
end


# this works but allocates
function compute_density!(ρ::AbstractVector{NTuple{N,_T}},P::Number,T::Number, s::AbstractVector{NTuple{N,AbstractDensity{_T}}}) where {N,_T}
    #@show s

    map!(x->compute_density.(P,T,x), ρ, s)

    return 
end


end
>>>>>>> origin/boris_refactorGeoUnit
