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
        Compressible_Density,   # Compressible density 
        No_Density,             # nothing defined
        AbstractDensity,
        fill_tup,       # to be moved
        max_length


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
@with_kw_noshow struct ConstantDensity{_T,U}   <: AbstractDensity{_T}
    equation::LaTeXString   =   L"\rho = cst"     
    ρ::GeoUnit{_T,U}        =   2900.0kg/m^3                # density
end
ConstantDensity(eq,args...) = ConstantDensity(eq, convert.(GeoUnit,args)...) 

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

function compute_density!(rho::Number,P::Number,T::Number, s::ConstantDensity{_T}) where _T
    @unpack_val ρ   = s
    #rho = ρ
    return ρ
end

function compute_density!(rho::AbstractArray{_T},P::_T,T::_T, s::ConstantDensity{_T}) where _T
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
@with_kw_noshow struct PT_Density{_T,U1,U2,U3,U4,U5} <: AbstractDensity{_T}
    equation::LaTeXString    =   L"\rho = \rho_0(1.0-\alpha (T-T_0) + \beta (P-P_0)"     
    ρ0 ::GeoUnit{_T,U1}      =   2900.0kg/m^3                # density
    α  ::GeoUnit{_T,U2}      =   3e-5/K                      # T-dependence of density
    β  ::GeoUnit{_T,U3}      =   1e-9/Pa                     # P-dependence of density
    T0 ::GeoUnit{_T,U4}      =   0.0C                        # Reference temperature
    P0 ::GeoUnit{_T,U5}      =   0.0MPa                      # Reference pressure
end
PT_Density(eq,args...) = PT_Density(eq, convert.(GeoUnit,args)...) 

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
#-------------------------------------------------------------------------


# Pressure-dependent density -------------------------------
"""
    Compressible_Density(ρ0=2900kg/m^3, β=1e-9/Pa, P₀=0MPa)
    
Set a pressure-dependent density:
```math  
    \\rho  = \\rho_0 \\exp(β*(P - P\\_0))  
```
where ``\\rho_0`` is the density [``kg/m^3``] at reference pressure ``P_0`` and ``\\beta`` the pressure dependence.

"""
@with_kw_noshow struct Compressible_Density{_T,U1,U2,U3} <: AbstractDensity{_T}
    equation::LaTeXString  =   L"\rho = \rho_0\exp(\beta*(P-P_0))"     
    ρ0::GeoUnit{_T,U1}     =   2900.0kg/m^3                # density
    β ::GeoUnit{_T,U2}     =   1e-9/Pa                     # P-dependence of density
    P0::GeoUnit{_T,U3}     =   0.0MPa                      # Reference pressure
end
Compressible_Density(eq,args...) = Compressible_Density(eq, convert.(GeoUnit,args)...) 

function compute_density(P::Number,T::Number, s::Compressible_Density{_T}) where _T
    @unpack_val ρ0,β,P0   = s
    
    ρ = ρ0*exp(β*(P - P0) )

    return ρ
end

function compute_density!(ρ::Number, P::Number,T::Number, s::Compressible_Density{_T}) where _T
    @unpack ρ0,β,P0   = s

    return ρ0*exp( β*(P-P0) )
end

function compute_density!(ρ::AbstractArray, P::Number,T::Number, s::Compressible_Density{_T}) where _T
    @unpack ρ0,β,P0   = s
    
    ρ .= ρ0*exp( β*(P-P0) )

    return nothing
end

# Print info 
function show(io::IO, g::Compressible_Density)  
    print(io, "Compressible density: ρ0=$(g.ρ0.val), β=$(g.β.val), P0=$(g.P0.val)")  
end
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# Phase diagrams
"""
    compute_density(P,T, s::PhaseDiagram_LookupTable)

Interpolates density as a function of `T,P` from a lookup table  
"""
function compute_density(P,T, s::PhaseDiagram_LookupTable)
    return s.Rho(T,P)
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


#-----------------------------------------------------------------------------------------------------------------#

# AbstractDensity

# Vector
#function compute_density!(ρ::Vector{_T}, P::Number,T::Number, s::Vector{AbstractDensity{_T}}) where {_T}
#    ρ .= compute_density.(P,T,s)
#end


# Tuple (this method allocates)
#function compute_density!(ρ::Vector{_T}, P::Number, T::Number, s::NTuple{N, AbstractDensity{_T}}) where {N, _T}
#    ρ .= compute_density.(P,T,s) 
#    #ρ .=  map(x->compute_density(P,T,x), s)
#end


# Tuple of Tuples
#function compute_density(P::Number, T::Number, s::NTuple{N, AbstractDensity{_T}}) where {N, _T}
#    map(x->compute_density(P,T,x), s)
#end

#function compute_density(ρ::NTuple{N,_T}, P::Number, T::Number, s::NTuple{N, AbstractDensity{_T}}) where {N, _T}
#    compute_density!.(ρ,P,T,s)
#end

#function compute_density!(rho::Vector{NTuple{M,_T}}, P::Number, T::Number, s::NTuple{N, NTuple{M, AbstractDensity{_T}}}) where {N,M,_T}
#    #rho .= compute_density!.(rho,P,T,s)
#    rho .= map(x->compute_density(P,T,x),s)
#end


#------------------------------------------------------------------------------------------------------------------#
# using AbstractMaterialParamsStruct


function compute_density(P::Number, T::Number, s::AbstractMaterialParamsStruct)
    return compute_density(P,T,s.Density[1])
end

#function compute_density!(rho::Number, P::Number, T::Number, s::AbstractMaterialParamsStruct)
#    return compute_density!(rho,P,T,s.Density[1])
#end

# with Tuple
function compute_density!(rho::Vector{_T}, P::Number, T::Number, MatParam::NTuple{N,AbstractMaterialParamsStruct}) where {N,_T}
    rho .= map(x->compute_density(P,T,x), MatParam)
end

function compute_density(P::Number, T::Number, MatParam::NTuple{N,AbstractMaterialParamsStruct}) where {N}
    map(x->compute_density(P,T,x), MatParam)
end

#function compute_density(P::Number, T::Number, MatParam::NTuple{N,AbstractMaterialParamsStruct}) where N
#    map(x->compute_density(P,T,x), MatParam)
#end

#now with Tuple of Tuples
#function compute_density!(rho::Vector{NTuple{M,_T}}, P::Number, T::Number, MatParam::NTuple{N, NTuple{M, AbstractMaterialParamsStruct}}) where {N,M,_T}
#    rho .= map(x->compute_density(P,T,x), MatParam)
#end

#-------------------------------------------------------------------------------------------------------------


# Automatically fill tuples with No_Density given a length n
function fill_tup(v::NTuple{N,Tuple{Vararg{AbstractMaterialParam}}}, n) where N
    ntuple(i->ntuple(j-> j <= length(v[i]) ? v[i][j] : No_Density(), Val(n)),Val(N))
end

# Find max element in a tuple
function find_max(t::NTuple{N,T}) where {N,T}
    max = t[1]
    @inbounds for i in 2:N
        if t[i] > max
            max = t[i]
        end
    end
    max
end

# Find inner tuple of maximum length
function max_length(t::NTuple{N, Tuple}) where N
    find_max(ntuple(x->length(t[x]), Val(N)))
end


"""
    compute_density!(rho::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractAbstractMaterialParamsStructStruct})

In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.

# Example
```julia
julia> MatParam = (SetMaterialParams(Name="Mantle", Phase=1,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pa*s)),
                        Density   = PT_Density()
                        ),
                    SetMaterialParams(Name="Crust", Phase=2,
                        CreepLaws= (PowerlawViscous(), LinearViscous(η=1e23Pas)),
                        Density   = ConstantDensity(ρ=2900kg/m^3))
                        );

julia> Phases = ones(Int64,400,400);
julia> Phases[:,20:end] .= 2
julia> rho     = zeros(size(Phases))
julia> T       =  ones(size(Phases))
julia> P       =  ones(size(Phases))*10
julia> compute_density!(rho, Phases, P,T, MatParam)
julia> rho
400×400 Matrix{Float64}:
2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  …  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  …  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  …  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
    ⋮                                            ⋮                                         ⋱     ⋮                                       ⋮                            
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  …  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  …  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
 2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91  2899.91     2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0  2900.0
```
The routine is made to minimize allocations:
```julia
julia> using BenchmarkTools
julia> @btime compute_density!(\$rho, \$Phases, \$P, \$T, \$MatParam)
    1.358 ms (28 allocations: 1.27 MiB)
```
"""

#= Many allocations
function compute_density!(rho::AbstractArray{_T, ndim}, Phases::AbstractArray{<:Integer, ndim}, P::AbstractArray{_T, ndim},T::AbstractArray{_T, ndim}, MatParam::NTuple{N,AbstractMaterialParamsStruct}) where {ndim,N,_T}
    # This is the general way to implement this.
    #  Remarks:
    #       1) It allocates a lot when used with a phase diagram
    #       2) We still have to check if this works with the way ParallelStencil operates
    
    # The phases in MatParam may be ordered differently (or start with zero)
    
    Phase_vec = [MatParam[i].Phase for i=1:N]
    Phase_ind = copy(Phases);
    for i=1:N
        Phase_ind[Phases.==Phase_vec[i]] .= i
    end
    
    rho_local = zeros(_T,N)
    for i in eachindex(Phase_ind)
        phase   =   Phase_ind[i]
        
        # This does not allocate but computes density for all phases simultaneously. 
        #  That is a bit of an overkill, but alternative approaches appear to allocate  
        compute_density!(rho_local,P[i],T[i], MatParam)     
        rho[i]  =   rho_local[phase]
    end   
    return
end
=#

#2 allocations that come from creating Phase_ind
function compute_density!(rho::AbstractArray{_T, ndim}, Phases::AbstractArray{<:Integer, ndim}, P::AbstractArray{_T, ndim},T::AbstractArray{_T, ndim}, MatParam::NTuple{N,AbstractMaterialParamsStruct}) where {ndim,N,_T}
    
    Phase_tup = ntuple(i->MatParam[i].Phase, Val(N))
    Phase_ind = Array{Int64}(undef, size(Phases))

    map!(x->find_ind(Phase_tup,x), Phase_ind, Phases) #stores Phases indexes in Phase_ind

    @inbounds for i in eachindex(Phase_ind)
        phase = Phase_ind[i]
        rho_tup = compute_density(P[i], T[i], MatParam)
        rho[i] = rho_tup[phase]
    end
end

function find_ind(x::NTuple{N,Int64}, k::Int64) where N
    @inbounds for i in 1:N
        if x[i] == k
            return i
        end
    end
    return 0
end

# OBSOLETE?
#=
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
=#

"""
    compute_density!(rho::AbstractArray{_T, N}, PhaseRatios::AbstractArray{_T, M}, P::AbstractArray{_T, N},T::AbstractArray{_T, N}, MatParam::NTuple{K,AbstractMaterialParamsStruct}) where {_T, N,M, K}
   
In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)

"""
function compute_density!(rho::AbstractArray{_T, N}, PhaseRatios::AbstractArray{_T, M}, P::AbstractArray{_T, N},T::AbstractArray{_T, N}, MatParam::NTuple{K,AbstractMaterialParamsStruct}) where {_T, N,M, K}
    if M!=(N+1)
        error("The PhaseRatios array should have one dimension more than the other arrays")
    end

    @inbounds for i in CartesianIndices(P)
        
        # hmm, I'm sure this can be generalized somehow..
        if N==1
            frac   =   view(PhaseRatios, i[1], : )
        elseif N==2
            frac   =   view(PhaseRatios, i[1], i[2], : )
        elseif N==3
            frac   =   view(PhaseRatios, i[1], i[2], i[3],  : )
        end

        # compute point-wise density:
        rho[i] = compute_density(frac, P[i], T[i], MatParam ) 
    end

    return 
end


function compute_density!(rho_vec::AbstractArray{_T, 1}, PhaseRatios::AbstractArray{_T, 1}, P::_T,T::_T, MatParam::NTuple{N,AbstractMaterialParamsStruct}) where {_T, N}

    compute_density!(rho_vec, P, T, MatParam ) 
    
    # Sum & multiply density with fraction
    rho = zero(_T)
    @inbounds for j=1:N             # @inbounds to not allocate
        rho += PhaseRatios[j]*rho_vec[j]
    end

    return rho 
end

function compute_density(PhaseRatios::AbstractArray{_T, 1}, P::_T,T::_T, MatParam::NTuple{N,AbstractMaterialParamsStruct}) where {_T, N}

    dens_tup = compute_density(P,T,MatParam)
    
    # Sum & multiply density with fraction
    rho = zero(_T)
    @inbounds for j=1:N             # @inbounds to not allocate
        rho += PhaseRatios[j]*dens_tup[j]
    end

    return rho 
end


#=
# OLD function (which actually works better in case of phase diagrams as it can do computations at omce)
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
=#


end
