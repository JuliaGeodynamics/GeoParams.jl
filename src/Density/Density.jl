module Density

# This implements different methods to compute density
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, Unitful, LaTeXStrings
using ..Units
using ..PhaseDiagrams
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
using ..MaterialParameters: No_MaterialParam, MaterialParamsInfo
import Base.show, GeoParams.param_info

include("../Utils.jl")
include("../Computations.jl")

abstract type AbstractDensity{T} <: AbstractMaterialParam end

export  compute_density,        # calculation routines
        compute_density!,       # in place calculation
        param_info,             # info about the parameters
        AbstractDensity,        
        ConstantDensity,        # constant
        PT_Density,             # P & T dependent density
        Compressible_Density    # Compressible density 
       
# Define "empty" computational routines in case nothing is defined
compute_density!(rho::_T,s::No_MaterialParam{_T}; P::_T=zero(_T),T::_T=zero(_T)) where {_T} = zero(_T)
compute_density(s::No_MaterialParam{_T}; P::_T=zero(_T),T::_T=zero(_T)) where {_T} = zero(_T)

# Constant Density -------------------------------------------------------
"""
    ConstantDensity(ρ=2900kg/m^3)
    
Set a constant density:
```math  
    \\rho  = cst
```
where ``\\rho`` is the density [``kg/m^3``].
"""
@with_kw_noshow struct ConstantDensity{_T,U} <: AbstractDensity{_T}
    ρ::GeoUnit{_T,U} = 2900.0kg/m^3 # density
end
ConstantDensity(args...) = ConstantDensity(convert.(GeoUnit,args)...) 

(ρ::ConstantDensity)(; args...) =  ρ.ρ.val
(ρ::ConstantDensity)(args) =  ρ(; args...)
compute_density(s::ConstantDensity{_T}, args) where _T = s(; args...)
compute_density(s::ConstantDensity{_T}) where _T = s()

# This assumes that density always has a single parameter. If that is not the case, we will have to extend this (to be done)
function param_info(s::ConstantDensity) # info about the struct
    return MaterialParamsInfo(Equation = L"\rho = cst")
end

# Calculation routines
function compute_density!(rho::AbstractArray{_T}, s::ConstantDensity{_T}; kwargs...) where _T
    @unpack_val ρ   = s
    rho[:] .= ρ
    return nothing
end

compute_density!(rho::AbstractArray{_T}, s::ConstantDensity{_T}, args) where _T = compute_density!(rho, s; args...) 

# Print info 
function show(io::IO, g::ConstantDensity)  
    print(io, "Constant density: ρ=$(UnitValue(g.ρ))")  
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
    ρ0 ::GeoUnit{_T,U1}      =   2900.0kg/m^3                # density
    α  ::GeoUnit{_T,U2}      =   3e-5/K                      # T-dependence of density
    β  ::GeoUnit{_T,U3}      =   1e-9/Pa                     # P-dependence of density
    T0 ::GeoUnit{_T,U4}      =   0.0C                        # Reference temperature
    P0 ::GeoUnit{_T,U5}      =   0.0MPa                      # Reference pressure
end
PT_Density(args...) = PT_Density(convert.(GeoUnit,args)...) 

function param_info(s::PT_Density)  # info
    return MaterialParamsInfo(Equation = L"\rho = \rho_0(1.0-\alpha (T-T_0) + \beta (P-P_0)" )
end

# Calculation routine in case units are provided
function (ρ::PT_Density)(; P::Number, T::Number, kwargs...)
    @unpack_val ρ0, α, β, P0, T0 = ρ
    return ρ0*(1.0 - α*(T - T0) + β*(P - P0) )
end

(ρ::PT_Density)(args) = ρ(; args...)

compute_density(s::PT_Density{_T}, args) where _T = s(args)
compute_density(s::PT_Density{_T}, P::AbstractArray, T::AbstractArray) where _T = s(P=P, T=T)

function compute_density!(ρ::AbstractArray, s::PT_Density{_T}; P::_T, T::_T, kwargs...) where _T
    @unpack ρ0,α,β,P0, T0   = s
    
    ρ .= ρ0*(1.0 - α*(T-T0) + β*(P-P0) )

    return nothing
end

compute_density!(ρ::AbstractArray, s::PT_Density{_T}, args) where _T = compute_density!(ρ, s; args...) 

# Print info 
function show(io::IO, g::PT_Density)  
    print(io, "P/T-dependent density: ρ0=$(UnitValue(g.ρ0)), α=$(UnitValue(g.α)), β=$(UnitValue(g.β)), T0=$(UnitValue(g.T0)), P0=$(UnitValue(g.P0))")  
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
    ρ0::GeoUnit{_T,U1}     =   2900.0kg/m^3                # density
    β ::GeoUnit{_T,U2}     =   1e-9/Pa                     # P-dependence of density
    P0::GeoUnit{_T,U3}     =   0.0MPa                      # Reference pressure
end
Compressible_Density(args...) = Compressible_Density(convert.(GeoUnit,args)...) 

function param_info(s::Compressible_Density) # info about the struct
    return MaterialParamsInfo(Equation = L"\rho = \rho_0\exp(\beta*(P-P_0))"     )
end

function (s::Compressible_Density{_T})(; P::_T=zero(_T), kwargs...) where _T
    @unpack_val ρ0, β, P0   = s
    return ρ0*exp(β*(P - P0) )
end

(s::Compressible_Density{_T})(args) where _T = s(; args...)
compute_density(s::Compressible_Density{_T}, args) where _T = s(; args...)

function compute_density!(ρ::_T, s::Compressible_Density{_T}; P::_T, kwargs...) where _T
    # function compute_density!(ρ::_T, s::Compressible_Density{_T}, P::_T=zero(_T),T::_T=zero(_T)) where _T
    @unpack ρ0,β,P0   = s

    return ρ0*exp( β*(P-P0) )
end

compute_density!(ρ::_T, s::Compressible_Density{_T}, P::_T, kwargs...) where _T = compute_density!(ρ, s; P, kwargs)

function compute_density!(ρ::AbstractArray, s::Compressible_Density{_T}; P::_T, kwargs...) where _T
    @unpack ρ0,β,P0   = s
    Threads.@threads for i in eachindex(P)
        @inbounds ρ[i] = ρ0*exp(β*(P[i]-P0))
    end
    return nothing
end

compute_density!(ρ::AbstractArray, s::Compressible_Density{_T}, args) where _T = compute_density!(ρ, s; args...)

# Print info 
function show(io::IO, g::Compressible_Density)  
    print(io, "Compressible density: ρ0=$(UnitValue(g.ρ0)), β=$(UnitValue(g.β)), P0=$(UnitValue(g.P0))")  
end
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# Phase diagrams
function param_info(s::PhaseDiagram_LookupTable) # info about the struct
    return MaterialParamsInfo(Equation = L"\rho = f_{PhaseDiagram}(T,P))" )
end

"""
    compute_density(P,T, s::PhaseDiagram_LookupTable)
Interpolates density as a function of `T,P` from a lookup table  
"""
function (s::PhaseDiagram_LookupTable)(; P, T, kwargs...)
    fn = s.Rho
    return fn(T,P)
end
(s::PhaseDiagram_LookupTable)(args) = s(; args...)
compute_density(s::PhaseDiagram_LookupTable, args) = s(; args...)
compute_density(s::PhaseDiagram_LookupTable; P, T) = s(; P=P, T=T)

"""
    compute_density!(rho::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::PhaseDiagram_LookupTable)
In-place computation of density as a function of `T,P`, in case we are using a lookup table.    
"""
function compute_density!(rho::AbstractArray{_T}, s::PhaseDiagram_LookupTable; P::AbstractArray{_T}=[zero(_T)],T::AbstractArray{_T}=[zero(_T)], kwargs...) where _T
    rho[:] = s.Rho.(T,P)
    return nothing
end

compute_density!(rho::AbstractArray, s::PhaseDiagram_LookupTable, args) = compute_density!(rho, s, args...)

#------------------------------------------------------------------------------------------------------------------#
# Computational routines needed for computations with the MaterialParams structure 

# This assumes that density always has a single parameter. If that is not the case, we will have to extend this (to be done)
# function compute_density(s::AbstractMaterialParamsStruct, args...) where {_T}
#     return compute_density(s.Density[1], args...)
# end
function compute_density(s::AbstractMaterialParamsStruct, args)
    return s.Density[1](args)
end

# these routines may come handy when we have >1 field for Density
#function compute_density(s::NTuple{N,AbstractDensity{_T}}, P::_T=zero(_T), T::_T=zero(_T)) where {N,_T}
#    compute_density.(s,P,T)
#end

#now with Tuple of Tuples
#function compute_density!(rho::Vector{NTuple{M,_T}}, P::Number, T::Number, MatParam::NTuple{N, NTuple{M, AbstractMaterialParamsStruct}}) where {N,M,_T}
#    rho .= map(x->compute_density(P,T,x), MatParam)
#end

#-------------------------------------------------------------------------------------------------------------

"""
    compute_density!(rho::AbstractArray{_T, ndim}, MatParam::NTuple{N,AbstractMaterialParamsStruct}, Phases::AbstractArray{_I, ndim}; P=nothing, T=nothing) where {ndim,N,_T,_I<:Integer}
    
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
julia> compute_density!(rho, MatParam, Phases, P, T)
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
julia> @btime compute_density!(\$rho, \$MatParam, \$Phases, P=\$P, T=\$T)
    203.468 μs (0 allocations: 0 bytes)
```
_____________________________________________________________________________________________________________________________   
    
    compute_density!(rho::AbstractArray{_T, N}, MatParam::NTuple{K,AbstractMaterialParamsStruct}, PhaseRatios::AbstractArray{_T, M}, P=nothing, T=nothing)

In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)
"""
compute_density!(args...) = compute_param!(compute_density, args...) #Multiple dispatch to rest of routines found in Computations.jl
compute_density(args...) = compute_param(compute_density, args...)

#=
function compute_density!(rho::AbstractArray{_T, ndim}, MatParam::NTuple{N,AbstractMaterialParamsStruct}, Phases::AbstractArray{_I, ndim}, P=nothing, T=nothing) where {ndim,N,_T,_I<:Integer}
    
    Phase_tup = ntuple(i->MatParam[i].Phase, Val(N))

    # check if we compute with P/T or use default values    
    isnothing(T) ? compute_T = false : compute_T = true
    isnothing(P) ? compute_P = false : compute_P = true
    
    Tval = zero(_T)
    Pval = zero(_T)

    @inbounds for I in eachindex(Phases)
        phase   = find_ind(Phase_tup, Phases[I])

        # Extract relevant value if requested 
        if compute_T; Tval = T[I]; end
        if compute_P; Pval = P[I]; end
       
        # this computes density for ALL phases, which is a bit of an overkill as we only need density for a single phase
        rho_tup = compute_density(MatParam, Pval, Tval)    
        rho[I]  = rho_tup[phase]
    end
end


function compute_density!(rho::AbstractArray{_T, N}, MatParam::NTuple{K,AbstractMaterialParamsStruct}, PhaseRatios::AbstractArray{_T, M}, P=nothing, T=nothing) where {_T<:AbstractFloat, N,M, K}
    if M!=(N+1)
        error("The PhaseRatios array should have one dimension more than the other arrays")
    end

    # check if we compute with P/T or use default values    
    isnothing(T) ? compute_T = false : compute_T = true
    isnothing(P) ? compute_P = false : compute_P = true

    Tval = zero(_T)
    Pval = zero(_T)
    if compute_P
        Rindex = CartesianIndices(P)
    else
        Rindex = CartesianIndices(T)
    end     
    
    @inbounds for I in Rindex
        frac    = view(PhaseRatios, Tuple(I)..., 1:K)    # fraction of each phase @ point I 
        
        # Extract relevant value if requested 
        if compute_T; Tval = T[I]; end
        if compute_P; Pval = P[I]; end

        # compute point-wise density:
        rho[I] = compute_density_times_frac(frac,MatParam, Pval, Tval) 
   
    end

    return 
end

# Multiplies density with the fraction of a phase
function compute_density_times_frac(PhaseRatios::AbstractArray{_T, 1}, MatParam::NTuple{N,AbstractMaterialParamsStruct}, P::_T, T::_T) where {_T, N}

    value_tup = compute_density(MatParam, P, T)
    
    # Sum & multiply density with fraction
    val = zero(_T)
    @inbounds for j=1:N             
        val += PhaseRatios[j]*value_tup[j]
    end

    return val 
end
=#

end


# OBSOLETE routines; I leave them in for now (commented), in case we want to reuse pieces of that @ a later stage.
# much of what we have above is more efficient (apart when using the Phase Diagram lookup tables)
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
# OBSOLETE?
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
# Many allocations
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
