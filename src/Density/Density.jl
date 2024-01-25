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

include("../Computations.jl")

abstract type AbstractDensity{T} <: AbstractMaterialParam end

export compute_density,     # calculation routines
    compute_density!,       # in place calculation
    compute_density_ratio, 
    param_info,             # info about the parameters
    AbstractDensity,
    ConstantDensity,        # constant
    PT_Density,             # P & T dependent density
    Compressible_Density    # Compressible density 

# Define "empty" computational routines in case nothing is defined
function compute_density!(
    rho::_T, s::No_MaterialParam{_T}; P::_T=zero(_T), T::_T=zero(_T)
) where {_T}
    return zero(_T)
end
function compute_density(s::No_MaterialParam{_T}; P::_T=zero(_T), T::_T=zero(_T)) where {_T}
    return zero(_T)
end

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
    ρ::GeoUnit{_T,U} = 2900.0kg / m^3 # density
end
ConstantDensity(args...) = ConstantDensity(convert.(GeoUnit, args)...)

@inline (ρ::ConstantDensity)(; args...)                          = ρ.ρ.val
@inline (ρ::ConstantDensity)(args)                               = ρ(; args...)
@inline compute_density(s::ConstantDensity{_T}, args) where {_T} = s(; args...)
@inline compute_density(s::ConstantDensity{_T})       where {_T} = s()

# This assumes that density always has a single parameter. If that is not the case, we will have to extend this (to be done)
function param_info(s::ConstantDensity) # info about the struct
    return MaterialParamsInfo(; Equation=L"\rho = cst")
end

# Calculation routines
function compute_density!(rho::AbstractArray, s::ConstantDensity; kwargs...)
    @unpack_val ρ = s
    rho[:] .= ρ
    return nothing
end

function compute_density!(rho::AbstractArray, s::ConstantDensity, args)
    return compute_density!(rho, s; args...)
end

# Print info 
function show(io::IO, g::ConstantDensity)
    return print(io, "Constant density: ρ=$(UnitValue(g.ρ))")
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
    ρ0::GeoUnit{_T,U1} = 2900.0kg / m^3 # density
    α::GeoUnit{_T,U2}  = 3e-5 / K       # T-dependence of density
    β::GeoUnit{_T,U3}  = 1e-9 / Pa      # P-dependence of density
    T0::GeoUnit{_T,U4} = 0.0C           # Reference temperature
    P0::GeoUnit{_T,U5} = 0.0MPa         # Reference pressure
end
PT_Density(args...) = PT_Density(convert.(GeoUnit, args)...)

function param_info(s::PT_Density) # info
    return MaterialParamsInfo(;
        Equation=L"\rho = \rho_0(1.0-\alpha (T-T_0) + \beta (P-P_0)"
    )
end

# Calculation routine 
@inline function (ρ::PT_Density)(; P::Number, T::Number, kwargs...)
    if T isa Quantity
        @unpack_units ρ0, α, β, P0, T0 = ρ
    else
        @unpack_val   ρ0, α, β, P0, T0 = ρ
    end

    # fma version of: ρ0 * (1.0 - α * (T - T0) + β * (P - P0))
    return ρ0 * fma(β, P - P0, fma(-α, T - T0, 1))
end

@inline (ρ::PT_Density)(args)                = ρ(; args...)
@inline compute_density(s::PT_Density, args) = s(args)
@inline function compute_density(s::PT_Density, P::AbstractArray, T::AbstractArray)
    s(; P=P, T=T)
end

# Print info 
function show(io::IO, g::PT_Density)
    return print(
        io,
        "P/T-dependent density: ρ0=$(UnitValue(g.ρ0)), α=$(UnitValue(g.α)), β=$(UnitValue(g.β)), T0=$(UnitValue(g.T0)), P0=$(UnitValue(g.P0))",
    )
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
    ρ0::GeoUnit{_T,U1} = 2900.0kg / m^3 # density
    β::GeoUnit{_T,U2}  = 1e-9 / Pa      # P-dependence of density
    P0::GeoUnit{_T,U3} = 0.0MPa         # Reference pressure
end
Compressible_Density(args...) = Compressible_Density(convert.(GeoUnit, args)...)

function param_info(s::Compressible_Density) # info about the struct
    return MaterialParamsInfo(; Equation = L"\rho = \rho_0\exp(\beta*(P-P_0))")
end

function (s::Compressible_Density{_T})(; P::_T=zero(_T), kwargs...) where {_T}
    if P isa Quantity
        @unpack_units ρ0, β, P0 = s
    else
        @unpack_val   ρ0, β, P0 = s
    end

    return ρ0 * exp(β * (P - P0))
end

@inline (s::Compressible_Density)(args)                = s(; args...)
@inline compute_density(s::Compressible_Density, args) = s(; args...)

# Print info 
function show(io::IO, g::Compressible_Density)
    return print(
        io,
        "Compressible density: ρ0=$(UnitValue(g.ρ0)), β=$(UnitValue(g.β)), P0=$(UnitValue(g.P0))",
    )
end
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Phase diagrams
function param_info(s::PhaseDiagram_LookupTable) # info about the struct
    return MaterialParamsInfo(; Equation=L"\rho = f_{PhaseDiagram}(T,P))")
end

"""
    compute_density(P,T, s::PhaseDiagram_LookupTable)
Interpolates density as a function of `T,P` from a lookup table  
"""
@inline function compute_density(s::PhaseDiagram_LookupTable; P, T, kwargs...)
    fn = s.Rho
    return fn(T, P)
end
@inline compute_density(s::PhaseDiagram_LookupTable, args) = compute_density(s; args...)

"""
    compute_density!(rho::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::PhaseDiagram_LookupTable)
In-place computation of density as a function of `T,P`, in case we are using a lookup table.    
"""
# function compute_density!(rho::AbstractArray{_T}, s::PhaseDiagram_LookupTable; P::AbstractArray{_T}=[zero(_T)],T::AbstractArray{_T}=[zero(_T)], kwargs...) where _T end

#------------------------------------------------------------------------------------------------------------------#
# Computational routines needed for computations with the MaterialParams structure 
function compute_density(s::AbstractMaterialParamsStruct, args)
    return compute_density(s.Density[1], args)
end
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
julia> args = (P=P, T=T)
julia> compute_density!(rho, MatParam, Phases, args)
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
@inline compute_density!(args::Vararg{Any, N})      where N = compute_param!(compute_density, args...) #Multiple dispatch to rest of routines found in Computations.jl
@inline compute_density(args::Vararg{Any, N})       where N = compute_param(compute_density, args...)
@inline compute_density_ratio(args::Vararg{Any, N}) where N = compute_param_times_frac(compute_density, args...)

end
