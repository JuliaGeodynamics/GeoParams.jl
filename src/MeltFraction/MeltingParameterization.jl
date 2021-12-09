module MeltingParam

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, PhaseDiagram_LookupTable, AbstractMaterialParamsStruct
import Base.show

abstract type AbstractMeltingParam <: AbstractMaterialParam end

export  ComputeMeltingParam, ComputeMeltingParam!,   # calculation routines
        MeltingParam_Caricchi                        # constant
        

# Constant  -------------------------------------------------------
"""
    MeltingParam_Caricchi()
    
Implements the T-dependent melting parameterisation used by Caricchi et al 
```math  
    \\theta = (800.0 .- (T + 273.15))./23.0 
```
```math  
    \\phi_{solid} = 1.0 - {1.0 \\over (1.0 + e^\\theta)}; 
```

Note that T is in Kelvin.

"""
@with_kw_noshow mutable struct MeltingParam_Caricchi <: AbstractMeltingParam
    equation::LaTeXString   =   L"\phi = {1 \over 1 + \exp( {800-T[^oC] \over 23})}"     
    a::GeoUnit              =   800.0K              
    b::GeoUnit              =   23.0K
    c::GeoUnit              =   273.15K # shift from C to K
end

# Calculation routine
function ComputeMeltingParam(P,T, p::MeltingParam_Caricchi)
    @unpack a,b,c   = p

    θ       =   (Value(a) .- (T .- Value(c)))./Value(b)
    ϕ       =   1.0./(1.0 .+ exp.(θ))

    return ϕ
end

function ComputeMeltingParam!(ϕ, P,T, p::MeltingParam_Caricchi)
    @unpack a,b,c   = p

    θ       =   (Value(a) .- (T .- Value(c)))./Value(b)
    ϕ      .=   1.0./(1.0 .+ exp.(θ))

    return nothing
end

function ComputeMeltingParam!(ϕ::Array{<:AbstractFloat}, P::Array{<:AbstractFloat},T::Array{<:AbstractFloat}, p::MeltingParam_Caricchi)
    @unpack a,b,c   = p
    a = ustrip(Value(a))
    b = ustrip(Value(b))
    c = ustrip(Value(c))

    θ       =   (a .- (T .- c))./b
    ϕ      .=   1.0./(1.0 .+ exp.(θ))

    return nothing
end

function ComputeMeltingParam!(ϕ::SubArray{<:AbstractFloat}, P::SubArray{<:AbstractFloat},T::SubArray{<:AbstractFloat}, p::MeltingParam_Caricchi)
    @unpack a,b,c   = p
    a = ustrip(Value(a))
    b = ustrip(Value(b))
    c = ustrip(Value(c))

    θ       =   (a .- (T .- c))./b
    ϕ      .=   1.0./(1.0 .+ exp.(θ))

    return nothing
end

# Print info 
function show(io::IO, g::MeltingParam_Caricchi)  
    print(io, "Caricchi et al. melting parameterization")  
end
#-------------------------------------------------------------------------


"""
    ComputeMeltingParam(P,T, p::AbstractPhaseDiagramsStruct)

Computes melt fraction in case we use a phase diagram lookup table. The table should have the collum `:meltFrac` specified.
"""
function ComputeMeltingParam(P,T, p::PhaseDiagram_LookupTable)
   return p.meltFrac.(T,P)
end

"""
    ComputeMeltingParam!(ϕ::Array{Float64}, P::Array{Float64},T:Array{Float64}, p::PhaseDiagram_LookupTable)

In-place computation of melt fraction in case we use a phase diagram lookup table. The table should have the collum `:meltFrac` specified.
"""
function ComputeMeltingParam!(ϕ::Array{<:AbstractFloat}, P::Array{<:AbstractFloat},T::Array{<:AbstractFloat}, p::PhaseDiagram_LookupTable)
    ϕ[:]    =   p.meltFrac.(T,P)

    return nothing
end

"""
    ComputeMeltingParam!(ϕ::SubArray{Float64}, P::SubArray{Float64},T::SubArray{Float64}, p::PhaseDiagram_LookupTable)

In-place computation of melt fraction in case we use a phase diagram lookup table. The table should have the collum `:meltFrac` specified.
"""
function ComputeMeltingParam!(ϕ::SubArray{<:AbstractFloat}, P::SubArray{<:AbstractFloat},T::SubArray{<:AbstractFloat}, p::PhaseDiagram_LookupTable)
    ϕ[:]    =   p.meltFrac.(T,P)

    return nothing
end



"""
    ComputeMeltingParam!(ϕ::Array{Float64}, Phases::Array{Int64}, P::Array{Float64},T::Array{Float64}, MatParam::Array{<:AbstractMaterialParamsStruct})

In-place computation of density `rho` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
"""
function ComputeMeltingParam!(ϕ::Array{<:AbstractFloat, N}, Phases::Array{<:Integer, N}, P::Array{<:AbstractFloat, N},T::Array{<:AbstractFloat, N}, MatParam::Array{<:AbstractMaterialParamsStruct, 1}) where N


    for i = 1:length(MatParam)

        if !isnothing(MatParam[i].Melting)

            # Create views into arrays (so we don't have to allocate)
            ind = Phases .== i;
            ϕ_local   =   view(ϕ, ind )
            P_local   =   view(P  , ind )
            T_local   =   view(T  , ind )

            ComputeMeltingParam!(ϕ_local, P_local, T_local, MatParam[i].Melting[1] ) 
        end

    end

    return nothing
end





end