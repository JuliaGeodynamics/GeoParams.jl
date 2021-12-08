module MeltingParam

# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, PhaseDiagram_LookupTable
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
    ComputeMeltingParam!(ϕ, P,T, p::AbstractPhaseDiagramsStruct)

Inplace computation of melt fraction in case we use a phase diagram lookup table. The table should have the collum `:meltFrac` specified.
"""
function ComputeMeltingParam!(ϕ, P,T, p::PhaseDiagram_LookupTable)
    ϕ    .=   p.meltFrac.(T,P)

    return nothing
end



# Help info for the calculation routines
"""
    ϕ  = ComputeMeltingParam(P,T, s:<MeltingParam_Caricchi)

Returns the melt fraction `ϕ`, for a given melting parameterization as a function of `T`,`P`

"""
ComputeMeltingParam()

"""
     ComputeMeltingParam!(ϕ, P,T, s:<MeltingParam_Caricchi)

In-place routine that computes melt fraction `ϕ`, for a given melting parameterization as a function of `T`,`P`

"""
ComputeMeltingParam!()



end