module SeismicVelocity

# This implements different methods to compute seismic velocities 
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)


using Parameters, LaTeXStrings, Unitful
using ..Units
using ..PhaseDiagrams
#using ..MaterialParameters
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
import Base.show


abstract type AbstractSeismicVelocity <: AbstractMaterialParam end

export  ComputePwaveVelocity,  ComputeSwaveVelocity,    # calculation routines
        ComputePwaveVelocity!, ComputeSwaveVelocity!,   # in place calculation
        ConstantSeismicVelocity                         # constant

# Constant Velocity -------------------------------------------------------
"""
    ConstantSeismicVelocity(Vp=8.1 km/s, Vs=4.5km/s)
    
Set a constant seismic P and S-wave velocity:
```math  
    V_p = cst \\
    V_s = cst
```
where ``V_p, V_s`` are the P-wave and S-wave velocities [``km/s``].
"""
@with_kw_noshow mutable struct ConstantSeismicVelocity <: AbstractSeismicVelocity
    equation::LaTeXString   =   L"v_p = cst \\ v_s = cst"     
    Vp::GeoUnit             =  8.1km/s               # P-wave velocity
    Vs::GeoUnit             =  4.5km/s               # S-wave velocity
end

# Calculation routines
function ComputePwaveVelocity(P,T, s::ConstantSeismicVelocity)
    @unpack Vp   = s
    if length(T)>1
        return Value(Vp).*ones(size(T))
    else
        return NumValue(Vp)
    end
end

function ComputeSwaveVelocity(P,T, s::ConstantSeismicVelocity)
    @unpack Vs   = s
    if length(T)>1
        return Value(Vs).*ones(size(T))
    else
        return NumValue(Vs)
    end
end

function ComputePwaveVelocity!(Vp_array::AbstractArray{<:AbstractFloat,N},P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, s::ConstantSeismicVelocity) where N
    @unpack Vp   = s
    
    Vp_array .= NumValue(Vp)
    
    return nothing
end

function ComputeSwaveVelocity!(Vs_array::AbstractArray{<:AbstractFloat,N},P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, s::ConstantSeismicVelocity) where N
    @unpack Vs   = s
    
    Vs_array .= NumValue(Vs)
    
    return nothing
end

# Print info 
function show(io::IO, g::ConstantSeismicVelocity)  
    print(io, "Constant seismic velocity: Vp=$(Value(g.Vp)), Vs=$(Value(g.Vs))")  
end
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# Phase diagrams
"""
    ComputePwaveVelocity(P,T, s::PhaseDiagram_LookupTable)

Interpolates `Vp` as a function of `T,P`   
"""
function ComputePwaveVelocity(P,T, s::PhaseDiagram_LookupTable)
    return s.Vp.(T,P)
end

"""
    ComputeSwaveVelocity(P,T, s::PhaseDiagram_LookupTable)

Interpolates `Vs` as a function of `T,P`   
"""
function ComputeSwaveVelocity(P,T, s::PhaseDiagram_LookupTable)
    return s.Vs.(T,P)
end

"""
    ComputePwaveVelocity!(Vp_array::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::PhaseDiagram_LookupTable)

In-place computation of P-wave velocity as a function of `T,P`, in case we are using a lookup table.    
"""
function ComputePwaveVelocity!(Vp_array::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::PhaseDiagram_LookupTable)
    Vp_array[:] = s.Vp.(T,P)
    return nothing
end

"""
    ComputeSwaveVelocity!(Vs_array::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::PhaseDiagram_LookupTable)

In-place computation of S-wave velocity as a function of `T,P`, in case we are using a lookup table.    
"""
function ComputeSwaveVelocity!(Vs_array::AbstractArray{<:AbstractFloat}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, s::PhaseDiagram_LookupTable)
    Vs_array[:] = s.Vs.(T,P)
    return nothing
end

#-------------------------------------------------------------------------

"""
    ComputePwaveVelocity!(Vp_array::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of P-wave velocity `Vp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.

"""
function ComputePwaveVelocity!(Vp_array::AbstractArray{<:AbstractFloat, N}, Phases::AbstractArray{<:Integer, N}, P::AbstractArray{<:AbstractFloat, N},T::AbstractArray{<:AbstractFloat, N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct, 1}) where N

    for i = 1:length(MatParam)

        if !isnothing(MatParam[i].SeismicVelocity)
            # Create views into arrays (so we don't have to allocate)
            ind = Phases .== MatParam[i].Phase
            Vp_local    =   view(Vp_array, ind )
            P_local     =   view(P  , ind )
            T_local     =   view(T  , ind )

            ComputePwaveVelocity!(Vp_local, P_local, T_local, MatParam[i].SeismicVelocity[1] ) 
        end
        
    end

end

"""
    ComputeSwaveVelocity!(Vs_array::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of S-wave velocity `Vp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.

"""
function ComputeSwaveVelocity!(Vs_array::AbstractArray{<:AbstractFloat, N}, Phases::AbstractArray{<:Integer, N}, P::AbstractArray{<:AbstractFloat, N},T::AbstractArray{<:AbstractFloat, N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct, 1}) where N

    for i = 1:length(MatParam)

        if !isnothing(MatParam[i].SeismicVelocity)
            # Create views into arrays (so we don't have to allocate)
            ind = Phases .== MatParam[i].Phase
            Vs_local    =   view(Vs_array, ind )
            P_local     =   view(P  , ind )
            T_local     =   view(T  , ind )

            ComputeSwaveVelocity!(Vs_local, P_local, T_local, MatParam[i].SeismicVelocity[1] ) 
        end
        
    end

end

"""
    ComputePwaveVelocity!(Vp_array::AbstractArray{<:AbstractFloat,N}, PhaseRatios::AbstractArray{<:AbstractFloat, M}, P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of seismic P-wave velocity `Vp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)

"""
function ComputePwaveVelocity!(Vp_array::AbstractArray{<:AbstractFloat, N}, PhaseRatios::AbstractArray{<:AbstractFloat, M}, P::AbstractArray{<:AbstractFloat, N},T::AbstractArray{<:AbstractFloat, N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct, 1}) where {N,M}
    
    if M!=(N+1)
        error("The PhaseRatios array should have one dimension more than the other arrays")
    end

    Vp_array .= 0.0;
    for i = 1:length(MatParam)
        
        Vp_local    = zeros(size(Vp_array))
        Fraction    = selectdim(PhaseRatios,M,i);
        if (maximum(Fraction)>0.0) & (!isnothing(MatParam[i].SeismicVelocity))

            ComputePwaveVelocity!(Vp_local, P, T, MatParam[i].SeismicVelocity[1] ) 

            Vp_array .= Vp_array .+ Vp_local.*Fraction
        end
        
    end

end

"""
    ComputeSwaveVelocity!(Vs_array::AbstractArray{<:AbstractFloat,N}, PhaseRatios::AbstractArray{<:AbstractFloat, M}, P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of seismic S-wave velocity `Vs` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)

"""
function ComputeSwaveVelocity!(Vs_array::AbstractArray{<:AbstractFloat, N}, PhaseRatios::AbstractArray{<:AbstractFloat, M}, P::AbstractArray{<:AbstractFloat, N},T::AbstractArray{<:AbstractFloat, N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct, 1}) where {N,M}
    
    if M!=(N+1)
        error("The PhaseRatios array should have one dimension more than the other arrays")
    end

    Vs_array .= 0.0;
    for i = 1:length(MatParam)
        
        Vs_local    = zeros(size(Vs_array))
        Fraction    = selectdim(PhaseRatios,M,i);
        if (maximum(Fraction)>0.0) & (!isnothing(MatParam[i].SeismicVelocity))

            ComputeSwaveVelocity!(Vs_local, P, T, MatParam[i].SeismicVelocity[1] ) 

            Vs_array .= Vs_array .+ Vs_local.*Fraction
        end
        
    end

end


end