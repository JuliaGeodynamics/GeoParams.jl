module HeatCapacity

# This implements different methods to specify heat capacity
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
import Base.show

abstract type AbstractHeatCapacity <: AbstractMaterialParam end

export  ComputeHeatCapacity,                # calculation routines
        ComputeHeatCapacity!,               # in-place routine
        ConstantHeatCapacity,               # constant
        T_HeatCapacity_Whittacker           # T-dependent heat capacity

# Constant Heat Capacity -------------------------------------------------------
"""
    ConstantHeatCapacity(cp=1050J/mol/kg)
    
Set a constant heat capacity:
```math  
    cp  = cst
```
where ``cp`` is the thermal heat capacity [``J/kg/K``].
"""
@with_kw_noshow mutable struct ConstantHeatCapacity <: AbstractHeatCapacity
    equation::LaTeXString   =   L"cp = cst"     
    cp::GeoUnit             =   1050J/kg/K                # heat capacity
end

# Calculation routine
function ComputeHeatCapacity(P, T, s::ConstantHeatCapacity)
    @unpack cp   = s
    if length(T)>1
        return Value(cp)*ones(size(T))
    else
        return Value(cp)
    end

end

# Calculation routine
function ComputeHeatCapacity!(P, T, s::ConstantHeatCapacity)
    @unpack cp   = s
    if length(T)>1
        return Value(cp)*ones(size(T))
    else
        return Value(cp)
    end

end

function ComputeHeatCapacity!(cp_array::AbstractArray{<:AbstractFloat,N},P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, s::ConstantHeatCapacity) where N
    @unpack cp   = s
    
    cp_array .= NumValue(cp)
    
    return nothing
end

# Print info 
function show(io::IO, g::ConstantHeatCapacity)  
    print(io, "Constant heat capacity: cp=$(g.cp.val)")  
end
#-------------------------------------------------------------------------


# Temperature dependent heat capacity -------------------------------
"""
    T_HeatCapacity_Whittacker()
    
Sets a temperature-dependent heat capacity following the parameterization of Whittacker et al. (2009), Nature:
```math  
    Cp = (a + b T - c/T^2)/m 
```

where ``Cp`` is the heat capacity [``J/kg/K``], and ``a,b,c`` are parameters that dependent on the temperature `T`:
- a = 199.50 J/mol/K    if T<= 846 K
- a = 199.50 J/mol/K    if T> 846 K
- b = 0.0857J/mol/K^2   if T<= 846 K
- b = 0.0323J/mol/K^2   if T> 846 K
- c = 5e6J/mol*K        if T<= 846 K
- c = 47.9e-6J/mol*K    if T> 846 K
- molmass =   0.22178kg/mol 

Note that this is slightly different than the equation in the manuscript, as Cp is in J/kg/K (rather than ``J/mol/K`` as in eq.3/4 of the paper)
"""
@with_kw_noshow mutable struct T_HeatCapacity_Whittacker <: AbstractHeatCapacity
    # Note: the resulting curve was visually compared with Fig. 2 of the paper
    equation::LaTeXString   =   L"cp = (a + b*T - c/T^2)/m"     
    a0::GeoUnit             =   199.5J/mol/K                # prefactor for low T       (T<= 846 K)
    a1::GeoUnit             =   229.32J/mol/K               # prefactor for high T      (T>  846 K)
    b0::GeoUnit             =   0.0857J/mol/K^2             # linear term for low T     (T<= 846 K)
    b1::GeoUnit             =   0.0323J/mol/K^2             # linear term for high T    (T>  846 K)
    c0::GeoUnit             =   5e6J/mol*K                  # quadratic term for low T  (T<= 846 K)
    c1::GeoUnit             =   47.9e-6J/mol*K              # quadratic term for high T (T>  846 K)
    molmass::GeoUnit        =   0.22178kg/mol               # average molar mass 
    Tcutoff::GeoUnit        =   846K                        # cutoff temperature
end

# Calculation routine
function ComputeHeatCapacity(P,T, s::T_HeatCapacity_Whittacker)
    @unpack a0,a1,b0,b1,c0,c1, molmass, Tcutoff   = s
    
    cp = zeros(size(T)).*Value(a0)./Value(molmass)

    for i in eachindex(T)
        if T[i] <= Value(Tcutoff)
            a,b,c = Value(a0),Value(b0),Value(c0)
        else
            a,b,c = Value(a1),Value(b1),Value(c1)
        end
       
        cp[i] = (a + b*T[i] - c/T[i]^2)/molmass 
    end

    return cp
end

"""
    ComputeHeatCapacity!(cp_array::AbstractArray{<:AbstractFloat,N}, T::AbstractArray{<:AbstractFloat,N}, s::T_HeatCapacity_Whittacker) where N
 
Computes T-dependent heat capacity in-place    
"""
function ComputeHeatCapacity!(cp_array::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, s::T_HeatCapacity_Whittacker) where N
    @unpack a0,a1,b0,b1,c0,c1, molmass, Tcutoff   = s
    a0 = NumValue(a0)
    a1 = NumValue(a1)
    b0 = NumValue(b0)
    b1 = NumValue(b1)
    c0 = NumValue(c0)
    c1 = NumValue(c1)
    molmass = NumValue(molmass)
    Tcutoff = NumValue(Tcutoff)

    ind             =   (T .<= Tcutoff)
    T_local         =   view(T, ind )
    cp_array[ind]   =   (a0 .+ b0*T_local - c0./T_local.^2)/molmass 

    ind             =   (T .> Tcutoff)
    T_local         =   view(T, ind )
    cp_array[ind]   =   (a1 .+ b1*T_local - c1./T_local.^2)/molmass 

    return nothing
end

"""
    ComputeHeatCapacity!(cp_array::AbstractArray{<:AbstractFloat,N},P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, s::T_HeatCapacity_Whittacker) where N

Computes T-dependent heat capacity in-place      
"""
function ComputeHeatCapacity!(cp_array::AbstractArray{<:AbstractFloat,N},P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, s::T_HeatCapacity_Whittacker) where N

    ComputeHeatCapacity!(cp_array,T, s)

    return nothing
end

# Print info 
function show(io::IO, g::T_HeatCapacity_Whittacker)  
    print(io, "T-dependent heat capacity: cp/$(g.molmass.val)=$(g.a0.val) + $(g.b0.val)*T - $(g.c0.val)/T^2 (for T<=$(g.Tcutoff.val)); ");
    print(io, " cp/$(g.molmass.val)=$(g.a1.val) + $(g.b1.val)*T - $(g.c1.val)/T^2 (for T>$(g.Tcutoff.val)) \n");
end
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Heat capacity from phase diagram

# to be implemented - see density implementation

#-------------------------------------------------------------------------


# Help info for the calculation routines
"""
    Cp = ComputeHeatCapacity(P, T, s:<AbstractHeatCapacity)

Returns the heat capacity `Cp` at any temperature `T` and pressure `P` using any of the heat capacity laws implemented.

Currently available:
- ConstantHeatCapacity
- T\\_HeatCapacity_Whittacker

# Example 
Using dimensional units
```julia
julia> T  = (250:100:1250)*K;
julia> cp = T_HeatCapacity_Whittacker()
julia> Cp = ComputeHeatCapacity(0,T,cp)
11-element Vector{Unitful.Quantity{Float64, 𝐋² 𝚯⁻¹ 𝐓⁻², Unitful.FreeUnits{(kg⁻¹, J, K⁻¹), 𝐋² 𝚯⁻¹ 𝐓⁻², nothing}}}:
  635.4269997294616 J kg⁻¹ K⁻¹
  850.7470171764261 J kg⁻¹ K⁻¹
  962.0959598489883 J kg⁻¹ K⁻¹
 1037.5420433770641 J kg⁻¹ K⁻¹
 1097.3517921966488 J kg⁻¹ K⁻¹
 1149.2745563671706 J kg⁻¹ K⁻¹
 1157.7915050948404 J kg⁻¹ K⁻¹
 1172.3554874197264 J kg⁻¹ K⁻¹
 1186.9194697445964 J kg⁻¹ K⁻¹
  1201.483452069455 J kg⁻¹ K⁻¹
 1216.0474343943067 J kg⁻¹ K⁻¹
```


"""
ComputeHeatCapacity()


"""
    Cp = ComputeHeatCapacity(T::Any, s::AbstractHeatCapacity)

Computes heat capacity if only temperature (and not pressure) is specified
"""
ComputeHeatCapacity(T::Any, s::AbstractHeatCapacity) =  ComputeHeatCapacity(0,T, s)

"""
    Cp = ComputeHeatCapacity(s::ConstantHeatCapacity)

Returns heat capacity if we are sure that we will only employ constant heat capacity in the simulation
"""
ComputeHeatCapacity(s::ConstantHeatCapacity) =  ComputeHeatCapacity(0,0, s)


"""
    ComputeHeatCapacity!(Cp::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of heat capacity `Cp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.

"""
function ComputeHeatCapacity!(Cp::AbstractArray{<:AbstractFloat, N}, Phases::AbstractArray{<:Integer, N}, P::AbstractArray{<:AbstractFloat, N},T::AbstractArray{<:AbstractFloat, N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct, 1}) where N

    for i = 1:length(MatParam)
        
        if !isnothing(MatParam[i].HeatCapacity)
            # Create views into arrays (so we don't have to allocate)
            ind = Phases .== MatParam[i].Phase;
            cp_local    =   view(Cp, ind )
            P_local     =   view(P , ind )
            T_local     =   view(T , ind )

            ComputeHeatCapacity!(cp_local, P_local, T_local, MatParam[i].HeatCapacity[1] ) 
        end
        
    end

end

"""
    ComputeHeatCapacity!(Cp::AbstractArray{<:AbstractFloat}, Phases::AbstractArray{<:Integer}, T::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of heat capacity `Cp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.

"""
function ComputeHeatCapacity!(Cp::AbstractArray{<:AbstractFloat, N}, Phases::AbstractArray{<:Integer, N}, T::AbstractArray{<:AbstractFloat, N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct, 1}) where N

    for i = 1:length(MatParam)
        if !isnothing(MatParam[i].HeatCapacity)
            # Create views into arrays (so we don't have to allocate)
            ind = Phases .== i;
            cp_local    =   view(Cp, ind )
            T_local     =   view(T , ind )

            ComputeHeatCapacity!(cp_local, T_local, MatParam[i].HeatCapacity[1] ) 
        end
        
    end

end


"""
    ComputeDensity!(Cp::AbstractArray{<:AbstractFloat,N}, PhaseRatios::AbstractArray{<:AbstractFloat, M}, P::AbstractArray{<:AbstractFloat,N},T::AbstractArray{<:AbstractFloat,N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct})

In-place computation of heat capacity `Cp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `PhaseRatio` of every point is specified as an Integer in the `PhaseRatios` array, which has one dimension more than the data arrays (and has a phase fraction between 0-1)

"""
function ComputeHeatCapacity!(Cp::AbstractArray{<:AbstractFloat, N}, PhaseRatios::AbstractArray{<:AbstractFloat, M}, P::AbstractArray{<:AbstractFloat, N},T::AbstractArray{<:AbstractFloat, N}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct, 1}) where {N,M}
    
    if M!=(N+1)
        error("The PhaseRatios array should have one dimension more than the other arrays")
    end

    Cp .= 0.0;
    for i = 1:length(MatParam)
        
        Cp_local   = zeros(size(Cp))
        Fraction    = selectdim(PhaseRatios,M,i);
        if (maximum(Fraction)>0.0) & (!isnothing(MatParam[i].HeatCapacity))

            ComputeHeatCapacity!(Cp_local, P, T, MatParam[i].HeatCapacity[1] ) 

            Cp .= Cp .+ Cp_local.*Fraction
        end
        
    end

end


end