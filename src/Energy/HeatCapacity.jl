module HeatCapacity

# This implements different methods to specify heat capacity
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
import Base.show, GeoParams.param_info
using ..MaterialParameters: MaterialParamsInfo

abstract type AbstractHeatCapacity{T} <: AbstractMaterialParam end

export compute_heatcapacity,               # calculation routines
    compute_heatcapacity!,              # in-place routine
    ConstantHeatCapacity,               # constant
    T_HeatCapacity_Whittington,          # T-dependent heat capacity
    param_info

#include("../Utils.jl")
include("../Computations.jl") 

# Constant Heat Capacity -------------------------------------------------------
"""
    ConstantHeatCapacity(cp=1050J/mol/kg)
    
Set a constant heat capacity:
```math  
    cp  = cst
```
where ``cp`` is the thermal heat capacity [``J/kg/K``].
"""
@with_kw_noshow struct ConstantHeatCapacity{T,U} <: AbstractHeatCapacity{T}
    cp::GeoUnit{T,U} = 1050J / kg / K                # heat capacity
end
ConstantHeatCapacity(args...) = ConstantHeatCapacity(convert.(GeoUnit, args)...)

function param_info(s::ConstantHeatCapacity) # info about the struct
    return MaterialParamsInfo(; Equation=L"c_p = cst")
end

# Calculation routine
function (s::ConstantHeatCapacity)(; kwargs...)
    @unpack_val cp = s
    return cp
end

compute_heatcapacity(s::ConstantHeatCapacity{_T}; kwargs...) where {_T} = s()

# Print info 
function show(io::IO, g::ConstantHeatCapacity)
    return print(io, "Constant heat capacity: cp=$(UnitValue(g.cp))")
end
#-------------------------------------------------------------------------

# Temperature dependent heat capacity -------------------------------
"""
    T_HeatCapacity_Whittington()
    
Sets a temperature-dependent heat capacity following the parameterization of Whittington et al. (2009), Nature:
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
@with_kw_noshow struct T_HeatCapacity_Whittington{T,U1,U2,U3,U4,U5} <:
                       AbstractHeatCapacity{T}
    # Note: the resulting curve was visually compared with Fig. 2 of the paper
    a0::GeoUnit{T,U1} = 199.5J / mol / K                # prefactor for low T       (T<= 846 K)
    a1::GeoUnit{T,U1} = 229.32J / mol / K               # prefactor for high T      (T>  846 K)
    b0::GeoUnit{T,U2} = 0.0857J / mol / K^2             # linear term for low T     (T<= 846 K)
    b1::GeoUnit{T,U2} = 0.0323J / mol / K^2             # linear term for high T    (T>  846 K)
    c0::GeoUnit{T,U3} = 5e6J / mol * K                  # quadratic term for low T  (T<= 846 K)
    c1::GeoUnit{T,U3} = 47.9e-6J / mol * K              # quadratic term for high T (T>  846 K)
    molmass::GeoUnit{T,U4} = 0.22178kg / mol               # average molar mass 
    Tcutoff::GeoUnit{T,U5} = 846K                        # cutoff temperature
end
T_HeatCapacity_Whittington(args...) = T_HeatCapacity_Whittington(convert.(GeoUnit, args)...)

function param_info(s::T_HeatCapacity_Whittington) # info about the struct
    return MaterialParamsInfo(; Equation=L"c_p = (a + b*T - c/T^2)/m")
end

# Calculation routine
function (s::T_HeatCapacity_Whittington{_T})(; T::_T=zero(_T), kwargs...) where {_T}
    @unpack_val a0, a1, b0, b1, c0, c1, molmass, Tcutoff = s

    cp = a0 / molmass

    if T <= Tcutoff
        a, b, c = a0, b0, c0
    else
        a, b, c = a1, b1, c1
    end

    cp = (a + b * T - c / T^2) / molmass

    return cp
end

"""
    compute_heatcapacity!(cp_array::AbstractArray{_T, N},s::T_HeatCapacity_Whittington{_T}, T::_T=zero(_T), P::_T=zero(_T)) where {_T,N}
 
Computes T-dependent heat capacity in-place    
"""
# function compute_heatcapacity!(cp_array::AbstractArray{_T, N},s::T_HeatCapacity_Whittington{_T}; T::AbstractArray{_T, N}, kwargs...) where {_T,N} end

# add methods programatically
for myType in (:ConstantHeatCapacity, :T_HeatCapacity_Whittington)
    @eval begin
        (s::$(myType))(args) = s(; args...)
        compute_heatcapacity(s::$(myType), args) = s(args)
    end
end

# Print info 
function show(io::IO, g::T_HeatCapacity_Whittington)
    return print(
        io,
        "T-dependent heat capacity following Whittington et al. (2009) for average crust). \n",
    )
end
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Heat capacity from phase diagram

# to be implemented - see density implementation

#-------------------------------------------------------------------------

# Help info for the calculation routines
"""
    Cp = compute_heatcapacity(s:<AbstractHeatCapacity, P, T)

Returns the heat capacity `Cp` at any temperature `T` and pressure `P` using any of the heat capacity laws implemented.

Currently available:
- ConstantHeatCapacity
- T\\_HeatCapacity_Whittington

# Example 
Using dimensional units
```julia
julia> T  = (250:100:1250)*K;
julia> cp = T_HeatCapacity_Whittington()
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
compute_heatcapacity()

"""
    Cp = ComputeHeatCapacity(T::Any, s::AbstractHeatCapacity)

Computes heat capacity if only temperature (and not pressure) is specified
"""
# compute_heatcapacity(s::AbstractHeatCapacity, T::AbstractArray{_T}) where _T =  compute_heatcapacity(s,similar(T), T)
# compute_heatcapacity!(cp_array::AbstractArray{_T}, s::AbstractHeatCapacity, T::AbstractArray{_T}) where _T =  compute_heatcapacity!(cp_array,s,similar(T), T)

"""
    Cp = ComputeHeatCapacity(s::ConstantHeatCapacity)

Returns heat capacity if we are sure that we will only employ constant heat capacity in the simulation
"""

# Computational routines needed for computations with the MaterialParams structure 
function compute_heatcapacity(s::AbstractMaterialParamsStruct, args)
    return s.HeatCapacity[1](args)
end

"""
    compute_heatcapacity!(Cp::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat})

In-place computation of heat capacity `Cp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.
"""
compute_heatcapacity!()


compute_heatcapacity(args...) = compute_param(compute_heatcapacity, args...)
compute_heatcapacity!(args...) = compute_param!(compute_heatcapacity, args...)

# In case just temperature is provided
function compute_heatcapacity!(
    Cp::AbstractArray{_T,ndim},
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    Phases::AbstractArray{_I,ndim},
    T::AbstractArray{_T,ndim},
) where {_T,ndim,N,_I<:Integer}
    return compute_param!(compute_heatcapacity, Cp, MatParam, Phases, nothing, T)
end

end
