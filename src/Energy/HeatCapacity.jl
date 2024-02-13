module HeatCapacity

# This implements different methods to specify heat capacity
#
# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct, @extractors, add_extractor_functions
import Base.show, GeoParams.param_info
import ..Units: isdimensional
using ..MaterialParameters: MaterialParamsInfo

abstract type AbstractHeatCapacity{T} <: AbstractMaterialParam end

export compute_heatcapacity,               # calculation routines
    compute_heatcapacity!,               # in-place routine
    ConstantHeatCapacity,                # constant
    T_HeatCapacity_Whittington,          # T-dependent heat capacity
    Latent_HeatCapacity,                 # Implement latent heat by modifying heat capacity
    param_info

include("../Computations.jl")

# Constant Heat Capacity -------------------------------------------------------
"""
    ConstantHeatCapacity(Cp=1050J/mol/kg)

Set a constant heat capacity:
```math
    Cp  = cst
```
where ``Cp`` is the thermal heat capacity [``J/kg/K``].
"""
@with_kw_noshow struct ConstantHeatCapacity{T,U} <: AbstractHeatCapacity{T}
    Cp::GeoUnit{T,U} = 1050J / kg / K                # heat capacity
end
ConstantHeatCapacity(args...) = ConstantHeatCapacity(convert.(GeoUnit, args)...)

function param_info(s::ConstantHeatCapacity) # info about the struct
    return MaterialParamsInfo(; Equation=L"c_p = cst")
end

# Calculation routine
@inline function compute_heatcapacity(a::ConstantHeatCapacity;  kwargs...)
    @unpack_val Cp = a
    return Cp
end

# Print info
function show(io::IO, g::ConstantHeatCapacity)
    return print(io, "Constant heat capacity: Cp=$(UnitValue(g.Cp))")
end
#-------------------------------------------------------------------------

# Temperature dependent heat capacity -------------------------------
"""
    T_HeatCapacity_Whittington()

Sets a temperature-dependent heat capacity following the parameterization of Whittington et al. (2009), Nature:
```math
    Cp = (a + b T - c/T^2)/m
```

where ``Cp`` is the heat capacity [``J/kg/K``], and ``a,b,c`` are parameters that dependent on the temperature ``T`` [``K``]:
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
@inline function compute_heatcapacity(
    a::T_HeatCapacity_Whittington{_T};
    T = zero(_T),
    kwargs...
) where _T
    @unpack_val a0, a1, b0, b1, c0, c1, molmass, Tcutoff = a

    Cp = a0 / molmass

    if T <= Tcutoff
        a, b, c = a0, b0, c0
    else
        a, b, c = a1, b1, c1
    end

    Cp = (a + b * T - c / T^2) / molmass

    return Cp
end

# LatentHeat by modifying heat capacity  ---------------------------------
"""
    Latent_HeatCapacity(Cp=ConstantHeatCapacity(), Q_L=400kJ/kg)

This takes the effects of latent heat into account by modifying the heat capacity in the temperature equation:

```math
\\rho C_p^{\\textrm{new}} \\frac{\\partial T}{\\partial t}  = \\frac{\\partial }{\\partial x_i} \\left( k \\frac{\\partial T}{\\partial x_i} \\right)  + H_s
```

with
```math
C_p^{\\textrm{new}}  = C_p + \\frac{\\partial \\phi}{\\partial T} Q_L
```
where ``Q_L`` is the latent heat [``kJ/kg``], and ``\\frac{\\partial \\phi}{\\partial T}`` is the derivative of the melt fraction with respect to temperature

"""
@with_kw_noshow struct Latent_HeatCapacity{T,U,S1<:AbstractHeatCapacity} <: AbstractHeatCapacity{T}
    Cp::S1 = ConstantHeatCapacity()
    Q_L::GeoUnit{T,U} = 400kJ/kg                            # Latent heat
end
Latent_HeatCapacity(args...) = Latent_HeatCapacity(args[1], convert.(GeoUnit, args[2:end])...)
Latent_HeatCapacity(Cp::AbstractHeatCapacity, args...) = Latent_HeatCapacity(Cp, convert.(GeoUnit, args)...)
isdimensional(g::Latent_HeatCapacity) = isdimensional(g.Q_L)

function param_info(s::Latent_HeatCapacity) # info about the struct
    return MaterialParamsInfo(; Equation=L"Q_L = cst; c_p = Cp + \frac{\partial \phi}{\partial T}")
end

# Calculation routine
@inline function compute_heatcapacity(
    a::Latent_HeatCapacity{_T};
    dϕdT = zero(_T),
    kwargs...
) where _T
    @unpack_val Q_L = a

    Cp = compute_heatcapacity(a.Cp, kwargs)

    Cp += Q_L*dϕdT      # latent heat contribution

    return Cp
end

# Print info
function show(io::IO, g::Latent_HeatCapacity)
    return print(io, "Latent heat through modifying heat capacity: Cp = Cp + dϕdT*$(Value(g.Q_L)); Cp=$(g.Cp)")
end
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
"""
    compute_heatcapacity!(Cp_array::AbstractArray{_T, N},s::T_HeatCapacity_Whittington{_T}, T::_T=zero(_T), P::_T=zero(_T)) where {_T,N}

Computes T-dependent heat capacity in-place
"""
# function compute_heatcapacity!(Cp_array::AbstractArray{_T, N},s::T_HeatCapacity_Whittington{_T}; T::AbstractArray{_T, N}, kwargs...) where {_T,N} end

# add methods programmatically
for myType in (:ConstantHeatCapacity, :T_HeatCapacity_Whittington, :Latent_HeatCapacity)
    @eval begin
        #(s::$(myType))(args) = s(; args...)
        compute_heatcapacity(s::$(myType), args) = compute_heatcapacity(s; args...)
    end
end

# Print info
function show(io::IO, g::T_HeatCapacity_Whittington)
    return print(
        io,
        "T-dependent heat capacity following Whittington et al. (2009) for average crust. \n",
    )
end
#-------------------------------------------------------------------------



#-------------------------------------------------------------------------
# Heat capacity from phase diagram

# to be implemented - see density implementation
#

#-------------------------------------------------------------------------

# Help info for the calculation routines
"""
    Cp = compute_heatcapacity(s:<AbstractHeatCapacity, P, T)

Returns the heat capacity `Cp` at any temperature `T` and pressure `P` using any of the heat capacity laws implemented.

Currently available:
- ConstantHeatCapacity
- T\\_HeatCapacity_Whittington
- Latent_HeatCapacity

# Example
Using dimensional units
```julia
julia> T = 250.0:100:1250
julia> Cp2 = T_HeatCapacity_Whittington()
julia> Cp = similar(T)
julia> args = (; T=T)
julia> Cp =compute_heatcapacity!(Cp, Cp2, args)
11-element Vector nitful.Quantity{Float64, 𝐋² 𝚯⁻¹ 𝐓⁻², Unitful.FreeUnits{(kg⁻¹, J, K⁻¹), 𝐋² 𝚯⁻¹ 𝐓⁻², nothing}}}:
635.4269997294616 J kg⁻¹ K⁻¹
850.7470171764261 J kg⁻¹ K⁻¹
962.0959598489883 J kg⁻¹ K⁻¹
1037.542043377064 J kg⁻¹ K⁻¹
1097.351792196648 J kg⁻¹ K⁻¹
1149.274556367170 J kg⁻¹ K⁻¹
1157.791505094840 J kg⁻¹ K⁻¹
1172.355487419726 J kg⁻¹ K⁻¹
1186.919469744596 J kg⁻¹ K⁻¹
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
# compute_heatcapacity!(Cp_array::AbstractArray{_T}, s::AbstractHeatCapacity, T::AbstractArray{_T}) where _T =  compute_heatcapacity!(Cp_array,s,similar(T), T)

"""
    Cp = ComputeHeatCapacity(s::ConstantHeatCapacity)

Returns heat capacity if we are sure that we will only employ constant heat capacity in the simulation
"""

# Computational routines needed for computations with the MaterialParams structure
function compute_heatcapacity(s::AbstractMaterialParamsStruct, args)
    return compute_heatcapacity(s.HeatCapacity[1],args)
end

"""
    compute_heatcapacity!(Cp::AbstractArray{<:AbstractFloat}, MatParam::AbstractArray{<:AbstractMaterialParamsStruct}, Phases::AbstractArray{<:Integer}, P::AbstractArray{<:AbstractFloat},T::AbstractArray{<:AbstractFloat})

In-place computation of heat capacity `Cp` for the whole domain and all phases, in case a vector with phase properties `MatParam` is provided, along with `P` and `T` arrays.
This assumes that the `Phase` of every point is specified as an Integer in the `Phases` array.
"""
compute_heatcapacity!()

compute_heatcapacity(args::Vararg{Any, N}) where N = compute_param(compute_heatcapacity, args...)
compute_heatcapacity!(args::Vararg{Any, N}) where N = compute_param!(compute_heatcapacity, args...)

# In case just temperature is provided
#=
function compute_heatcapacity!(
    Cp::AbstractArray{_T,ndim},
    MatParam::NTuple{N,AbstractMaterialParamsStruct},
    Phases::AbstractArray{_I,ndim},
    T::AbstractArray{_T,ndim},
) where {_T,ndim,N,_I<:Integer}
    return compute_param!(compute_heatcapacity, Cp, MatParam, Phases, nothing, T)
end
=#

# extractor methods
for type in (ConstantHeatCapacity, T_HeatCapacity_Whittington, Latent_HeatCapacity)
    @extractors(type, :HeatCapacity)
end

end
