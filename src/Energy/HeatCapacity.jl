module HeatCapacity

# This implements different methods to specify heat capacity
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import Base.show

abstract type AbstractHeatCapacity <: AbstractMaterialParam end

export  ComputeHeatCapacity,                # calculation routines
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
    
    return cp*1.0
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
- m 

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
    m::GeoUnit              =   0.22178kg/mol               # average molar mass 
    Tcutoff::GeoUnit        =   846K                        # cutoff temperature
end

# Calculation routine
function ComputeHeatCapacity(P,T, s::T_HeatCapacity_Whittacker)
    @unpack a0,a1,b0,b1,c0,c1, m, Tcutoff   = s
    
    cp = zeros(size(T)).*Value(a0)./Value(m)

    for i in eachindex(T)
        if T[i] <= Value(Tcutoff)
            a,b,c = Value(a0),Value(b0),Value(c0)
        else
            a,b,c = Value(a1),Value(b1),Value(c1)
        end
       
        cp[i] = (a + b*T[i] - c/T[i]^2)/m 
    end

    return cp
end

# Print info 
function show(io::IO, g::T_HeatCapacity_Whittacker)  
    print(io, "T-dependent heat capacity: cp/$(g.m.val)=$(g.a0.val) + $(g.b0.val)*T - $(g.c0.val)/T^2 (for T<=$(g.Tcutoff.val)); ");
    print(io, " cp/$(g.m.val)=$(g.a1.val) + $(g.b1.val)*T - $(g.c1.val)/T^2 (for T>$(g.Tcutoff.val)) \n");
end
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


end