module Conductivity

# This implements different methods to specify conductivity of rocks
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import Base.show

abstract type AbstractConductivity <: AbstractMaterialParam end

export  ComputeConductivity,                # calculation routines
        ConstantConductivity,               # constant
        T_Conductivity_Whittacker           # T-dependent heat capacity

# Constant Conductivity -------------------------------------------------------
"""
    ConstantConductivity(k=3.0W/m/K)
    
Set a constant conductivity
```math  
    k  = cst
```
where ``k`` is the thermal conductivity [``W/m/K``].
"""
@with_kw_noshow mutable struct ConstantConductivity <: AbstractConductivity
    equation::LaTeXString   =   L"k = cst"     
    k::GeoUnit              =   3.0Watt/m/K               
end

# Calculation routine
function ComputeConductivity(P, T, s::ConstantConductivity)
    @unpack k   = s
    
    return k*1.0
end

# Print info 
function show(io::IO, g::ConstantConductivity)  
    print(io, "Constant conductivity: k=$(g.k.val)")  
end
#-------------------------------------------------------------------------


# Temperature dependent conductivity -------------------------------
"""
    T_Conductivity_Whittacker()
    
Sets a temperature-dependent conductivity following the parameterization of Whittacker et al. (2009), Nature. 
Their parameterization is originally given for the thermal diffusivity, together with a parameterization for thermal conductivity, which allows us 
```math  
    Cp = a + b T - c/T^2 
    \\kappa = d/T - e, if T<=846K
    \\kappa = f - g*T, if T>846K
    \\rho = 2700 kg/m3
    k = \\kappa*Cp*\\rho
```

where ``Cp`` is the heat capacity [``J/mol/K``], and ``a,b,c`` are parameters that dependent on the temperature `T`:
- a = 199.50 J/mol/K    if T<= 846 K
- a = 199.50 J/mol/K    if T> 846 K
- b = 0.0857J/mol/K^2   if T<= 846 K
- b = 0.0323J/mol/K^2   if T> 846 K
- c = 5e6J/mol*K        if T<= 846 K
- c = 47.9e-6J/mol*K    if T> 846 K
- d = 576.3m^2/s*K      
- e = 0.062m^2/s        
- f = 0.732m^2/s        
- g = 0.000135m^2/s/K 
"""
@with_kw_noshow mutable struct T_Conductivity_Whittacker <: AbstractConductivity
    # Note: the resulting curve was visually compared with Fig. 2 of the paper
    equation::LaTeXString   =   L"k = f(T) "     
    a0::GeoUnit             =   199.5J/mol/K                # prefactor for low T       (T<= 846 K)
    a1::GeoUnit             =   229.32J/mol/K               # prefactor for high T      (T>  846 K)
    b0::GeoUnit             =   0.0857J/mol/K^2             # linear term for low T     (T<= 846 K)
    b1::GeoUnit             =   0.0323J/mol/K^2             # linear term for high T    (T>  846 K)
    c0::GeoUnit             =   5e6J/mol*K                  # quadratic term for low T  (T<= 846 K)
    c1::GeoUnit             =   47.9e-6J/mol*K              # quadratic term for high T (T>  846 K)
    Tcutoff::GeoUnit        =   846K                        # cutoff temperature
    rho::GeoUnit              =   2700kg/m^3                  # Density they use for an average crust
    d::GeoUnit              =   576.3m^2/s*K                # diffusivity parameterization
    e::GeoUnit              =   0.062m^2/s                  # diffusivity parameterization
    f::GeoUnit              =   0.732m^2/s                  # diffusivity parameterization
    g::GeoUnit              =   0.000135m^2/s/K             # diffusivity parameterization
end

# Calculation routine
function ComputeConductivity(P,T, s::T_Conductivity_Whittacker)
    @unpack a0,a1,b0,b1,c0,c1,Tcutoff,rho,d,e,f,g   = s
    
    k  = zeros(size(T))*Watt/m/K
    ρ  = Value(rho)
    
    for i in eachindex(T)
        if T[i] <= Value(Tcutoff)
            a,b,c = Value(a0),Value(b0),Value(c0)
            κ     = Value(d)/T[i] - Value(e)  
        else
            a,b,c = Value(a1),Value(b1),Value(c1)
            κ     = Value(f) - Value(g)*T[i]
        end
       
        cp = a + b*T[i] - c/T[i]^2 # conductivity

        k[i] = κ*ρ*cp       # compute conductivity from diffusivity

    end

    return k
end

# Print info 
function show(io::IO, g::T_Conductivity_Whittacker)  
    print(io, "T-dependent conductivity following Whittacker et al. (2009) for average crust). \n");
end
#-------------------------------------------------------------------------



# Help info for the calculation routines
"""
    k = ComputeConductivity(P, T, s:<AbstractConductivity)

Returns the thermal conductivity `k` at any temperature `T` and pressure `P` using any of the parameterizations implemented.

Currently available:
- ConstantConductivity
- T\\_Conductivity_Whittacker

# Example 
Using dimensional units
```julia
julia> T  = (250:100:1250)*K;
julia> cp = T_HeatCapacity_Whittacker()
julia> Cp = ComputeHeatCapacity(0,T,cp)
```


"""
ComputeConductivity()


"""
    k = ComputeConductivity(T::Any, s::AbstractConductivity)

Computes conductivity if only temperature (and not pressure) is specified
"""
ComputeConductivity(T::Any, s::AbstractConductivity) =  ComputeConductivity(0,T, s)

"""
    k = ComputeConductivity(s::ConstantConductivity)

Returns conductivity if we are sure that we will only employ constant values throughout the simulation
"""
ComputeConductivity(s::ConstantConductivity) =  ComputeConductivity(0,0, s)


end