module Density

# This implements different methods to compute density
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)


using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import Base.show

abstract type AbstractDensity <: AbstractMaterialParam end

export  ComputeDensity,         # calculation routines
        ConstantDensity,        # constant
        PT_Density              # P & T dependent density

# Constant Density -------------------------------------------------------
"""
    ConstantDensity(ρ=2900kg/m3)
    
Set a constant density:
```math  
    \\rho  = cst
```
where ``\\rho`` is the density [``kg/m^3``].
"""
@with_kw_noshow mutable struct ConstantDensity <: AbstractDensity
    equation::LaTeXString   =   L"\rho = cst"     
    ρ::GeoUnit              =   2900kg/m^3                # density
end

# Calculation routine
function ComputeDensity(P,T, s::ConstantDensity)
    @unpack ρ   = s
    
    return ρ*1.0
end

# Print info 
function show(io::IO, g::ConstantDensity)  
    print(io, "Constant density: ρ=$(g.ρ.val)")  
end
#-------------------------------------------------------------------------

# Pressure & Temperature dependent density -------------------------------
"""
    PT_Density(ρ0=2900kg/m3, α=3e-5/K, β=1e-9/Pa, T0=0C, P=0MPa)
    
Set a pressure and temperature-dependent density:
```math  
    \\rho  = \\rho_0 (1.0 - \\alpha (T-T_0) + \\beta  (P-P_0) )  
```
where ``\\rho_0`` is the density [``kg/m^3``] at reference temperature ``T_0`` and pressure ``P_0``,
``\\alpha`` is the temperature dependence of density and ``\\beta`` the pressure dependence.

"""
@with_kw_noshow mutable struct PT_Density <: AbstractDensity
    equation::LaTeXString   =   L"\rho = \rho_0(1.0-\alpha (T-T_0) + \beta (P-P_0)"     
    ρ0::GeoUnit             =   2900kg/m^3                  # density
    α::GeoUnit              =   3e-5/K                      # T-dependence of density
    β::GeoUnit              =   1e-9/Pa                     # P-dependence of density
    T0::GeoUnit             =   0C                          # Reference temperature
    P0::GeoUnit             =   0MPa                        # Reference pressure
end

# Calculation routine
function ComputeDensity(P,T, s::PT_Density)
    @unpack ρ0,α,β,P0, T0   = s
    
    ρ = ρ0*(1.0 - α*(T-T0) + β*(P-P0) )

    return ρ
end

# Print info 
function show(io::IO, g::PT_Density)  
    print(io, "P/T-dependent density: ρ0=$(g.ρ0.val), α=$(g.α.val), β=$(g.β.val), T0=$(g.T0.val), P0=$(g.P0.val)")  
end
#-------------------------------------------------------------------------



# Help info for the calculation routines
"""
ComputeDensity(P,T, s:<AbstractDensity)

Returns the density ``ρ`` at a given pressure and temperature using any 
of the density EoS implemented.

"""
ComputeDensity



end