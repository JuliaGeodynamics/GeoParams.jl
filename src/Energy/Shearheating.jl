module Shearheating

# This implements different methods to specify Shearheating of rocks
#
# If you want to add a new method here, feel free to do so. 
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam
import Base.show

abstract type AbstractShearheating{T} <: AbstractMaterialParam end

export  ComputeShearheating, ComputeShearheating!,  # calculation routines
        ConstantShearheating                        # constant
        
# Constant Shearheating -------------------------------------------------------
"""
    ConstantShearheating(Χ=0.0NoUnits)
    
Set the shear heating efficiency [0-1] parameter
```math  
Χ  = cst
```
where ``\\Chi`` is the shear heating efficiency [NoUnits]

Shear heating is computed as
```math  
H_s = \\Chi \\cdot \\tau_{ij}(\\dot{\\varepsilon}_{ij} - \\dot{\\varepsilon}^{el}_{ij})
```

"""
@with_kw_noshow struct ConstantShearheating{T} <: AbstractShearheating{T}
    equation::LaTeXString   =   L"\H_s = \Chi \tau_{ij}(\dot{\varepsilon}_{ij} - \dot{\varepsilon}^{el}_{ij})"     
    Χ::GeoUnit{T}              =   0NoUnits               
end
ConstantShearheating(a...) = ConstantShearheating{Float64}(a...)

ConstantShearheating(x::Number) = ConstantShearheating(Χ = x*NoUnits)   # allow 

# In-place routine
function ComputeShearheating!(H_s, τ, ε, ε_el, s::ConstantShearheating)
    @unpack Χ   = s

    if isnothing(ε_el)
        H_s = NumValue(Χ)*sum( τ .* ε  )
    else
        H_s = NumValue(Χ)*sum( τ .* (ε .- ε_el) )
    end

end

# Calculation routine
function ComputeShearheating(τ, ε, ε_el, s::ConstantShearheating)
    @unpack Χ   = s
    
    if isnothing(ε_el)
        H_s = NumValue(Χ)*sum( τ .* ε  )
    else
        H_s = NumValue(Χ)*sum( τ .* (ε .- ε_el) )
    end

    return H_s
end

# Print info 
function show(io::IO, g::ConstantShearheating)  
    print(io, "Shear heating: H_s = $(g.Χ.val) τ_ij*(ε_ij - ε^el_ij)")   
end
#-------------------------------------------------------------------------

# Help info for the calculation routines
"""
    H_s = ComputeShearheating(τ, ε, ε_el,  s:<AbstractShearheating)

Computes the shear heating source term

```math  
H_s = \\Chi \\cdot \\tau_{ij} ( \\dot{\\varepsilon}_{ij} - \\dot{\\varepsilon}^{el}_{ij})
```

# Parameters
- ``\\Chi`` : The efficiency of shear heating (between 0-1)
- ``\\tau_{ij}`` : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
- ``\\dot{\\varepsilon}_{ij}`` : The full deviatoric strainrate tensor
- ``\\dot{\\varepsilon}^{el}_{ij}`` : The full elastic deviatoric strainrate tensor

"""
ComputeShearheating(τ, ε, ε_el, s::AbstractShearheating)


"""
    H_s = ComputeShearheating(τ, ε,  s:<AbstractShearheating)

Computes the shear heating source term when there is no elasticity

```math  
H_s = \\Chi \\cdot \\tau_{ij}  \\dot{\\varepsilon}_{ij} 
```

# Parameters
- ``\\Chi`` : The efficiency of shear heating (between 0-1)
- ``\\tau_{ij}`` : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
- ``\\dot{\\varepsilon}_{ij}`` : The full deviatoric strainrate tensor
"""
ComputeShearheating(τ::Any, ε::Any, s::AbstractShearheating) = ComputeShearheating(τ, ε, nothing, s::ConstantShearheating)


"""
    ComputeShearheating!(H_s, τ, ε, ε_el,  s:<AbstractShearheating)

Computes the shear heating source term in-place

```math  
H_s = \\Chi \\cdot \\tau_{ij} ( \\dot{\\varepsilon}_{ij} - \\dot{\\varepsilon}^{el}_{ij})
```

# Parameters
- ``\\Chi`` : The efficiency of shear heating (between 0-1)
- ``\\tau_{ij}`` : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
- ``\\dot{\\varepsilon}_{ij}`` : The full deviatoric strainrate tensor
- ``\\dot{\\varepsilon}^{el}_{ij}`` : The full elastic deviatoric strainrate tensor

*NOTE:*
The shear heating terms require the full deviatoric stress & strain rate tensors, i.e.:

```math  
2D: \\tau_{ij} = \\left(   
                \\begin{matrix}
                    \\tau_{xx} & \\tau_{xz} \\\\
                    \\tau_{zx} & \\tau_{zz} 
                \\end{matrix} 
            \\right)
```
Since ``\\tau_{zx}=\\tau_{xz}``, most geodynamic codes only take one of the terms into account; shear heating requires all components to be used! 
"""
ComputeShearheating!(H_s, τ, ε, ε_el, s::AbstractShearheating)

"""
    ComputeShearheating!(H_s, τ, ε,  s:<AbstractShearheating)

Computes the shear heating source term `H_s` in-place when there is no elasticity

```math  
H_s = \\Chi \\cdot \\tau_{ij}  \\dot{\\varepsilon}_{ij} 
```

# Parameters
- ``\\Chi`` : The efficiency of shear heating (between 0-1)
- ``\\tau_{ij}`` : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
- ``\\dot{\\varepsilon}_{ij}`` : The full deviatoric strainrate tensor
"""
ComputeShearheating!(H_s::Any, τ::Any, ε::Any, s::AbstractShearheating) = ComputeShearheating!(H_s, τ, ε, nothing, s::ConstantShearheating)

end