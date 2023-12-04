module Shearheating

# This implements different methods to specify Shearheating of rocks
#
# If you want to add a new method here, feel free to do so.
# Remember to also export the function name in GeoParams.jl (in addition to here)

using Parameters, LaTeXStrings, Unitful
using ..Units
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct
import Base.show, GeoParams.param_info
using ..MaterialParameters: MaterialParamsInfo

abstract type AbstractShearheating{T} <: AbstractMaterialParam end

export ConstantShearheating,   # constant
    compute_shearheating,   # calculation routines
    compute_shearheating!,
    param_info

include("../Computations.jl")
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
@with_kw_noshow struct ConstantShearheating{T,U} <: AbstractShearheating{T}
    Χ::GeoUnit{T,U} = 0.0 * NoUnits
end
ConstantShearheating(args...) = ConstantShearheating(convert.(GeoUnit, args)...)

function param_info(s::ConstantShearheating) # info about the struct
    return MaterialParamsInfo(;
        Equation=L"\H_s = \Chi \tau_{ij}(\dot{\varepsilon}_{ij} - \dot{\varepsilon}^{el}_{ij})",
    )
end

# In-place routine
function compute_shearheating!(H_s::AbstractArray, s::ConstantShearheating, τ::NTuple{N, AbstractArray}, ε::NTuple{N, AbstractArray}, ε_el::Union{Nothing, NTuple{N, AbstractArray}}) where N
    V = Val(N)

    @inline f(x, i) = ntuple(j -> x[j][i], V)
    @inline f(::Nothing, i) = nothing

    for i in eachindex(H_s)
        τ_i    = f(τ, i)
        ε_i    = f(ε, i)
        ε_el_i = f(ε_el, i)
        H_s[i]  = compute_shearheating(s, τ_i, ε_i, ε_el_i)
    end
end

# Calculation routine
@inline function compute_shearheating(s::ConstantShearheating, τ, ε, ε_el)
    @unpack_val Χ = s
    H_s = Χ * _compute_shearheating(τ, ε, ε_el)
    return H_s
end


@inline _compute_shearheating(τ, ε, ::Nothing) = sum(τi * εi for (τi, εi) in zip(τ, ε))
@inline _compute_shearheating(τ, ε, ε_el) = sum(τi * (εi - εi_el) for (τi, εi, εi_el) in zip(τ, ε, ε_el))
@inline _compute_shearheating(τ::NTuple{N,T}, ε::NTuple{N,T}, ε_el::NTuple{N,T}) where {N,T} = @.(τ * (ε - ε_el)) |> sum
@inline _compute_shearheating(τ::NTuple{N,T}, ε::NTuple{N,T}, ::Nothing) where {N,T} = @.(τ * ε) |> sum

# Print info
function show(io::IO, g::ConstantShearheating)
    return print(io, "Shear heating: H_s = $(UnitValue(g.Χ)) τ_ij*(ε_ij - ε^el_ij)")
end
#-------------------------------------------------------------------------

# Help info for the calculation routines
"""
    H_s = compute_shearheating(s:<AbstractShearheating, τ, ε, ε_el)

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
compute_shearheating(s::AbstractShearheating, τ, ε, ε_el)

"""
    H_s = ComputeShearheating(s:<AbstractShearheating, τ, ε)

Computes the shear heating source term when there is no elasticity

```math
H_s = \\Chi \\cdot \\tau_{ij}  \\dot{\\varepsilon}_{ij}
```

# Parameters
- ``\\Chi`` : The efficiency of shear heating (between 0-1)
- ``\\tau_{ij}`` : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
- ``\\dot{\\varepsilon}_{ij}`` : The full deviatoric strainrate tensor
"""
@inline function compute_shearheating(s::AbstractShearheating{_T}, τ::Any, ε::Any) where {_T}
    return compute_shearheating(s, τ, ε, nothing)
end

"""
    compute_shearheating!(H_s, s:<AbstractShearheating,  τ, ε, ε_el)

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
compute_shearheating!(H_s, s::AbstractShearheating, τ, ε, ε_el)

"""
    compute_shearheating!(H_s, s:<AbstractShearheating, τ, ε)

Computes the shear heating source term `H_s` in-place when there is no elasticity

```math
H_s = \\Chi \\cdot \\tau_{ij}  \\dot{\\varepsilon}_{ij}
```

# Parameters
- ``\\Chi`` : The efficiency of shear heating (between 0-1)
- ``\\tau_{ij}`` : The full deviatoric stress tensor [4 components in 2D; 9 in 3D]
- ``\\dot{\\varepsilon}_{ij}`` : The full deviatoric strainrate tensor
"""
@inline function compute_shearheating!(H_s, s::AbstractShearheating, τ, ε)
    return compute_shearheating!(H_s, s, τ, ε, nothing)
end

function compute_shearheating(s::AbstractMaterialParamsStruct, args::Vararg{Any,N}) where N
    if isempty(s.ShearHeat)
        return isempty(args) ? 0.0 : zero(typeof(args).types[1])  # return zero if not specified
    else
        return compute_shearheating(s.ShearHeat[1], args...)
    end
end

compute_shearheating(args::Vararg{Any,N}) where N = compute_param(compute_shearheating, args...)
compute_shearheating!(args::Vararg{Any,N}) where N = compute_param!(compute_shearheating!, args...)

end
