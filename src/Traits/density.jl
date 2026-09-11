"""
    DensityTrait

Supertype of the Holy traits classifying a density parameterization as either
`ConstantDensityTrait` or `NonConstantDensityTrait`, returned by [`isconstant`](@ref).
"""
abstract type DensityTrait end
"[`DensityTrait`](@ref) marking a constant density; returned by [`isconstant`](@ref)."
struct ConstantDensityTrait <: DensityTrait end
"[`DensityTrait`](@ref) marking a pressure/temperature-dependent density; returned by [`isconstant`](@ref)."
struct NonConstantDensityTrait <: DensityTrait end

# traits individual densities
"""
    isconstant(v) -> DensityTrait

Returns [`ConstantDensityTrait`](@ref) if the density parameterization `v` (or `MaterialParams`,
or tuple thereof) is pressure- and temperature-independent, and `NonConstantDensityTrait`
otherwise.
"""
@inline isconstant(::ConstantDensity) = ConstantDensityTrait()
@inline isconstant(::AbstractDensity) = NonConstantDensityTrait()
@inline isconstant(::AbstractPhaseDiagramsStruct) = NonConstantDensityTrait()
@inline isconstant(::T) where {T} = throw(ArgumentError("$T is an unsupported density type"))

# compares two densities and return constant trait only and if only both are constant
@inline isconstant(::ConstantDensityTrait, ::ConstantDensityTrait) = ConstantDensityTrait()
@inline isconstant(::DensityTrait, ::DensityTrait) = NonConstantDensityTrait()

# traits for MaterialParams
@inline isconstant(r::MaterialParams) = isconstant(r.Density...)

# recursively (pairwise, right-to-left) compare density traits of a tuple of material params
@inline isconstant(r::NTuple{N, MaterialParams}) where {N} = isconstant(isconstant(first(r)), isconstant(Base.tail(r)))
@inline isconstant(v::NTuple{1, MaterialParams}) = isconstant(v...)
