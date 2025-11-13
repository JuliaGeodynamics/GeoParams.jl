abstract type DensityTrait end
struct ConstantDensityTrait <: DensityTrait end
struct NonConstantDensityTrait <: DensityTrait end

# traits individual densities
@inline isconstant(::ConstantDensity) = ConstantDensityTrait()
@inline isconstant(::AbstractDensity) = NonConstantDensityTrait()
@inline isconstant(::PhaseDiagram_LookupTable) = NonConstantDensityTrait()
@inline isconstant(::T) where {T} = throw(ArgumentError("$T is an unsupported density type"))

# compares two densities and return constant trait only and if only both are constant
@inline isconstant(::ConstantDensityTrait, ::ConstantDensityTrait) = ConstantDensityTrait()
@inline isconstant(::DensityTrait, ::DensityTrait) = NonConstantDensityTrait()

# traits for MaterialParams
@inline isconstant(r::MaterialParams) = isconstant(r.Density...)

# recursively (pairwise, right-to-left) compare density traits of a tuple of material params
@inline isconstant(r::NTuple{N, MaterialParams}) where {N} = isconstant(isconstant(first(r)), isconstant(Base.tail(r)))
@inline isconstant(v::NTuple{1, MaterialParams}) = isconstant(v...)
