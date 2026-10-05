"""
    RheologyTrait

Supertype of the Holy traits used to classify a rheology. The concrete singletons come in
mutually exclusive pairs — `LinearRheologyTrait`/`NonLinearRheologyTrait`,
`ElasticRheologyTrait`/`NonElasticRheologyTrait`, `PlasticRheologyTrait`/`NonPlasticRheologyTrait` —
and are returned by [`islinear`](@ref), [`isviscoelastic`](@ref), and [`isplasticity`](@ref) to
enable branch-free dispatch on a rheology's properties.
"""
abstract type RheologyTrait end
"[`RheologyTrait`](@ref) marking a linear rheology; returned by [`islinear`](@ref)."
struct LinearRheologyTrait <: RheologyTrait end
"[`RheologyTrait`](@ref) marking a non-linear rheology; returned by [`islinear`](@ref)."
struct NonLinearRheologyTrait <: RheologyTrait end
"[`RheologyTrait`](@ref) marking a rheology that contains elasticity; returned by [`isviscoelastic`](@ref)."
struct ElasticRheologyTrait <: RheologyTrait end
"[`RheologyTrait`](@ref) marking a rheology with no elasticity; returned by [`isviscoelastic`](@ref)."
struct NonElasticRheologyTrait <: RheologyTrait end
"[`RheologyTrait`](@ref) marking a rheology that contains plasticity; returned by [`isplasticity`](@ref)."
struct PlasticRheologyTrait <: RheologyTrait end
"[`RheologyTrait`](@ref) marking a rheology with no plasticity; returned by [`isplasticity`](@ref)."
struct NonPlasticRheologyTrait <: RheologyTrait end

## LINEAR RHEOLOGY traits

# traits individual rheologies
"""
    islinear(v) -> RheologyTrait

Returns [`LinearRheologyTrait`](@ref) if the rheology `v` (a constitutive law, [`CompositeRheology`](@ref),
`MaterialParams`, or tuple thereof) is linear, and `NonLinearRheologyTrait` otherwise. When given
two arguments, returns the linear trait only if both are linear.
"""
@inline islinear(::AbstractElasticity) = LinearRheologyTrait()
@inline islinear(::LinearViscous) = LinearRheologyTrait()
@inline islinear(::AbstractConstitutiveLaw) = NonLinearRheologyTrait()
@inline islinear(::T) where {T} = throw(ArgumentError("$T is an unsupported rheology type"))

# compares two rheologies and return linear trait only and if only both are linear
@inline islinear(::LinearRheologyTrait, ::LinearRheologyTrait) = LinearRheologyTrait()
@inline islinear(::RheologyTrait, ::RheologyTrait) = NonLinearRheologyTrait()
@inline islinear(v1::AbstractConstitutiveLaw, v2::AbstractConstitutiveLaw) = islinear(islinear(v1), islinear(v2))

# traits for composite rheologies
@inline islinear(c::CompositeRheology) = islinear(c.elements)
# traits for MaterialParams
@inline islinear(r::MaterialParams) = islinear(r.CompositeRheology...)

# recursively (pairwise, right-to-left) compare rheology traits of a composite or tuple of material params
@inline islinear(r::NTuple{N, Union{AbstractConstitutiveLaw, MaterialParams}}) where {N} = islinear(islinear(first(r)), islinear(Base.tail(r)))
@inline islinear(v::NTuple{1, Union{AbstractConstitutiveLaw, MaterialParams}}) = islinear(v...)

## ELASTICITY RHEOLOGY traits

# traits individual rheologies
"""
    isviscoelastic(v) -> RheologyTrait

Returns [`ElasticRheologyTrait`](@ref) if the rheology `v` contains an elastic element, and
`NonElasticRheologyTrait` otherwise. When given two arguments, returns the elastic trait if either
is elastic.
"""
@inline isviscoelastic(::AbstractElasticity) = ElasticRheologyTrait()
@inline isviscoelastic(::AbstractConstitutiveLaw) = NonElasticRheologyTrait()
@inline isviscoelastic(::T) where {T} = throw(ArgumentError("$T is an unsupported rheology type"))

# compares two rheologies and return linear trait only and if only both are linear
@inline isviscoelastic(::RheologyTrait, ::ElasticRheologyTrait) = ElasticRheologyTrait()
@inline isviscoelastic(::ElasticRheologyTrait, ::RheologyTrait) = ElasticRheologyTrait()
@inline isviscoelastic(::ElasticRheologyTrait, ::ElasticRheologyTrait) = ElasticRheologyTrait()
@inline isviscoelastic(::RheologyTrait, ::RheologyTrait) = NonElasticRheologyTrait()
@inline isviscoelastic(v1::Union{AbstractConstitutiveLaw, AbstractPlasticity}, v2::Union{AbstractConstitutiveLaw, AbstractPlasticity}) = isviscoelastic(isviscoelastic(v1), isviscoelastic(v2))

# traits for composite rheologies
@inline isviscoelastic(c::CompositeRheology) = isviscoelastic(c.elements)
# traits for MaterialParams
@inline isviscoelastic(r::MaterialParams) = isviscoelastic(r.CompositeRheology...)

# recursively (pairwise, right-to-left) compare rheology traits of a composite or tuple of material params
@inline isviscoelastic(r::NTuple{N, Union{AbstractConstitutiveLaw, AbstractPlasticity, MaterialParams}}) where {N} = isviscoelastic(isviscoelastic(first(r)), isviscoelastic(Base.tail(r)))
@inline isviscoelastic(v::NTuple{1, Union{AbstractConstitutiveLaw, AbstractPlasticity, MaterialParams}}) = isviscoelastic(v...)

## PLASTIC RHEOLOGY TRAITS

"""
    isplasticity(v) -> RheologyTrait

Returns [`PlasticRheologyTrait`](@ref) if the rheology `v` contains a plastic element, and
`NonPlasticRheologyTrait` otherwise.
"""
@inline isplasticity(::AbstractPlasticity) = PlasticRheologyTrait()
@inline isplasticity(::AbstractConstitutiveLaw) = NonPlasticRheologyTrait()
@inline isplasticity(::T) where {T} = throw(ArgumentError("$T is an unsupported rheology type"))

# compares two rheologies and return linear trait only and if only both are linear
@inline isplasticity(::NonPlasticRheologyTrait, ::NonPlasticRheologyTrait) = NonPlasticRheologyTrait()
@inline isplasticity(::RheologyTrait, ::RheologyTrait) = PlasticRheologyTrait()
@inline isplasticity(v1::AbstractConstitutiveLaw, v2::AbstractConstitutiveLaw) = isplasticity(isplasticity(v1), isplasticity(v2))

# traits for composite rheologies
@inline isplasticity(c::CompositeRheology) = isplasticity(c.elements)
# traits for MaterialParams
@inline isplasticity(r::MaterialParams) = isplasticity(r.CompositeRheology...)

# recursively (pairwise, right-to-left) compare rheology traits of a composite or tuple of material params
@inline isplasticity(r::NTuple{N, Union{AbstractConstitutiveLaw, MaterialParams}}) where {N} = isplasticity(isplasticity(first(r)), isplasticity(Base.tail(r)))
@inline isplasticity(v::NTuple{1, Union{AbstractConstitutiveLaw, MaterialParams}}) = isplasticity(v...)
