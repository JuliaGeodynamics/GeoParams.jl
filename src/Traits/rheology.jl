abstract type RheologyTrait end
struct LinearRheologyTrait     <: RheologyTrait end
struct NonLinearRheologyTrait  <: RheologyTrait end
struct PlasticRheologyTrait   <: RheologyTrait end
struct NonPlasticRheologyTrait <: RheologyTrait end

## LINEAR RHEOLOGY traits

# traits individual rheologies
@inline islinear(::AbstractElasticity)      = LinearRheologyTrait()
@inline islinear(::LinearViscous)           = LinearRheologyTrait()
@inline islinear(::AbstractConstitutiveLaw) = NonLinearRheologyTrait()
@inline islinear(::T) where T               = throw(ArgumentError("$T is an unsupported rheology type"))

# compares two rheologies and return linear trait only and if only both are linear
@inline islinear(::LinearRheologyTrait, ::LinearRheologyTrait) = LinearRheologyTrait()
@inline islinear(::RheologyTrait, ::RheologyTrait)             = NonLinearRheologyTrait()
@inline islinear(v1::AbstractConstitutiveLaw, v2::AbstractConstitutiveLaw) = islinear(islinear(v1), islinear(v2))

# traits for composite rheologies
@inline islinear(c::CompositeRheology) = islinear(c.elements)
# traits for MaterialParams
@inline islinear(r::MaterialParams)    = islinear(r.CompositeRheology...)

# recursively (pairwise, right-to-left) compare rheology traits of a composite or tuple of material params
@inline islinear(r::NTuple{N, Union{AbstractConstitutiveLaw, MaterialParams}}) where N = islinear(islinear(first(r)), islinear(Base.tail(r)))    
@inline islinear(v::NTuple{1, Union{AbstractConstitutiveLaw, MaterialParams}})         = islinear(v...)

## PLASTIC RHEOLOGY TRAITS

@inline isplastic(::AbstractPlasticity)      = PlasticRheologyTrait()
@inline isplastic(::AbstractConstitutiveLaw) = NonPlasticRheologyTrait()
@inline isplastic(::T) where T               = throw(ArgumentError("$T is an unsupported rheology type"))

# compares two rheologies and return linear trait only and if only both are linear
@inline isplastic(::NonPlasticRheologyTrait, ::NonPlasticRheologyTrait) = NonPlasticRheologyTrait()
@inline isplastic(::RheologyTrait, ::RheologyTrait)                     = PlasticRheologyTrait()
@inline isplastic(v1::AbstractConstitutiveLaw, v2::AbstractConstitutiveLaw) = isplastic(isplastic(v1), isplastic(v2))

# traits for composite rheologies
@inline isplastic(c::CompositeRheology) = isplastic(c.elements)
# traits for MaterialParams
@inline isplastic(r::MaterialParams)    = isplastic(r.CompositeRheology...)

# recursively (pairwise, right-to-left) compare rheology traits of a composite or tuple of material params
@inline isplastic(r::NTuple{N, Union{AbstractConstitutiveLaw, MaterialParams}}) where N = isplastic(isplastic(first(r)), isplastic(Base.tail(r)))    
@inline isplastic(v::NTuple{1, Union{AbstractConstitutiveLaw, MaterialParams}})         = isplastic(v...)
