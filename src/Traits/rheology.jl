abstract type RheologyTrait end
struct LinearRheologyTrait     <: RheologyTrait end
struct NonLinearRheologyTrait  <: RheologyTrait end
struct ElasticRheologyTrait    <: RheologyTrait end
struct NonElasticRheologyTrait <: RheologyTrait end
struct PlasticRheologyTrait    <: RheologyTrait end
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

## ELASTICITY RHEOLOGY traits

# traits individual rheologies
@inline isviscoelastic(::AbstractElasticity)      = ElasticRheologyTrait()
@inline isviscoelastic(::AbstractConstitutiveLaw) = NonElasticRheologyTrait()
@inline isviscoelastic(::T) where T               = throw(ArgumentError("$T is an unsupported rheology type"))

# compares two rheologies and return linear trait only and if only both are linear
@inline isviscoelastic(::RheologyTrait, ::ElasticRheologyTrait) = ElasticRheologyTrait()
@inline isviscoelastic(::ElasticRheologyTrait, ::RheologyTrait) = ElasticRheologyTrait()
@inline isviscoelastic(::RheologyTrait, ::RheologyTrait)        = NonElasticRheologyTrait()
@inline isviscoelastic(v1::Union{AbstractConstitutiveLaw, AbstractPlasticity}, v2::Union{AbstractConstitutiveLaw, AbstractPlasticity}) = isviscoelastic(isviscoelastic(v1), isviscoelastic(v2))

# traits for composite rheologies
@inline isviscoelastic(c::CompositeRheology) = isviscoelastic(c.elements)
# traits for MaterialParams
@inline isviscoelastic(r::MaterialParams)    = isviscoelastic(r.CompositeRheology...)

# recursively (pairwise, right-to-left) compare rheology traits of a composite or tuple of material params
@inline isviscoelastic(r::NTuple{N, Union{AbstractConstitutiveLaw, AbstractPlasticity, MaterialParams}}) where N = isviscoelastic(isviscoelastic(first(r)), isviscoelastic(Base.tail(r)))    
@inline isviscoelastic(v::NTuple{1, Union{AbstractConstitutiveLaw, AbstractPlasticity, MaterialParams}})         = isviscoelastic(v...)

## PLASTIC RHEOLOGY TRAITS

@inline isplasticity(::AbstractPlasticity)      = PlasticRheologyTrait()
@inline isplasticity(::AbstractConstitutiveLaw) = NonPlasticRheologyTrait()
@inline isplasticity(::T) where T               = throw(ArgumentError("$T is an unsupported rheology type"))

# compares two rheologies and return linear trait only and if only both are linear
@inline isplasticity(::NonPlasticRheologyTrait, ::NonPlasticRheologyTrait) = NonPlasticRheologyTrait()
@inline isplasticity(::RheologyTrait, ::RheologyTrait)                     = PlasticRheologyTrait()
@inline isplasticity(v1::AbstractConstitutiveLaw, v2::AbstractConstitutiveLaw) = isplasticity(isplasticity(v1), isplasticity(v2))

# traits for composite rheologies
@inline isplasticity(c::CompositeRheology) = isplasticity(c.elements)
# traits for MaterialParams
@inline isplasticity(r::MaterialParams)    = isplasticity(r.CompositeRheology...)

# recursively (pairwise, right-to-left) compare rheology traits of a composite or tuple of material params
@inline isplasticity(r::NTuple{N, Union{AbstractConstitutiveLaw, MaterialParams}}) where N = isplasticity(isplasticity(first(r)), isplasticity(Base.tail(r)))    
@inline isplasticity(v::NTuple{1, Union{AbstractConstitutiveLaw, MaterialParams}})         = isplasticity(v...)
