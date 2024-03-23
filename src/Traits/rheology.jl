abstract type RheologyTrait end
struct LinearRheologyTrait    <: RheologyTrait end
struct NonLinearRheologyTrait <: RheologyTrait end

# traits individual rheologies
@inline islinear(::AbstractElasticity)      = LinearRheologyTrait()
@inline islinear(::LinearViscous)           = LinearRheologyTrait()
@inline islinear(::AbstractConstitutiveLaw) = NonLinearRheologyTrait()
@inline islinear(::T) where T               = throw(ArgumentError("$T is an unsupported rheology type"))

# compares two rheologies and return linear trait only and if only both are linear
@inline islinear(::LinearRheologyTrait, ::LinearRheologyTrait) = LinearRheologyTrait()
@inline islinear(::RheologyTrait, ::RheologyTrait)   = NonLinearRheologyTrait()
@inline islinear(v1::AbstractConstitutiveLaw, v2::AbstractConstitutiveLaw) = islinear(islinear(v1), islinear(v2))

# traits for composite rheologies
@inline islinear(c::CompositeRheology) = islinear(c.elements)
# traits for MaterialParams
@inline islinear(r::MaterialParams)    = islinear(r.CompositeRheology...)

# recursively (pairwise, right-to-left) compare rheology traits of a composite or tuple of material params
@inline islinear(r::NTuple{N, Union{AbstractConstitutiveLaw, MaterialParams}}) where N = islinear(islinear(first(r)), islinear(Base.tail(r)))    
@inline islinear(v::NTuple{1, Union{AbstractConstitutiveLaw, MaterialParams}})         = islinear(v...)
