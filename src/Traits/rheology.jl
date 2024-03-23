abstract type RheologyTrait end
struct LinearRheology <: RheologyTrait end
struct NonLinearRheology <: RheologyTrait end

# traits individual rheologies
islinear(::AbstractElasticity)      = LinearRheology()
islinear(::LinearViscous)           = LinearRheology()
islinear(::AbstractConstitutiveLaw) = NonLinearRheology()
islinear(::T) where T               = throw(ArgumentError("$T is an unsupported rheology type"))

# compares two rheologies and return linear trait only and if only both are linear
islinear(::LinearRheology, ::LinearRheology) = LinearRheology()
islinear(::RheologyTrait, ::RheologyTrait) = NonLinearRheology()
islinear(v1::AbstractConstitutiveLaw, v2::AbstractConstitutiveLaw) = islinear(islinear(v1), islinear(v2))

# traits for composite rheologies
islinear(c::CompositeRheology) = islinear(c.elements)
# traits for MaterialParams
islinear(r::MaterialParams) = islinear(r.CompositeRheology...)

# recursively (pairwise, right-to-left) compare rheology traits of a composite or tuple of material params
islinear(r::NTuple{N, Union{AbstractConstitutiveLaw, MaterialParams}}) where N = islinear(islinear(r[1]), islinear(Base.tail(r)))    
islinear(v::NTuple{1, Union{AbstractConstitutiveLaw, MaterialParams}}) = islinear(v[1])
