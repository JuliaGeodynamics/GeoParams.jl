# extract elements from composite rheology
@inline elements(v::CompositeRheology) = v.elements

# compute viscosity given second invariants of strain rate and deviatoric stress tensors
@inline _viscosity(τII, εII) = τII / (2  * εII)

# compute effective "creep" viscosity from strain rate tensor
function compute_viscosity_εII(v::AbstractCreepLaw, εII, args::Vararg{T, N}) where {T, N}
    τII = compute_τII(v, εII, args...)
    η = _viscosity(τII, εII)
    return η
end

function compute_viscosity_εij(v::AbstractCreepLaw, xx, yy, xy, args::Vararg{T, N}) where {T, N}
    εII = second_invariant(xx, yy, xy)
    η = compute_viscosity_εII(v, εII, args...)
    return η
end

function compute_viscosity_εij(v::AbstractCreepLaw, xx, yy, xy::NTuple, args::Vararg{T, N}) where {T, N}
    εII = second_invariant_staggered(xx, yy, xy)
    η = compute_viscosity_εII(v, εII, args...)
    return η
end
   
# compute effective "creep" viscosity from deviatoric stress tensor
function compute_viscosity_τII(v::AbstractCreepLaw, τII, args::Vararg{T, N}) where {T, N}
    εII = compute_εII(v, τII, args...)
    η = _viscosity(τII, εII)
    return η
end

function compute_viscosity_τij(v::AbstractCreepLaw, xx, yy, xy, args::Vararg{T, N}) where {T, N}
    τII = second_invariant(xx, yy, xy)
    η = compute_viscosity_τII(v, τII, args...)
    return η
end

function compute_viscosity_τij(v::AbstractCreepLaw, xx, yy, xy::NTuple, args::Vararg{T, N}) where {T, N}
    τII = second_invariant_staggered(xx, yy, xy)
    η = compute_viscosity_τII(v, τII, args...)
    return η
end

# compute effective "visco-elastic" viscosity
@inline compute_elastoviscosity(v::ConstantElasticity, η, dt) = (inv(η) + inv(v.G.val, dt)) |> inv
@inline compute_elastoviscosity(v::ConstantElasticity, η, args::NamedTuple) = compute_elastoviscosity(v, η, args.dt)

# compute effective "creep" viscosity from strain rate tensor given a composite rheology
@inline function compute_viscosity_εII(v::CompositeRheology, εII, args::Vararg{T, N}) where {T, N}
    e = elements(v)
    compute_viscosity_II(e, compute_viscosity_εII, εII, args...)
end

# compute effective "creep" viscosity from deviatoric stress tensor given a composite rheology
@inline function compute_viscosity_τII(v::CompositeRheology, τII, args::Vararg{T, N}) where {T, N}
    e = elements(v)
    compute_viscosity_II(e, compute_viscosity_τII, τII, args...)
end

# compute effective "creep" for a composite rheology where elements are in series
@generated function compute_viscosity_II(v::NTuple{N1, AbstractCreepLaw}, fun::F, II, args::Vararg{T, N2}) where {T, N1, N2, F}
    quote
        Base.@_inline_meta
        η = 0.0
        Base.@nexprs $N1 i -> η += inv(fun(v[i], II, args...))
        return inv(η)
    end
end