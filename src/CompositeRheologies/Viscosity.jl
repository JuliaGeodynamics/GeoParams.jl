@inline _viscosity(τII, εII) = τII / (2  * εII)

# compute effective "creep" viscosity from strain rate tensor
function compute_viscosity_εII(v::AbstractCreepLaw, εII, args::Vararg{T, N}) where {T, N}
    τII = compute_τII(v, εII, args)
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
    εII = compute_εII(v, εII, args)
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
@inline compute_elastoviscosity_τII(v::ConstantElasticity, η, dt) = inv(inv(η) + inv(v.G.val, dt))
@inline compute_elastoviscosity_τII(v::ConstantElasticity, η, args::NamedTuple) = compute_elastoviscosity_τII(v, η, args.dt)