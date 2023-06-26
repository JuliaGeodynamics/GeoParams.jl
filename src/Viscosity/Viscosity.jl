export compute_viscosity_εII,
    compute_viscosity_τII,
    compute_elastoviscosity

# extract elements from composite rheology
@inline elements(v::Union{CompositeRheology, Parallel}) = v.elements

# compute viscosity given second invariants of strain rate and deviatoric stress tensors
@inline _viscosity(τII, εII) = τII / (2  * εII)

# compute effective "creep" viscosity from strain rate tensor
function compute_viscosity_εII(v::AbstractConstitutiveLaw, εII, args)
    τII = compute_τII(v, εII, args)
    η = _viscosity(τII, εII)
    return η
end

function compute_viscosity_εII(v::AbstractConstitutiveLaw, xx, yy, xy, args)
    εII = second_invariant(xx, yy, xy)
    η = compute_viscosity_εII(v, εII, args)
    return η
end
   
# compute effective "creep" viscosity from deviatoric stress tensor
function compute_viscosity_τII(v::AbstractConstitutiveLaw, τII, args)
    εII = compute_εII(v, τII, args)
    η = _viscosity(τII, εII)
    return η
end

function compute_viscosity_τII(v::AbstractConstitutiveLaw, xx, yy, xy, args)
    τII = second_invariant(xx, yy, xy)
    η = compute_viscosity_τII(v, τII, args)
    return η
end

# single phase versions
compute_viscosity_εII(v::MaterialParams, args::Vararg{Any, N}) where {N} = compute_viscosity_εII(v.CompositeRheology[1], args...)
compute_viscosity_τII(v::MaterialParams, args::Vararg{Any, N}) where {N} = compute_viscosity_τII(v.CompositeRheology[1], args...)

# multi-phase versions
@generated function compute_viscosity_εII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase, args::Vararg{Any, N2}) where {N1, N2}
    quote
        Base.@nexprs $N1 i -> i == phase && (return compute_viscosity_εII(v.CompositeRheology[i], args...))
        return 0.0
    end
end

@generated function compute_viscosity_τII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase, args::Vararg{Any, N2}) where {N1, N2}
    quote
        Base.@nexprs $N1 i -> i == phase && (return compute_viscosity_τII(v.CompositeRheology[i], args...))
        return 0.0
    end
end

# compute effective "creep" viscosity from strain rate tensor given a composite rheology
@inline function compute_viscosity_εII(v::CompositeRheology, εII, args::Vararg{T, N} where {T, N})
    e = elements(v)
    compute_viscosity_II(e, compute_viscosity_εII, εII, args...)
end

# compute effective "creep" viscosity from deviatoric stress tensor given a composite rheology
@inline function compute_viscosity_τII(v::CompositeRheology, τII, args::Vararg{T, N} where {T, N})
    e = elements(v)
    compute_viscosity_II(e, compute_viscosity_τII, τII, args...)
end

# compute effective "creep" for a composite rheology where elements are in series
@generated function compute_viscosity_II(v::NTuple{N, AbstractConstitutiveLaw}, fn::F, II, args) where {F, N}
    quote
        Base.@_inline_meta
        η = 0.0
        Base.@nexprs $N i -> η += inv(fn(v[i], II, args)) * !iselastic(v[i]) * !isplastic(v[i])
        return inv(η)
    end
end

# compute effective "creep" viscosity from strain rate tensor given a composite rheology
@inline function compute_viscosity_εII(v::Parallel, εII, args)
    e = elements(v)
    compute_viscosity_II_parallel(e, compute_viscosity_εII, εII, args)
end

# compute effective "creep" viscosity from deviatoric stress tensor given a composite rheology
@inline function compute_viscosity_τII(v::Parallel, τII, args)
    e = elements(v)
    compute_viscosity_II_parallel(e, compute_viscosity_τII, τII, args)
end

# compute effective "creep" for a composite rheology where elements are in parallel
@generated function compute_viscosity_II_parallel(v::NTuple{N, AbstractConstitutiveLaw}, fn::F, II, args) where {F, N}
    quote
        Base.@_inline_meta
        η = 0.0
        Base.@nexprs $N i -> η += fn(v[i], II, args)
        return η
    end
end

# compute effective "visco-elastic" viscosity
@inline compute_elastoviscosity(v::ConstantElasticity, η, dt) = (inv(η) + inv(v.G.val * dt)) |> inv
@inline compute_elastoviscosity(G, η, dt) = (inv(η) + inv(G * dt)) |> inv
@inline compute_elastoviscosity(v::ConstantElasticity, η, args::NamedTuple) = compute_elastoviscosity(v, η, args.dt)
@inline compute_elastoviscosity(G, η, args::NamedTuple) = compute_elastoviscosity(G, η, args.dt)