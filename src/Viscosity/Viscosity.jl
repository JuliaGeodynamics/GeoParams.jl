# extract elements from composite rheology
@inline elements(v::Union{CompositeRheology, Parallel}) = v.elements

# compute effective "creep" viscosity from strain rate tensor
"""
    compute_viscosity_εII(s::AbstractConstitutiveLaw, εII, kwargs...)

Compute effective viscosity given a 2nd invariant of the deviatoric strain rate tensor, extra parameters are passed as a named tuple, e.g., (;T=T)
"""
@inline compute_viscosity_εII(v::AbstractConstitutiveLaw, εII, args) = compute_τII(v, εII, args) / (2 * εII)
@inline compute_viscosity_εII(v::LinearViscous, εII, args) = v.η.val
@inline compute_viscosity_εII(v::ConstantElasticity, εII, args) = v.G * args.dt
@inline compute_viscosity_εII(v::HerschelBulkley, εII, args) = compute_hb_viscosity_εII(v, εII; args...)


# compute effective "creep" viscosity from deviatoric stress tensor
"""
    compute_viscosity_τII(s::AbstractConstitutiveLaw, τII, kwargs...)

Compute effective viscosity given a 2nd invariant of the deviatoric stress tensor and, extra parameters are passed as a named tuple, e.g., (;T=T)
"""
@inline compute_viscosity_τII(v::AbstractConstitutiveLaw, τII, args) = τII / (2 * compute_εII(v, τII, args))
@inline compute_viscosity_τII(v::LinearViscous, τII, args) = v.η.val
@inline compute_viscosity_τII(v::ConstantElasticity, τII, args) = v.G * args.dt
@inline compute_viscosity_τII(v::HerschelBulkley, εII, args) = compute_hb_viscosity_τII(v, εII; args...)

for fn in (:compute_viscosity_εII, :compute_viscosity_τII)
    @eval begin

        @inline $fn(v::AbstractConstitutiveLaw, xx, yy, xy, args) = $fn(v, second_invariant(xx, yy, xy), args)

        # For single phase versions MaterialParams
        @inline $fn(v::MaterialParams, args::Vararg{Any, N}) where {N} = $fn(v.CompositeRheology[1], args...)

        # compute effective "creep" viscosity from strain rate tensor given a composite rheology
        @inline $fn(v::CompositeRheology, II, args...) = compute_viscosity_II(elements(v), $fn, II, args...)

    end
end


# For multi phases given the i-th phase
@generated function compute_viscosity_εII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase::Int, args::Vararg{Any, N2}) where {N1, N2}
    return quote
        Base.@_inline_meta
        Base.@nexprs $N1 i -> i == phase && (return compute_viscosity_εII(v[i], args...))
        return 0.0
    end
end

@generated function compute_viscosity_τII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase::Int, args::Vararg{Any, N2}) where {N1, N2}
    return quote
        Base.@_inline_meta
        Base.@nexprs $N1 i -> i == phase && (return compute_viscosity_τII(v[i], args...))
        return 0.0
    end
end

# For multi phases given phase ratios
@generated function compute_viscosity_εII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase_ratio::Union{NTuple{N1}, SVector{N1}}, args::Vararg{Any, N2}) where {N1, N2}
    return quote
        Base.@_inline_meta
        val = 0.0
        Base.@nexprs $N1 i -> val += compute_viscosity_εII(v[i], args...) * phase_ratio[i]
        return val
    end
end

@generated function compute_viscosity_τII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase_ratio::Union{NTuple{N1}, SVector{N1}}, args::Vararg{Any, N2}) where {N1, N2}
    return quote
        Base.@_inline_meta
        val = 0.0
        Base.@nexprs $N1 i -> val += compute_viscosity_τII(v[i], args...) * phase_ratio[i]
        return val
    end
end

# compute effective "creep" for a composite rheology where elements are in series
@generated function compute_viscosity_II(v::NTuple{N, AbstractConstitutiveLaw}, fn::F, II::T, args) where {F, N, T}
    return quote
        Base.@_inline_meta
        η = zero(T)
        Base.@nexprs $N i -> (
            v_i = v[i];
            !isplastic(v_i) && !iselastic(v_i) && (η += inv(fn(v_i, II, args)))
        )
        return inv(η)
    end
end

# compute effective "creep" viscosity from strain rate/stress tensor given a parallel rheology
@inline compute_viscosity_εII(v::Parallel, εII, args) = compute_viscosity_II_parallel(elements(v), compute_viscosity_εII, εII, args)
@inline compute_viscosity_τII(v::Parallel, τII, args) = compute_viscosity_II_parallel(elements(v), compute_viscosity_τII, τII, args)

# compute effective "creep" for a composite rheology where elements are in parallel
@generated function compute_viscosity_II_parallel(v::NTuple{N, AbstractConstitutiveLaw}, fn::F, II, args) where {F, N}
    return quote
        Base.@_inline_meta
        η = 0.0
        Base.@nexprs $N i -> η += fn(v[i], II, args)
        return η
    end
end

# compute effective "visco-elastic" viscosity
@inline compute_elastoviscosity(v::ConstantElasticity, η, dt) = compute_elastoviscosity(v.G, η, dt)
@inline compute_elastoviscosity(G, η, dt) = (inv(η) + inv(G * dt)) |> inv
@inline compute_elastoviscosity(v::ConstantElasticity, η, args::NamedTuple) = compute_elastoviscosity(v.G, η, args.dt)
@inline compute_elastoviscosity(G, η, args::NamedTuple) = compute_elastoviscosity(G, η, args.dt)

for fn in (:compute_elastoviscosity_εII, :compute_elastoviscosity_τII)
    @eval begin
        # single phase versions
        @inline $fn(v::MaterialParams, args::Vararg{Any, N}) where {N} = $fn(v.CompositeRheology[1], args...)

        # multi-phase versions
        @generated function $fn(v::NTuple{N, AbstractMaterialParamsStruct}, phase::Int, args::Vararg{Any, N2}) where {N, N2}
            return quote
                Base.@_inline_meta
                Base.@nexprs $N i -> i == phase && (return $$fn(v[i].CompositeRheology[1], args...))
                return 0.0
            end
        end

        # For multi phases given phase ratios
        @generated function $fn(v::NTuple{N1, AbstractMaterialParamsStruct}, phase_ratio::Union{NTuple{N1, T}, SVector{N1, T}}, args::Vararg{Any, N2}) where {N1, N2, T}
            return quote
                Base.@_inline_meta
                val = 0.0
                Base.@nexprs $N1 i -> val += $$fn(v[i].CompositeRheology[1], args...) * phase_ratio[i]
                return val
            end
        end
    end
end

@inline compute_elastoviscosity_εII(v::CompositeRheology, εII, args) = compute_elastoviscosity_II(elements(v), compute_viscosity_εII, εII, args)
@inline compute_elastoviscosity_τII(v::CompositeRheology, τII, args) = compute_elastoviscosity_II(elements(v), compute_viscosity_τII, τII, args)

@generated function compute_elastoviscosity_II(v::NTuple{N, AbstractConstitutiveLaw}, fn::F, II, args) where {F, N}
    return quote
        Base.@_inline_meta
        η = 0.0
        Base.@nexprs $N i -> !isplastic(v[i]) && (η += inv(fn(v[i], II, args)))
        return inv(η)
    end
end

# special cases for constant viscosity and elasticity

@inline compute_viscosity(v::LinearViscous; kwargs...) = v.η.val
@inline compute_viscosity(v::ConstantElasticity; dt = 0.0, kwargs...) = v.G * dt
@inline compute_viscosity(v::Union{LinearViscous, ConstantElasticity}, kwargs) = compute_viscosity(v; kwargs...)
@inline compute_viscosity(v, kwargs) = throw("compute_viscosity only works for linear rheologies")

@inline compute_viscosity(v::CompositeRheology, args) = compute_viscosity(elements(v), args)

@generated function compute_viscosity(v::NTuple{N, AbstractConstitutiveLaw}, args) where {N}
    return quote
        Base.@_inline_meta
        η = 0.0
        Base.@nexprs $N i -> !isplastic(v[i]) && (η += inv(compute_viscosity(v[i], args)))
        return inv(η)
    end
end

# single phase versions
@inline compute_viscosity(v::MaterialParams, args::Vararg{Any, N}) where {N} = compute_viscosity(v.CompositeRheology[1], args...)

# multi-phase versions
@generated function compute_viscosity(v::NTuple{N, AbstractMaterialParamsStruct}, phase, args) where {N}
    return quote
        Base.@_inline_meta
        Base.@nexprs $N i -> i == phase && (return compute_viscosity(v[i].CompositeRheology[1], args))
        return 0.0
    end
end

# For multi phases given phase ratios
@generated function compute_viscosity(v::NTuple{N1, AbstractMaterialParamsStruct}, phase_ratio::Union{NTuple{N1, T}, SVector{N1, T}}, args::Vararg{Any, N2}) where {N1, N2, T}
    return quote
        Base.@_inline_meta
        val = 0.0
        Base.@nexprs $N1 i -> val += compute_viscosity(v[i].CompositeRheology[1], args...) * phase_ratio[i]
        return val
    end
end

@inline compute_elasticviscosity(v::CompositeRheology, args) = compute_elasticviscosity(elements(v), args)

@generated function compute_elasticviscosity(v::NTuple{N, AbstractMaterialParamsStruct}, phase, args) where {N}
    return quote
        Base.@_inline_meta
        Base.@nexprs $N i -> i == phase && (return compute_elasticviscosity(v[i].CompositeRheology[1], args))
        return 0.0
    end
end

@generated function compute_elasticviscosity(v::NTuple{N, AbstractConstitutiveLaw}, args) where {N}
    return quote
        Base.@_inline_meta
        Base.@nexprs $N i -> iselastic(v[i]) && return compute_viscosity(v[i], args)
        return 0.0
    end
end

@generated function compute_elasticviscosity(v::NTuple{N1, AbstractMaterialParamsStruct}, phase_ratio::Union{NTuple{N1}, SVector{N1}}, args::Vararg{Any, N2}) where {N1, N2}
    return quote
        Base.@_inline_meta
        val = 0.0
        Base.@nexprs $N1 i -> val += compute_elasticviscosity(v[i].CompositeRheology[1], args...) * phase_ratio[i]
        return val
    end
end
