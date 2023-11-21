export compute_viscosity_εII,
    compute_viscosity_τII,
    compute_viscosity_II,
    compute_elastoviscosity

# extract elements from composite rheology
@inline elements(v::Union{CompositeRheology, Parallel}) = v.elements

# compute viscosity given second invariants of strain rate and deviatoric stress tensors
@inline _viscosity(τII, εII) = τII / (2  * εII)

# compute effective "creep" viscosity from strain rate tensor
"""
    compute_viscosity_εII(s::AbstractConstitutiveLaw, εII, kwargs...)

Compute effective viscosity given a 2nd invariant of the deviatoric strain rate tensor, extra parameters are passed as a named tuple, e.g., (;T=T) 
"""
@inline function compute_viscosity_εII(v::AbstractConstitutiveLaw, εII, args)

    τII = compute_τII(v, εII, args)
    η = _viscosity(τII, εII)
    return η
end

# compute effective "creep" viscosity from deviatoric stress tensor
"""
    compute_viscosity_τII(s::AbstractConstitutiveLaw, τII, kwargs...)

Compute effective viscosity given a 2nd invariant of the deviatoric stress tensor and, extra parameters are passed as a named tuple, e.g., (;T=T) 
"""
@inline function compute_viscosity_τII(v::AbstractConstitutiveLaw, τII, args)

    εII = compute_εII(v, τII, args)
    η = _viscosity(τII, εII)
    return η
end

for fn in (:compute_viscosity_εII, :compute_viscosity_τII)
    @eval begin

        @inline function $fn(v::AbstractConstitutiveLaw, xx, yy, xy, args)

            II = second_invariant(xx, yy, xy)
            η = $fn(v, II, args)
            return η
        end
                
        # For single phase versions MaterialParams
        @inline $fn(v::MaterialParams, args::Vararg{Any, N}) where {N} = $fn(v.CompositeRheology[1], args...)

       # compute effective "creep" viscosity from strain rate tensor given a composite rheology
        @inline function $fn(v::CompositeRheology, II, args::Vararg{T, N} where {T, N})
            e = elements(v)
            compute_viscosity_II(e, $fn, II, args...)
        end

    end
end


# For multi phases given the i-th phase
@generated function compute_viscosity_εII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase::Int, args::Vararg{Any, N2}) where {N1, N2}
    quote
        Base.@_inline_meta
        Base.@nexprs $N1 i -> i == phase && (return compute_viscosity_εII(v[i], args...))
        return 0.0
    end
end

@generated function compute_viscosity_τII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase::Int, args::Vararg{Any, N2}) where {N1, N2}
    quote
        Base.@_inline_meta
        Base.@nexprs $N1 i -> i == phase && (return compute_viscosity_τII(v[i], args...))
        return 0.0
    end
end
        
# For multi phases given phase ratios
@generated function compute_viscosity_εII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase_ratio::Union{NTuple{N1,T}, SVector{N1,T}}, args::Vararg{Any, N2}) where {N1, N2, T}
    quote
        Base.@_inline_meta
        val = 0.0
        Base.@nexprs $N1 i -> val += compute_viscosity_εII(v[i], args...) * phase_ratio[i]
        return val
    end
end

@generated function compute_viscosity_τII(v::NTuple{N1, AbstractMaterialParamsStruct}, phase_ratio::Union{NTuple{N1,T}, SVector{N1,T}}, args::Vararg{Any, N2}) where {N1, N2, T}
    quote
        Base.@_inline_meta
        val = 0.0
        Base.@nexprs $N1 i -> val += compute_viscosity_τII(v[i], args...) * phase_ratio[i]
        return val
    end
end

# compute effective "creep" for a composite rheology where elements are in series
@generated function compute_viscosity_II(v::NTuple{N, AbstractConstitutiveLaw}, fn::F, II::T, args)where {F, N, T}
    quote
        Base.@_inline_meta
        η = zero(T)
        Base.@nexprs $N i -> (
            v_i = v[i];
            !isplastic(v_i) && !iselastic(v_i) && (η += inv(fn(v_i, II, args)))
        )
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
        Base.@nexprs $N i -> η += fn(v[i], II, args...)
        return η
    end
end

# compute effective "visco-elastic" viscosity
@inline compute_elastoviscosity(v::ConstantElasticity, η, dt) = (inv(η) + inv(v.G.val * dt)) |> inv
@inline compute_elastoviscosity(G, η, dt) = (inv(η) + inv(G * dt)) |> inv
@inline compute_elastoviscosity(v::ConstantElasticity, η, args::NamedTuple) = compute_elastoviscosity(v, η, args.dt)
@inline compute_elastoviscosity(G, η, args::NamedTuple) = compute_elastoviscosity(G, η, args.dt)

for fn in (:compute_elastoviscosity_εII, :compute_elastoviscosity_τII)
    @eval begin
        # single phase versions
        @inline $fn(v::MaterialParams, args::Vararg{Any, N}) where {N} = $fn(v.CompositeRheology[1], args...)

        # multi-phase versions
        @generated function $fn(v::NTuple{N, AbstractMaterialParamsStruct}, phase, args) where N
            quote
                Base.@_inline_meta
                Base.@nexprs $N i -> i == phase && (return $fn(v.CompositeRheology[i], args))
                return 0.0
            end
        end

        # For multi phases given the i-th phase
        @generated function $fn(v::NTuple{N1, AbstractMaterialParamsStruct}, phase_ratio::Union{NTuple{N1,T}, SVector{N1, T}}, args::Vararg{Any, N2}) where {N1, N2, T}
            quote
                Base.@_inline_meta
                val = 0.0
                Base.@nexprs $N1 i -> val += $fn(v[i].CompositeRheology[1], args...) * phase_ratio[i]
                return 0.0
            end
        end
    end
end

# compute effective "creep" viscosity from strain rate tensor given a composite rheology
@inline function compute_elastoviscosity_εII(v::CompositeRheology, εII, args)
    e = elements(v)
    compute_elastoviscosity_II(e, compute_viscosity_εII, εII, args)
end

# compute effective "creep" viscosity from deviatoric stress tensor given a composite rheology
@inline function compute_elastoviscosity_τII(v::CompositeRheology, τII, args)
    e = elements(v)
    compute_elastoviscosity_II(e, compute_viscosity_τII, τII, args)
end

@generated function compute_elastoviscosity_II(v::NTuple{N, AbstractConstitutiveLaw}, fn::F, II, args) where {F, N}
    quote
        Base.@_inline_meta
        η = 0.0
        Base.@nexprs $N i -> !isplastic(v[i]) && (η += inv(fn(v[i], II, args)))
        return inv(η)
    end
end
