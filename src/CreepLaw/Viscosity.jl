export strain_rate_circuit,
    computeViscosity_TauII,
    computeViscosity_EpsII,
    computeViscosity_TauII!,
    computeViscosity_EpsII!,
    dεII_dτII

"""
    computeViscosity_EpsII(εII, v::AbstractCreepLaw, args)
    
    Compute viscosity given strain rate 2nd invariant.


"""
# method for a material with single rheology law
@inline function computeViscosity_EpsII(εII, v::AbstractCreepLaw, args)
    τII = computeCreepLaw_TauII(εII, v, args)
    η = 0.5 * τII / εII
    return η
end

"""
    computeViscosity_EpsII(εII, v::NTuple{N, AbstractCreepLaw}, args)
    
    Compute viscosity given strain rate 2nd invariant

"""
# method for a material an arbitrary rheology circuit
function computeViscosity_EpsII(
    εII, v::NTuple{N,AbstractCreepLaw}, args; tol=1e-6
) where {N}
    return local_iterations_EpsII(εII, v, args; tol=tol)
end

@inline function local_iterations_EpsII(
    εII, v::NTuple{N,AbstractCreepLaw}, args; tol=1e-6
) where {N}
    # Initial guess
    η_ve = computeViscosity(computeViscosity_EpsII, εII, v, args) # viscosity guess
    τII = 2 * η_ve * εII # deviatoric stress guess

    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    τII_prev = τII
    # while ϵ > tol
    for _ in 1:3
        iter += 1
        f = εII - strain_rate_circuit(τII, v, args)
        dfdτII = -dεII_dτII(τII, v, args)
        τII -= f / dfdτII

        ϵ = abs(τII - τII_prev) / τII
        τII_prev = τII
    end

    return 0.5 * τII / εII
end

"""
    compute viscosity given strain rate 2nd invariant

    τ = 2ηε -> η = τ/2/ε
"""
@inline function computeViscosity_TauII(τII, v, args)
    εII = computeCreepLaw_EpsII(τII, v, args)
    η = 0.5 * τII / εII
    return η
end

@inline function computeViscosity_TauII(
    τII::T, v::NTuple{N,AbstractCreepLaw}, args
) where {T,N}
    return computeViscosity(computeViscosity_TauII, τII, v, args)
end

@generated function computeViscosity(
    fn::F, CII::T, v::NTuple{N,AbstractCreepLaw}, args::NamedTuple; n=-1
) where {F,T,N}
    quote
        Base.@_inline_meta
        η = zero(T)
        Base.Cartesian.@nexprs $N i ->
            η += if v[i] isa Tuple
                computeViscosity(fn, CII, v[i], args; n=1) # viscosities in parallel → ηeff = 1/(η1 + η2)
            else
                fn(CII, v[i], args)^n # viscosities in series → ηeff = (1/η1 + 1/η2)^-1
            end
        return 1 / η
    end
end

@inline function computeViscosity_TauII!(
    η::AbstractArray{T,nDim},
    τII::AbstractArray{T,nDim},
    v::NTuple{N,AbstractCreepLaw},
    args;
    cutoff=(1e16, 1e25),
) where {T,nDim,N}
    Threads.@threads for I in eachindex(τII)
        η[I] = max(
            cutoff[1],
            min(
                cutoff[2],
                computeViscosity(
                    computeViscosity_TauII,
                    τII[I],
                    v,
                    (; zip(keys(args), getindex.(values(args), I))...),
                ),
            ),
        )
    end
end

@inline function computeViscosity_EpsII!(
    η::AbstractArray{T,nDim},
    τII::AbstractArray{T,nDim},
    v::NTuple{N,AbstractCreepLaw},
    args;
    cutoff=(1e16, 1e25),
) where {T,nDim,N}
    Threads.@threads for I in eachindex(τII)
        η[I] = max(
            cutoff[1],
            min(
                cutoff[2],
                computeViscosity(
                    computeViscosity_EpsII,
                    τII[I],
                    v,
                    (; zip(keys(args), getindex.(values(args), I))...),
                ),
            ),
        )
    end
end

@inline @generated function strain_rate_circuit(
    TauII, v::NTuple{N,AbstractCreepLaw}, args; n=1
) where {N}
    quote
        c = 0.0
        Base.Cartesian.@nexprs $N i ->
            c += if v[i] isa Tuple
                1 / strain_rate_circuit(TauII, v[i], args; n=-1)
            else
                computeCreepLaw_EpsII(TauII, v[i], args)^n
            end
        return c
    end
end

@inline function viscosityCircuit_TauII(τII, v, args)
    εII = strain_rate_circuit(τII, v, args)
    return 0.5 * τII / εII
end

@generated function dεII_dτII(τII::T, v::NTuple{N,AbstractCreepLaw}, args) where {T,N}
    quote
        Base.@_inline_meta
        val = zero(T)
        Base.Cartesian.@nexprs $N i -> val += dεII_dτII(τII, v[i], args)
        return val
    end
end
