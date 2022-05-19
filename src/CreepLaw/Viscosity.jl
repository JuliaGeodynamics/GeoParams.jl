export strain_rate_circuit, 
    computeViscosity_TauII, 
    computeViscosity_EpsII, 
    computeViscosity_TauII!, 
    computeViscosity_EpsII!,
    dεII_dτII


"""
    compute viscosity given strain rate 2nd invariant

    τ = 2ηε -> η = τ/2/ε
"""
function computeViscosity_EpsII(εII, v::AbstractCreepLaw, args)
    τII = computeCreepLaw_TauII(εII, v, args) # gives 
    η = 0.5 * τII / εII
    return η
end

function computeViscosity_EpsII(εII, v::Tuple, args; tol = 1e-6)
    return local_iterations_EpsII(εII, v, args; tol = tol)
end

function local_iterations_EpsII(εII, v::Tuple, args; tol = 1e-6)
    # Initial guess
    η_ve =  computeViscosity(computeViscosity_EpsII, εII, v, args, Val(length(v))) # viscosity guess
    τII = 2 * η_ve * εII # deviatoric stress guess
    
    # Local Iterations
    iter = 0
    ϵ = 2*tol
    τII_prev = τII
    while ϵ > tol  # Newton
        iter += 1
        f = εII - strain_rate_circuit(τII, v, args)
        dfdτII =  - dεII_dτII(rand(), v, args)
        τII = τII - f / dfdτII

        ϵ = abs(τII-τII_prev)/τII
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

@inline function computeViscosity_TauII(τII::T, v::Tuple, args) where {T}
    return computeViscosity(computeViscosity_TauII, τII, v, args, Val(length(v)))
end

@generated function computeViscosity(
    fn::F, CII::T, v::Tuple, args::NamedTuple, ::Val{N}; n = -1
) where {F,T,N}
    quote
        Base.@_inline_meta
        η = 0.0
        Base.Cartesian.@nexprs $N i -> 
            η += if v[i] isa Tuple 
                computeViscosity(fn, CII, v[i], args; n=1) # viscosities in parallel → ηeff = 1/(η1 + η2)
            else
                fn(CII, v[i], args)^n # viscosities in series → ηeff = (1/η1 + 1/η2)^-1
            end
        return 1 / η
    end
end

@inline function computeViscosity_TauII!(η::AbstractArray{T,nDim}, τII::AbstractArray{T,nDim}, v::Tuple, args; cutoff=(1e16, 1e25)) where {T, nDim}
    Threads.@threads for I in eachindex(τII)
        η[I] = max(
            cutoff[1], min(cutoff[2], computeViscosity(
                    computeViscosity_TauII,
                    τII[I],
                    v,
                    (; zip(keys(args), getindex.(values(args), I))...),
                    Val(length(v)),
                )
            )
        )
    end
end

@inline function computeViscosity_EpsII!(η::AbstractArray{T,nDim}, τII::AbstractArray{T,nDim}, v::Tuple, args; cutoff=(1e16, 1e25)) where {T, nDim}
    Threads.@threads for I in eachindex(τII)
        η[I] =  max(
            cutoff[1], min(cutoff[2], computeViscosity(
                computeViscosity_EpsII,
                    τII[I],
                    v,
                    (; zip(keys(args), getindex.(values(args), I))...),
                    Val(length(v)),
                )
            )
        )
    end
end

strain_rate_circuit(TauII, v, args) = strain_rate_circuit(TauII, v, args, Val(length(v)))

@inline @generated function strain_rate_circuit(TauII, v, args, ::Val{N}; n=1) where {N}
    quote
        c = 0.0
        Base.Cartesian.@nexprs $N i ->
            c += if v[i] isa Tuple
                1 / strain_rate_circuit(TauII, v[i], args, Val(length(v[i])); n=-1)
            else
                computeCreepLaw_EpsII(TauII, v[i], args)^n
            end
        return c
    end
end

@inline function viscosityCircuit_TauII(τII::T, v, args) where {T}
    εII = strain_rate_circuit(τII, v, args)
    return η = 0.5 * τII / εII
end

@generated function dεII_dτII(τII::T, v::NTuple{N, AbstractCreepLaw}, args) where {T, N}
    quote
        Base.@_inline_meta
        val = zero(T)
        Base.Cartesian.@nexprs $N i -> val += dεII_dτII(τII, v[i], args)
        return val 
    end
end

# Plot
# using Plots
# x        = -20:.1:-11
# E2nd_vec = 10.0.^x;
# T        = (700+273.15)*ones(size(x));
# P        = 1e9   *ones(size(x));
# τII     = ones(size(x))
# args = (P=P, T=T)


# v = (SetDiffusionCreep("Dry Plagioclase | Bürgmann & Dresen (2008)"),)
 
