export strain_rate_circuit,
    computeViscosity_τII,
    computeViscosity_εII,
    computeViscosity_τII!,
    computeViscosity_εII!,
    compute_τII,
    dεII_dτII,
    local_iterations_εII,
    computeViscosity

"""
    computeViscosity_EpsII(v::AbstractConstitutiveLaw, εII, args)
    
Compute viscosity given strain rate 2nd invariant for a given rheological element
"""
@inline function computeViscosity_εII(v::AbstractConstitutiveLaw, εII,  args)
    τII = compute_τII(v, εII, args)
    η = 0.5 * τII / εII
    return η
end

"""
    computeViscosity_εII(εII, v::NTuple{N, AbstractConstitutiveLaw}, args)
    
Compute viscosity given strain rate 2nd invariant
"""
function computeViscosity_εII(
    v::NTuple{N,AbstractConstitutiveLaw}, εII, args; tol=1e-6
) where {N}
        τII = local_iterations_εII(v, εII, args; tol=tol)
        η = 0.5 * τII / εII
    return  η
end


"""
    compute_τII(εII, v::NTuple{N, AbstractConstitutiveLaw}, args)
    
Compute deviatoric stress invariant given strain rate 2nd invariant

"""
function compute_τII(
    v::NTuple{N,AbstractConstitutiveLaw}, εII, args; tol=1e-6
) where {N}
        τII = local_iterations_εII(v, εII, args; tol=tol)
    return τII
end

"""

Performs local iterations versus stress
"""
@inline function local_iterations_εII(
    v::NTuple{N,AbstractConstitutiveLaw}, εII, args; tol=1e-6, verbose=false
) where {N}
    # Initial guess
    η_ve = computeViscosity(computeViscosity_εII, v, εII, args) # viscosity guess
    τII = 2 * η_ve * εII # deviatoric stress guess

    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    τII_prev = τII
    while ϵ > tol
        iter += 1
        f   = εII - strain_rate_circuit(τII, v, args)
        dfdτII = -dεII_dτII(v, τII, args)
        τII -= f / dfdτII

        ϵ = abs(τII - τII_prev) / τII
        τII_prev = τII
        if verbose
            println(" iter $(iter) $ϵ")
        end
    end
    if verbose
        println("---")
    end
    return τII
end

"""
    compute viscosity given strain rate 2nd invariant

    τ = 2ηε -> η = τ/2/ε
"""
@inline function computeViscosity_τII(v, τII, args)
    εII = compute_εII(v, τII, args)
    η = 0.5 * τII / εII
    return η
end

@inline function computeViscosity_τII(
    v::NTuple{N,AbstractConstitutiveLaw}, τII::_T, args
) where {_T,N}
    return computeViscosity(computeViscosity_τII, v, τII, args)
end

@generated function computeViscosity(
    fn::F, v::NTuple{N,AbstractConstitutiveLaw}, CII::T, args::NamedTuple; n=-1
) where {F,T,N}
    quote
        Base.@_inline_meta
        η = zero(T)
        Base.Cartesian.@nexprs $N i ->
            η += if v[i] isa Tuple
                computeViscosity(fn, v[i], CII, args; n=1) # viscosities in parallel → ηeff = 1/(η1 + η2)
            else
                fn(v[i], CII, args)^n # viscosities in series → ηeff = (1/η1 + 1/η2)^-1
            end
        return 1 / η
    end
end

@inline function computeViscosity_τII!(
    η::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    τII::AbstractArray{T,nDim},
    args;
    cutoff=(1e16, 1e25),
) where {T,nDim,N}
    Threads.@threads for I in eachindex(τII)
        η[I] = max(
            cutoff[1],
            min(
                cutoff[2],
                computeViscosity(
                    computeViscosity_τII,
                    v,
                    τII[I],
                    (; zip(keys(args), getindex.(values(args), I))...),
                ),
            ),
        )
    end
end

@inline function computeViscosity_εII!(
    η::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    τII::AbstractArray{T,nDim},
    args;
    cutoff=(1e16, 1e25),
) where {T,nDim,N}
    Threads.@threads for I in eachindex(τII)
        η[I] = max(
            cutoff[1],
            min(
                cutoff[2],
                computeViscosity(
                    computeViscosity_εII,
                    v,
                    τII[I],
                    (; zip(keys(args), getindex.(values(args), I))...),
                ),
            ),
        )
    end
end


@inline function compute_τII!(
    τII::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    εII::AbstractArray{T,nDim},
    args;
) where {T,nDim,N}
    for I in eachindex(τII)
        τII[I] =  compute_τII(  v, 
                                εII[I],
                                (; zip(keys(args), getindex.(values(args), I))...)
                             )
    end
end


@inline @generated function strain_rate_circuit(
    TauII, v::NTuple{N,AbstractConstitutiveLaw}, args; n=1
) where {N}
    quote
        c = 0.0
        Base.Cartesian.@nexprs $N i ->
            c += if v[i] isa Tuple
                1 / strain_rate_circuit(v[i], TauII, args; n=-1)
            else
                compute_εII(v[i], TauII, args)^n
            end
        return c
    end
end

@inline function viscosityCircuit_τII(v, τII, args)
    εII = strain_rate_circuit(v, τII, args)
    return 0.5 * τII / εII
end

@generated function dεII_dτII(v::NTuple{N,AbstractConstitutiveLaw}, τII::T, args) where {T,N}
    quote
        Base.@_inline_meta
        val = zero(T)
        Base.Cartesian.@nexprs $N i -> val += dεII_dτII(v[i], τII, args)
        return val
    end
end
