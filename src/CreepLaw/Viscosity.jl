# export strain_rate_circuit,
#     computeViscosity_τII,
#     computeViscosity_εII,
#     computeViscosity_τII!,
#     computeViscosity_εII!,
#     compute_τII,
#     compute_τII!,
#     compute_εII!,
#     compute_εII,
#     dεII_dτII,
#     local_iterations_εII,
#     computeViscosity,
#     InverseCreepLaw,
#     KelvinVoigt

struct InverseCreepLaw{N} <: AbstractConstitutiveLaw{Float64}
    v::NTuple{N,AbstractConstitutiveLaw}
    
    function InverseCreepLaw(v::Vararg{AbstractConstitutiveLaw, N}) where N
        new{N}(ntuple(i->v[i], Val(N)))
    end

    function InverseCreepLaw(v::NTuple{N,AbstractConstitutiveLaw}) where N
        new{N}(v)
    end
end


"""
    struct KelvinVoigt{N, V1, V2} <: AbstractConstitutiveLaw{Float64}
        v_el::V1
        v_vis::V2
    end

    Elastic spring and viscous dashpot in parallel. 

    τ = 2Gε + η̇ε
"""
struct KelvinVoigt{N, V1, V2} <: AbstractConstitutiveLaw{Float64}
    spring::V1
    dashpot::V2

    function KelvinVoigt(v1::AbstractConstitutiveLaw, v2::AbstractConstitutiveLaw)
        T1 = typeof(v1)
        T2 = typeof(v2)
        new{2, T1, T2}(v1, v2)
    end
end

"""
    τ = 2Gε + η̇ε =  2G ̇ε Δt + η̇ε
    ̇ε = τ/(2GΔt + η)
"""
function compute_εII(v::KelvinVoigt, τII, args)
    η = computeViscosity_τII(v.viscous, τII,  args)
    G = v.elastic.G
    εII = τII/(2*G*args.Δt + η)

    return εII
end

"""
    τ = 2Gε + η̇ε =  2G ̇ε Δt + η̇ε
"""
function compute_τII(v::KelvinVoigt, εII, args)
    η = computeViscosity_εII(v.viscous, εII,  args)
    G = v.elastic.G
    τII = εII*(2*G*args.Δt + η)

    return τII
end

"""
    computeViscosity_εII(v::AbstractConstitutiveLaw, εII, args)
    
Compute viscosity given strain rate 2nd invariant for a given rheological element
"""
@inline function computeViscosity_εII(v::AbstractConstitutiveLaw, εII,  args)
    τII = compute_τII(v, εII, args)
    η = 0.5 * τII / εII
    return η
end

# special case for linar viscous rheology
computeViscosity_εII(v::LinearViscous, args...) = v.η.val

"""
    computeViscosity_τII(v::AbstractConstitutiveLaw, τII, args; tol=1e-6, verbose=false)
    
Compute viscosity given stress 2nd invariant for a given rheological element
"""
@inline function computeViscosity_τII(v::AbstractConstitutiveLaw, τII,  args)
    εII = compute_εII(v, τII, args)
    η = 0.5 * τII / εII
    return η
end

# special case for linar viscous rheology
computeViscosity_τII(v::LinearViscous, args...) = v.η.val

"""
    computeViscosity_εII(v::NTuple{N, AbstractConstitutiveLaw}, εII, args)
    
Compute viscosity given strain rate 2nd invariant
"""
function computeViscosity_εII(
    v::NTuple{N,AbstractConstitutiveLaw}, εII, args; tol=1e-6, verbose=false
) where {N}
        τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose)
        η = 0.5 * τII / εII
    return  η
end

"""
    computeViscosity_τII(v::NTuple{N, AbstractConstitutiveLaw}, τII, args; tol=1e-6, verbose=false)
    
Compute viscosity given deviatoric stress 2nd invariant
"""
function computeViscosity_τII(
    v::NTuple{N,AbstractConstitutiveLaw}, τII, args; tol=1e-6, verbose=false
) where {N}
        εII = local_iterations_τII(v, τII, args; tol=tol, verbose=verbose)
        η = 0.5 * τII / εII
    return  η
end

"""
    compute_τII(v::NTuple{N, AbstractConstitutiveLaw}, εII, args; tol=1e-6, verbose=false)
    
Compute deviatoric stress invariant given strain rate 2nd invariant

"""
function compute_τII(
    v::NTuple{N,AbstractConstitutiveLaw}, εII, args; tol=1e-6, verbose=false
) where {N}

       τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose)
    return τII
end



"""
    compute_εII(v::NTuple{N, AbstractConstitutiveLaw},τII,  args)
    
Compute deviatoric strain rate given deviatoric stress invariant

"""
function compute_εII(
    v::NTuple{N,AbstractConstitutiveLaw}, τII, args; tol=1e-6, verbose=false, n=1
) where {N}
        εII = local_iterations_τII(v, τII, args; tol=tol, verbose=verbose, n=n)
    return εII
end


"""

Performs local iterations versus stress for a given strain rate 
"""
@inline function local_iterations_εII(
    v::NTuple{N,AbstractConstitutiveLaw}, εII, args; tol=1e-6, verbose=false
) where {N}
    # Initial guess
    η_ve = computeViscosity(computeViscosity_εII, v, εII, args) # viscosity guess
    @show η_ve
    τII = 2 * η_ve * εII # deviatoric stress guess

    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    τII_prev = τII
    while ϵ > tol
        iter += 1
        f = εII - strain_rate_circuit(v, τII, args)
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

Performs local iterations versus strain rate for a given stress
"""
@inline function local_iterations_τII(
    v::NTuple{N,AbstractConstitutiveLaw}, τII, args; tol=1e-6, verbose=false, n=1
) where {N}
    # Initial guess
    η_ve = computeViscosity(computeViscosity_τII, v, τII, args) # viscosity guess
  
    εII = τII / (2 * η_ve )  # deviatoric strain rate guess
    @show η_ve, εII

    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    εII_prev = εII
    while ϵ > tol
        iter += 1
        f = τII - stress_circuit(v, εII, args, n=n)
        dfdεII = -dτII_dεII(v, εII, args)
        εII -= f / dfdεII

        ϵ = abs(εII - εII_prev) / εII
        εII_prev = εII
        if verbose
            println(" iter $(iter) $ϵ")
        end
    end
    if verbose
        println("---")
    end
    return εII
end

"""
    compute viscosity given strain rate 2nd invariant

    τ = 2ηε -> η = τ/2/ε
"""
@inline function computeViscosity_τII(v, τII, args; tol=1e-6, verbose=false, n=1)
    εII = compute_εII(v, τII, args; tol=tol, verbose=verbose)
    η = 0.5 * τII / εII
    return η
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
    εII::AbstractArray{T,nDim},
    args;
    cutoff=(1e16, 1e25),
) where {T,nDim,N}
    Threads.@threads for I in eachindex(εII)
        η[I] = max(
            cutoff[1],
            min(
                cutoff[2],
                computeViscosity(
                    computeViscosity_εII,
                    v,
                    εII[I],
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


@inline function compute_εII!(
    εII::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    τII::AbstractArray{T,nDim},
    args;
) where {T,nDim,N}
    for I in eachindex(εII)
        εII[I] =  compute_εII(  v, 
                                τII[I],
                                (; zip(keys(args), getindex.(values(args), I))...)
                             )
    end
end


@inline @generated function strain_rate_circuit(
    v::NTuple{N,AbstractConstitutiveLaw},  TauII, args
) where {N}
    quote
        c = 0.0
        Base.Cartesian.@nexprs $N i ->
            c += if v[i] isa InverseCreepLaw
                strain_rate_circuit(v[i], TauII, args)
            else
                compute_εII(v[i], TauII, args)
            end
        return c
    end
end

@generated function strain_rate_circuit(v_ice::InverseCreepLaw{N}, τII::_T, args) where {_T,N}
    quote
        Base.@_inline_meta
        c = zero(_T)
        Base.Cartesian.@nexprs $N i -> c += 1/compute_εII(v_ice.v[i], τII, args)
        return 1/c
    end
end

@inline @generated function stress_circuit(
 v::NTuple{N,AbstractConstitutiveLaw}, EpsII, args; n=1
) where {N}
    quote
        c = 0.0
        Base.Cartesian.@nexprs $N i ->
            c += if v[i] isa Tuple
                1 / stress_circuit(v[i], TauII, args; n=-1)
            else
                compute_τII(v[i], EpsII, args)^n
            end
        return c
    end
end

@inline function viscosityCircuit_τII(v, τII, args)
    εII = strain_rate_circuit(v, τII, args)
    return 0.5 * τII / εII
end

@inline function viscosityCircuit_εII(v, εII, args)
    τII = stress_circuit(v, εII, args)
    return 0.5 * τII / εII
end

@generated function dεII_dτII(v::NTuple{N,AbstractConstitutiveLaw}, τII::_T, args) where {_T,N}
    quote
        Base.@_inline_meta
        val = zero(_T)
        Base.Cartesian.@nexprs $N i -> val += dεII_dτII(v[i], τII, args)
        return val
    end
end

@generated function dτII_dεII(v::NTuple{N,AbstractConstitutiveLaw}, εII::_T, args) where {_T,N}
    quote
        Base.@_inline_meta
        val = zero(_T)
        Base.Cartesian.@nexprs $N i -> val += dτII_dεII(v[i], εII, args)
        return val
    end
end