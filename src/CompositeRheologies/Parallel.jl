"""
Put rheological elements in parallel 
"""
struct Parallel{T, N, Nplast, is_plastic, Nvol, is_vol} <: AbstractConstitutiveLaw{T}
    elements::T
end

function Parallel(v::T) where {T}
    v = tuple(v...)
    n = length(v)

    is_plastic = isa.(v, AbstractPlasticity)     # Is one of the elements a plastic element?
    Nplast = count(is_plastic)

    is_vol = isvolumetric.(v)
    Nvol = count(is_vol)

    return Parallel{typeof(v), n, Nplast, is_plastic, Nvol, is_vol}(v)
end
Parallel(a, b...) = Parallel((a, b...))

@generated function getindex(p::Parallel{T, N}, I::Int64) where {T, N}
    return quote
        Base.@_inline_meta
        @assert I ≤ $N
        Base.Cartesian.@nexprs $N i -> I == i && return p.elements[i]
    end
end


function show(io::IO, a::Parallel)
    println(io, "Parallel:   ")

    # Compose a string with rheological elements, so we have an overview in the REPL
    str = print_rheology_matrix(a)
    println.(str)

    return nothing
end

isplastic(v::Parallel{T, N, 0, is_plastic}) where {T, N, is_plastic} = false;
isplastic(v::Parallel{T, N, Nplast, is_plastic}) where {T, N, Nplast, is_plastic} = true;
isvolumetric(v::Parallel{T, N, Nplast, is_plastic, 0, is_vol}) where {T, N, Nplast, is_plastic, is_vol} = false;
isvolumetric(v::Parallel{T, N, Nplast, is_plastic, Nvol, is_vol}) where {T, N, Nplast, is_plastic, Nvol, is_vol} = true;

# COMPUTE STRAIN RATE
"""
    compute_εII(v::Parallel{T,N}, τII, args; tol=1e-6, verbose=false, n=1)

Computing `εII` as a function of `τII` for a Parallel elements is (usually) a nonlinear problem
"""
function compute_εII(
        v::Parallel{T, N},
        τII::_T,
        args;
        tol = 1.0e-6, verbose = false, n = 1
    ) where {T, N, _T}
    εII = local_iterations_τII(v, τII, args; tol = tol, verbose = verbose, n = n)
    return εII
end

# Here we do need to do iterations
function compute_εII_AD(v::Parallel, τII, args; tol = 1.0e-6, verbose = false)
    return local_iterations_τII_AD(v, τII, args; tol = tol, verbose = verbose)
end

# Sum τII of parallel elements (with & w/out quantities)
compute_τII(v::Parallel{T, N}, εII::_T, args; tol = 1.0e-6, verbose = false) where {T, _T, N} = nreduce(vi -> first(compute_τII(vi, εII, args)), v.elements)
compute_τII(v::Parallel{T, N}, εII::Quantity, args; tol = 1.0e-6, verbose = false) where {T, N} = nreduce(vi -> first(compute_τII(vi, εII, args)), v.elements)

compute_τII_AD(v::Parallel{T, N}, εII::_T, args; tol = 1.0e-6, verbose = false) where {T, N, _T} = first(compute_τII(v, εII, args))

# sum P for parallel elements:
compute_p(v::Parallel{T, N}, εvol::_T, args; tol = 1.0e-6, verbose = false) where {T, _T, N} = nreduce(vi -> first(compute_p(vi, εvol, args)), v.elements)


@generated function ∂Q∂τII(
        v::Parallel{T, N, Nplast, is_plastic}, τ::_T, args
    ) where {_T, T, N, Nplast, is_plastic}
    return quote
        Base.@_inline_meta
        Base.Cartesian.@nexprs $N i -> is_plastic[i] == true && return ∂Q∂τII(v[i], τ, args)
    end
end

@generated function ∂Q∂P(
        v::Parallel{T, N, Nplast, is_plastic}, P::_T, args
    ) where {_T, T, N, Nplast, is_plastic}
    return quote
        Base.@_inline_meta
        Base.Cartesian.@nexprs $N i -> is_plastic[i] == true && return ∂Q∂P(v[i], P, args)
    end
end

@generated function compute_yieldfunction(
        v::Parallel{T, N, Nplast, is_plastic}, args
    ) where {T, N, Nplast, is_plastic}
    return quote
        Base.@_inline_meta
        Base.Cartesian.@nexprs $N i -> is_plastic[i] == true && return compute_yieldfunction(v[i], args)
    end
end

compute_yieldfunction(v::Parallel{T, N, 0, is_plastic}, args) where {T, N, is_plastic} = NaN

function dεII_dτII(
        v::Parallel{T, N}, τII::_T, args
    ) where {T, N, _T}
    ε = compute_εII(v, τII, args)
    return inv(dτII_dεII(v, ε, args))
end

"""
    dτII_dεII(v::Parallel{T,N}, TauII::_T, args)

Computes the derivative of `τII` vs `εII` for parallel elements   
"""
dτII_dεII(v::Parallel{T, N}, τII::_T, args) where {T, _T, N} = nreduce(vi -> first(dτII_dεII(vi, τII, args)), v.elements)

"""
    dτII_dεII_nonplastic(v::Parallel{T,N}, TauII::_T, args)

Computes the derivative of `τII` vs `εII` for parallel elements that are non-plastic  
"""
dτII_dεII_nonplastic(v::Parallel{T, N}, τII::_T, args) where {T, _T, N} = nreduce(vi -> first(dτII_dεII_nonplastic(vi, τII, args)), v.elements)
dτII_dεII_nonplastic(v, TauII, args) = dτII_dεII(v, TauII, args)
dτII_dεII_nonplastic(v::AbstractPlasticity, TauII, args) = 0.0

"""
    compute_τII_nonplastic(v::Parallel, EpsII, args)

Harmonic average of stress of all elements in a `CompositeRheology` structure that are not || elements
"""
compute_τII_nonplastic(v::Parallel{T, N}, τII::_T, args) where {T, _T, N} = nreduce(vi -> first(_compute_τII_nonplastic(vi, τII, args)), v.elements)
_compute_τII_nonplastic(v, EpsII, args) = first(compute_τII(v, EpsII, args))
_compute_τII_nonplastic(v::AbstractPlasticity, EpsII, args) = 0.0
