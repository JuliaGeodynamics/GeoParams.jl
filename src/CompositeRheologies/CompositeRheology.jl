# All related to the CompositeRheology struct 

"""
    Structure that holds composite rheologies (e.g., visco-elasto-viscoplastic),
    but also indicates (in the name) whether we need to perform non-linear iterations.
"""
struct CompositeRheology{T, N, 
                        Npar, is_parallel, 
                        Nplast, is_plastic, 
                        Nvol, is_vol,
                        vol_plastic
                        } <: AbstractComposite
    elements::T
end

# Defines tuples of composite rheologies, while also checking which type of iterations need to be performed
function CompositeRheology(v::T) where {T}

    # determine if we have parallel elements & if yes: where
    n = length(v)
    is_parallel = isa.(v,Parallel)
    Npar = count(is_parallel)

    # determine if we have plastic elements 
    is_plastic = isplastic.(v) 
    Nplast = count(is_plastic)

    # determine if we have elements that have volumetric deformation
    is_vol = isvolumetric.(v);
    Nvol   =   count(is_vol);

    # determine if we have a volumetric plastic element
    vol_plastic = any(is_plastic .& is_vol)
    
    return CompositeRheology{typeof(v), n, Npar, is_parallel, Nplast, is_plastic, Nvol, is_vol, vol_plastic}(v)
end
CompositeRheology(a,b...) = CompositeRheology( (a,b...,)) 
CompositeRheology(a::Parallel) = CompositeRheology( (a,)) 

@generated function getindex(p::CompositeRheology{T, N}, I::Int64) where {T,N}
    quote
        Base.@_inline_meta
        @assert I ≤ $N
        Base.Cartesian.@nexprs $N i -> I == i && return p.elements[i]
        return 0
    end
end

# Print info in the REPL
include("CompositeRheologies_print.jl")

function show(io::IO, g::AbstractComposite)

    # Compose a string with rheological elements, so we have an overview in the REPL
    str = print_rheology_matrix(g)
    println.(str)
    
    return nothing
end


# HELPER FUNCTIONS

# determine if 3 element is plastic or not
@inline isplastic(v) = false;
@inline isplastic(v::AbstractPlasticity) = true;
@inline isplastic(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic}) where {T, N,  Npar, is_parallel, Nplast, is_plastic} = true;
@inline isplastic(v::CompositeRheology{T, N,  Npar, is_parallel, 0, is_plastic}) where {T, N,  Npar, is_parallel, is_plastic} = false;

@inline isvolumetric(v) = false;
@inline isvolumetric(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic, 0, is_vol}) where {T,N,Npar, is_parallel, Nplast,is_plastic, is_vol} = false;
@inline isvolumetric(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic, Nvol, is_vol}) where {T,N,Npar, is_parallel, Nplast,is_plastic, Nvol, is_vol} = true;

@inline isvolumetricplastic(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic, Nvol, is_vol, volumetricplasticity}) where {T,N,Npar, is_parallel, Nplast,is_plastic, Nvol, is_vol, volumetricplasticity} = volumetricplasticity;

"""
    compute_εII(v::CompositeRheology{T,N}, τII, args; tol=1e-6, verbose=false, n=1)

Computing `εII` as a function of `τII` for a composite element is the sum of the individual contributions
"""
@inline compute_εII(v::CompositeRheology, τII, args; tol=1e-6, verbose=false) = nreduce(vi -> first(compute_εII(vi, τII, args)), v.elements)
@inline compute_εII(v::CompositeRheology, τII::Quantity, args; tol=1e-6, verbose=false) = nreduce(vi -> first(compute_εII(vi, τII, args)), v.elements)
@inline compute_εII(v::AbstractMaterialParamsStruct, τII, args) = compute_εII(v.CompositeRheology[1], τII, args)

# compute strain rate partitioning
@inline compute_elements_εII(v::CompositeRheology{T,N}, τII, args) where {T, N} = ntuple(i -> first(compute_εII( v.elements[i], τII, args)), Val(N));

"""
    ε_part = compute_elements_εII(v::NTuple{N1,AbstractMaterialParamsStruct}, τII::NTuple{N2,T}, args, phase::I)
This returns the individual strainrate components of the `CompositeRheology` defined in `v`
"""
function compute_elements_εII(
    v::NTuple{N,AbstractMaterialParamsStruct},
    τII,
    args,
    phase::Int,
) where N

    ε_part = nphase(vi -> compute_elements_εII(vi.CompositeRheology[1], τII, args), phase, v)

    return ε_part
end

# As we don't do iterations, this is the same
function compute_εII_AD(v::CompositeRheology, τII, args; tol=1e-6, verbose=false)
    return  compute_εII(v, τII, args)
end
compute_εII_AD(v::AbstractMaterialParamsStruct, τII, args) = compute_εII_AD(v.CompositeRheology[1], τII, args)

#COMPUTE VOLUMETRIC STRAIN-RATE
"""
    compute_εvol(v::CompositeRheology{T,N}, p, args; tol=1e-6, verbose=false, max_iter=1000)

Computing `εvol` as a function of `p` for a composite element is the sum of the individual contributions
"""
compute_εvol(v::CompositeRheology, p, args; tol=1e-6, verbose=false) = nreduce(vi -> first(compute_εvol(vi, p, args)), v.elements)
compute_εvol(v::CompositeRheology, p::Quantity, args; tol=1e-6, verbose=false) = nreduce(vi -> first(compute_εvol(vi, p, args)), v.elements)

# COMPUTE DEVIATORIC STRESS AND PRESSURE
function compute_τII(v::CompositeRheology{T,N,0}, εII, args; tol=1e-6, verbose=false, max_iter=1000, AD=false) where {T,N}
    # A composite rheology case with no parallel element; iterations for τII
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose, max_iter=max_iter, AD=AD)
    return τII
end

"""
    p,τII = compute_p_τII(v::CompositeRheology, εII, εvol, args; tol=1e-6, verbose=false, max_iter=1000, AD=false) 

This updates pressure `p` and deviatoric stress invariant `τII` in case the composite rheology structure has volumetric components, but does not contain plastic or parallel elements.
The 'old' pressure should be stored in `args` as `args.P_old`   
"""
function compute_p_τII(
        v::CompositeRheology{T,N,
                    Npar,is_parallel,
                    Nplastic,is_plastic,
                    Nvol,is_vol,
                    false}, 
        εII, 
        εvol,
        args; 
        tol=1e-6, verbose=false, max_iter=1000, AD=false
    ) where {T, N, Npar, is_parallel, Nplastic, is_plastic, Nvol, is_vol}

    # A composite rheology case that may have volumetric elements, but the are not 
    # tightly coupled, so we do NOT perform coupled iterations.
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose, AD=false)
    P   = local_iterations_εvol(v, εvol, args; tol=tol, verbose=verbose, AD=false)

    return P,τII
end


"""
    p,τII = compute_p_τII(v::CompositeRheology, εII, εvol, args; tol=1e-6, verbose=false, AD=false) 

This updates pressure `p` and deviatoric stress invariant `τII` in case the composite rheology structure has no volumetric elements, but may contain plastic or parallel elements. 
In that case, pressure is not updated (`args.P` is used instead).     
"""
function compute_p_τII(
        v::CompositeRheology{T,N,
                    Npar,is_parallel,
                    Nplast,is_plastic,
                    0,is_vol, 
                    false}, 
        εII,
        εvol,
        args;
        tol=1e-6, verbose=false, AD=false
    ) where {T, N, Npar, is_parallel, Nplast, is_plastic, is_vol}
    
    # A composite rheology case with no parallel element; iterations for τII
    τII, = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose)
    
    P = any(keys(args) .=== :P_old) ? args.P_old : 0.0

    return P, τII
end

"""
    p,τII = compute_p_τII(v::CompositeRheology, εII, εvol, args; tol=1e-6, verbose=false) 

This updates pressure `p` and deviatoric stress invariant `τII` in case the composite rheology structure has volumetric components and has volumetric plasticity
The 'old' pressure should be stored in `args` as `args.P_old`   
"""
function compute_p_τII(
        v::CompositeRheology{T,N,
                        Npar,is_parallel,
                        Nplastic,is_plastic,
                        Nvol,is_vol,
                        true}, 
        εII, 
        εvol,
        args; 
        tol=1e-6, verbose=false
    ) where {T, N, Npar, is_parallel, Nplastic, is_plastic, Nvol, is_vol}

    # A composite rheology case that may have volumetric elements, but the are not 
    # tightly coupled, so we do NOT perform coupled iterations.
    out = local_iterations_εvol_εII(v, εII, εvol, args; tol=tol, verbose=verbose)

    τII = out[1]
    P = out[end]

    return P,τII, out[2:N]
end

# COMPUTE STRAIN RATE
"""
    τII = compute_τII(v::CompositeRheology{T,N}, εII, args; tol=1e-6, verbose=false, τ_initial=nothing, ε_init=nothing, AD=false)
    
"""
function compute_τII(v::CompositeRheology, εII, args; tol=1e-6, verbose=false, τ_initial=nothing, ε_init=nothing, AD=false)
    # A composite rheology case with parallel elements
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose, τ_initial=τ_initial, ε_init=ε_init, AD=AD)
    return τII
end

function compute_τII_AD(v::CompositeRheology, εII, args; tol=1e-6, verbose=false)
     τII = local_iterations_εII_AD(v, εII, args; tol=tol, verbose=verbose)
     return τII
end

compute_τII_AD(v::AbstractMaterialParamsStruct, εII, args) = compute_τII_AD(v.CompositeRheology[1], εII, args)
compute_τII(v::AbstractMaterialParamsStruct, εII, args) = compute_τII(v.CompositeRheology[1], εII, args)

# STRESS AND STRAIN RATE DERIVATIVES
dεII_dτII(v::CompositeRheology, τII, args; tol=1e-6, verbose=false) = nreduce(vi -> first(dεII_dτII(vi, τII, args)), v.elements)

"""
    dεII_dτII_AD(v::Union{Parallel,CompositeRheology}, τII, args) 

Uses AD to compute the derivative of `εII` vs. `τII`
"""
dεII_dτII_AD(v::Union{Parallel,CompositeRheology}, τII, args) = ForwardDiff.derivative(x->compute_εII(v, x, args), τII)
dεII_dτII_nonplastic_AD(v::Union{Parallel,CompositeRheology}, τII, args) = ForwardDiff.derivative(x->compute_εII_nonplastic(v, x, args), τII)

# Computes sum of dεII/dτII for all elements that are NOT parallel and NOT plastic elements
"""
    dεII_dτII_elements(v::CompositeRheology, TauII, args)

Sums the derivative ∂εII/∂τII (strainrate vs. stress) of all non-parallel elements in a `CompositeRheology` structure. Internally used for jacobian iterations.
"""
dεII_dτII_elements(v::CompositeRheology, τII, args; tol=1e-6, verbose=false) = nreduce(vi -> first(dεII_dτII_nonparallel(vi, τII, args)), v.elements)
dεII_dτII_nonparallel(v::Any, TauII, args) =   dεII_dτII(v, TauII, args)
dεII_dτII_nonparallel(v::Parallel, TauII, args) = zero(TauII)
dεII_dτII_nonparallel(v::AbstractPlasticity, TauII, args) = zero(TauII)

dεvol_dp(v::CompositeRheology, p, args) = nreduce(vi -> first(_dεvol_dp(vi, p, args)), v.elements)
function _dεvol_dp(v, p, args)
    if isvolumetric(v)
        val = dεvol_dp_nonparallel_nonplastic(v, p, args)
    else
        val = 0.0
    end
    return val
end
dεvol_dp_nonparallel_nonplastic(v, P, args) = dεvol_dp(v, P, args)
dεvol_dp_nonparallel_nonplastic(v::Parallel, P, args) = zero(P)
dεvol_dp_nonparallel_nonplastic(v::AbstractPlasticity, P, args) = zero(P)


"""
    dτII_dεII(v::CompositeRheology, TauII::_T, args)

Computes the derivative of `τII` vs `εII` for `CompositeRheology`   
"""
function dτII_dεII(v::CompositeRheology, εII, args)
    τ,  = compute_τII(v, εII, args)
    return inv(dεII_dτII(v, τ, args))
end


@generated  function dτII_dεII_i(
    v::CompositeRheology{T,N}, 
    εII, 
    args, I::Int64;
    tol=1e-6, verbose=false
) where {T,N}
    quote
        Base.@_inline_meta
        @assert I ≤ $N
        Base.Cartesian.@nexprs $N i -> I == i && return dτII_dεII(v.elements[i], εII, args)
    end
end

dτII_dεII_AD(v::Union{Parallel,CompositeRheology}, εII, args) = ForwardDiff.derivative(x->compute_τII_AD(v, x, args), εII)


# AVERAGES (mostly used as initial guesses)

"""
    compute_εII_nonplastic(v::CompositeRheology, TauII, args)

Harmonic average of stress of all elements in a `CompositeRheology` structure that are not plastic elements
"""
compute_εII_nonplastic(v::CompositeRheology, τII, args; tol=1e-6, verbose=false) = nreduce(vi -> first(_compute_εII_nonplastic(vi, τII, args)), v.elements)
_compute_εII_nonplastic(v, TauII, args) = first(compute_εII(v, TauII, args))
_compute_εII_nonplastic(v::AbstractPlasticity, TauII, args) = 0.0

"""
    compute_τII_harmonic(v::CompositeRheology, EpsII, args)

Harmonic average of stress of all elements in a `CompositeRheology` structure that are not || elements
"""
function compute_τII_harmonic(v::CompositeRheology, EpsII, args; tol=1e-6, verbose=false)
    return inv(nreduce(vi -> first(_compute_τII_harmonic_element(vi, EpsII, args)), v.elements))
end
_compute_τII_harmonic_element(v, EpsII, args) = inv(first(compute_τII(v, EpsII, args)))
_compute_τII_harmonic_element(v::AbstractPlasticity, EpsII, args) = 0.0
_compute_τII_harmonic_element(v::Parallel{T, N,  Nplast, is_plastic}, EpsII, args) where {T, N,  Nplast, is_plastic}  = 0.0

"""
    compute_p_harmonic(v::CompositeRheology, EpsVol, args)

Harmonic average of stress of all elements in a `CompositeRheology` structure that are not || elements
"""
function compute_p_harmonic(v::CompositeRheology, EpsVol, args; tol=1e-6, verbose=false)
    out = inv(nreduce(vi -> first(_compute_p_harmonic_element(vi, EpsVol, args)), v.elements))
    return out
end
function _compute_p_harmonic_element(v, EpsVol, args) 
    if isvolumetric(v)
       out = inv(compute_p(v, EpsVol, args))
    else
        out = 0.0
    end
    return out
end
_compute_p_harmonic_element(v::AbstractPlasticity, EpsVol, args) = 0.0
_compute_p_harmonic_element(v::Parallel, EpsVol, args) = 0.0

"""
    compute_εII_harmonic(v::Parallel{T,N}, TauII::_T, args)

Computes the harmonic average of strainrate for a parallel element
"""
function compute_εII_harmonic(v::Union{Parallel, CompositeRheology}, TauII, args; tol=1e-6, verbose=false)
    out = inv(nreduce(vi -> inv(first(compute_εII(vi, TauII, args))), v.elements))
    return out
end

@generated  function compute_εII_harmonic_i(
    v::CompositeRheology{T,N}, 
    TauII, 
    args, I::Int64;
    tol=1e-6, verbose=false
) where {T, N}
    quote
        Base.@_inline_meta
        @assert I ≤ $N
        Base.Cartesian.@nexprs $N i -> I == i && return compute_εII(v.elements[i], TauII, args)
    end
end


"""
    compute_εII_elements(v::CompositeRheology, TauII, args)

Sums the strainrate of all non-parallel and non-plastic elements in a `CompositeRheology` structure. Mostly internally used for jacobian iterations.
"""
compute_εII_elements(v::CompositeRheology, τII, args; tol=1e-6, verbose=false) = nreduce(vi -> first(_compute_εII_nonparallel(vi, τII, args)), v.elements)
_compute_εII_nonparallel(v, TauII, args) = compute_εII(v, TauII, args)
_compute_εII_nonparallel(v::Parallel, TauII, args) = zero(TauII)
_compute_εII_nonparallel(v::AbstractPlasticity, TauII, args) = zero(TauII)


"""
    compute_εvol_elements(v::CompositeRheology, TauII, args)

Sums the strainrate of all non-parallel and non-plastic elements in a `CompositeRheology` structure. Mostly internally used for jacobian iterations.
"""
compute_εvol_elements(v::CompositeRheology, p, args; tol=1e-6, verbose=false) = nreduce(vi -> first(_compute_εvol_elements(vi, p, args)), v.elements)
compute_εvol(v::Any, P; kwargs...)= zero(P) # in case nothing more specific is defined
_compute_εvol_elements(v, P, args) = compute_εvol(v, P, args)
_compute_εvol_elements(v::Parallel, P, args) = zero(P)
_compute_εvol_elements(v::AbstractPlasticity, P, args) = zero(P)
