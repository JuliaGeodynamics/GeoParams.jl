# all related to the CompositeRheology struct 



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
isplastic(v) = false;
isplastic(v::AbstractPlasticity) = true;
isplastic(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic}) where {T, N,  Npar, is_parallel, Nplast, is_plastic} = true;
isplastic(v::CompositeRheology{T, N,  Npar, is_parallel, 0, is_plastic}) where {T, N,  Npar, is_parallel, is_plastic} = false;

isvolumetric(v) = false;
isvolumetric(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic, 0, is_vol}) where {T,N,Npar, is_parallel, Nplast,is_plastic, is_vol} = false;
isvolumetric(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic, Nvol, is_vol}) where {T,N,Npar, is_parallel, Nplast,is_plastic, Nvol, is_vol} = true;

isvolumetricplastic(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic, Nvol, is_vol, volumetricplasticity}) where {T,N,Npar, is_parallel, Nplast,is_plastic, Nvol, is_vol, volumetricplasticity} = volumetricplasticity;

"""
    compute_εII(v::CompositeRheology{T,N}, τII, args; tol=1e-6, verbose=false, n=1)

Computing `εII` as a function of `τII` for a composite element is the sum of the individual contributions
"""
@generated  function compute_εII(
    v::CompositeRheology{T,N}, 
    τII::_T, 
    args; 
    tol=1e-6, verbose=false
) where {T,_T,N}
    quote
        Base.@_inline_meta
        εII = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            εII += compute_εII(v.elements[i], τII, args)
    end
end

@generated  function compute_εII(
    v::CompositeRheology{T,N}, 
    τII::Quantity, 
    args; 
    tol=1e-6, verbose=false
) where {T,N}
    quote
        Base.@_inline_meta
        εII = 0/s
        Base.Cartesian.@nexprs $N i ->
            εII += compute_εII(v.elements[i], τII, args)
    end
end

# As we don't do iterations, this is the same
function compute_εII_AD(v::CompositeRheology, τII, args; tol=1e-6, verbose=false)
    return  compute_εII(v, τII, args)
end

#COMPUTE VOLUMETRIC STRAIN-RATE
"""
    compute_εvol(v::CompositeRheology{T,N}, p, args; tol=1e-6, verbose=false, n=1)

Computing `εvol` as a function of `p` for a composite element is the sum of the individual contributions
"""
@generated  function compute_εvol(
    v::CompositeRheology{T,N,
                        Npar,is_par,            # no ||
                        Nplast, is_plastic,     # with plasticity
                        Nvol,is_vol}, 
    p::_T, 
    args; 
    tol=1e-6, verbose=false
) where {T,_T,N, Npar,is_par, Nplast, is_plastic, Nvol, is_vol}
    quote
        Base.@_inline_meta
        εvol = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            if is_vol[i]
                εvol += compute_εvol(v.elements[i], p, args)
            end    
        return εvol
    end
end

@generated  function compute_εvol(
    v::CompositeRheology{_T,N,
                        Npar,is_par,            # no ||
                        Nplast, is_plastic,     # with plasticity
                        Nvol,is_vol},  
    p::Quantity, 
    args; 
    tol=1e-6, verbose=false
) where {_T,N, Npar,is_par, Nplast, is_plastic, Nvol, is_vol}
    quote
        Base.@_inline_meta
        εvol = zero(_T)/s
        Base.Cartesian.@nexprs $N i ->
            if is_vol[i]
                εvol += compute_εvol(v.elements[i], p, args)
            end
        return εvol
    end
end



# COMPUTE DEVIATORIC STRESS AND PRESSURE
function compute_τII(v::CompositeRheology{T,N,0}, εII, args; tol=1e-6, verbose=false) where {T,N}
    # A composite rheology case with no parallel element; iterations for τII
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose)
    return τII
end

"""
    p,τII = compute_p_τII(v::CompositeRheology, εII, εvol, args; tol=1e-6, verbose=false) 

This updates pressure `p` and deviatoric stress invariant `τII` in case the composite rheology structure has volumetric components, but does not contain plastic or parallel elements.
The 'old' pressure should be stored in `args` as `args.P_old`   
"""
function compute_p_τII(
        v::CompositeRheology{T,N,
                    Npar,is_parallel,
                    Nplastic,is_plastic,
                    Nvol,is_vol,
                    false}, 
        εII::_T, 
        εvol::_T,
        args; 
        tol=1e-6, verbose=false
    ) where {T, N, _T, Npar, is_parallel, Nplastic, is_plastic, Nvol, is_vol}

    # A composite rheology case that may have volumetric elements, but the are not 
    # tightly coupled, so we do NOT perform coupled iterations.
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose)
    P   = local_iterations_εvol(v, εvol, args; tol=tol, verbose=verbose)

    return P,τII
end


"""
    p,τII = compute_p_τII(v::CompositeRheology, εII, εvol, args; tol=1e-6, verbose=false) 

This updates pressure `p` and deviatoric stress invariant `τII` in case the composite rheology structure has no volumetric elemnts, but may contain plastic or parallel elements. 
In that case, pressure is not updated (`args.P` is used instead).     
"""
function compute_p_τII(
        v::CompositeRheology{T,N,
                    Npar,is_parallel,
                    Nplast,is_plastic,
                    0,is_vol, 
                    false}, 
        εII::_T, 
        εvol::_T,
        args;
        tol=1e-6, verbose=false
    ) where {T, N, _T, Npar, is_parallel, Nplast, is_plastic, Nvol, is_vol}
    
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
        εII::_T, 
        εvol::_T,
        args; 
        tol=1e-6, verbose=false
    ) where {T, N, _T, Npar, is_parallel, Nplastic, is_plastic, Nvol, is_vol}

    # A composite rheology case that may have volumetric elements, but the are not 
    # tightly coupled, so we do NOT perform coupled iterations.
    out = local_iterations_εvol_εII(v, εII, εvol, args; tol=tol, verbose=verbose)

    τII = out[1]
    P = out[end]

    return P,τII, out[2:N]
end


# COMPUTE STRAIN RATE


"""
    τII = compute_τII(v::CompositeRheology{T,N}, εII, args; tol=1e-6, verbose=false)
    
"""
function compute_τII(v::CompositeRheology, εII, args; tol=1e-6, verbose=false, τ_initial=nothing, ε_init=nothing)
    # A composite rheology case with parallel elements
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose, τ_initial=τ_initial, ε_init=ε_init)
    return τII
end

function compute_τII_AD(v::CompositeRheology, εII, args; tol=1e-6, verbose=false)
     τII = local_iterations_εII_AD(v, εII, args; tol=tol, verbose=verbose)
     return τII
end

# STRESS AND STRAIN RATE DERIVATIVES
@generated function dεII_dτII(
    v::CompositeRheology{T,N}, τII::_T, args
) where {T,_T,N}
    quote
        Base.@_inline_meta
        val = zero(_T)
        Base.Cartesian.@nexprs $N i -> val += dεII_dτII(v.elements[i], τII, args)
        return val
    end
end

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
@inline @generated function dεII_dτII_elements(
    v::CompositeRheology{T,N}, 
    TauII::_T, 
    args
) where {T, N, _T}
    quote
        out = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            out += dεII_dτII_nonparallel(v.elements[i], TauII, args)
    end
end
dεII_dτII_nonparallel(v::Any, TauII, args) =   dεII_dτII(v, TauII, args)
dεII_dτII_nonparallel(v::Parallel, TauII::_T, args) where _T =    zero(_T)
dεII_dτII_nonparallel(v::AbstractPlasticity, TauII::_T, args) where _T =    zero(_T)

@generated function dεvol_dp(
    v::CompositeRheology{T,N}, p::_T, args
) where {T,_T,N}
    quote
        Base.@_inline_meta
        val = zero(_T)
        Base.Cartesian.@nexprs $N i -> 
        if isvolumetric(v.elements[i])
            val += dεvol_dp_nonparallel_nonplastic(v.elements[i], p, args)
        end
        return val
    end
end
dεvol_dp_nonparallel_nonplastic(v::Any, P, args) =   dεvol_dp(v, P, args)
dεvol_dp_nonparallel_nonplastic(v::Parallel, P::_T, args) where _T =    zero(_T)
dεvol_dp_nonparallel_nonplastic(v::AbstractPlasticity, P::_T, args) where _T =    zero(_T)


"""
    dτII_dεII(v::CompositeRheology, TauII::_T, args)

Computes the derivative of `τII` vs `εII` for `CompositeRheology`   
"""
function dτII_dεII(
    v::CompositeRheology{T,N}, εII::_T, args
) where {T,N,_T}
    τ,  = compute_τII(v, εII, args)
    return inv(dεII_dτII(v, τ, args))
end


@generated  function dτII_dεII_i(
    v::CompositeRheology{T,N}, 
    εII::_T, 
    args, I::Int64;
    tol=1e-6, verbose=false
) where {T,_T,N}
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
@inline @generated function compute_εII_nonplastic(
    v::CompositeRheology{T,N}, 
    TauII::_T, 
    args
) where {T,N, _T}
    quote
        out = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            out += _compute_εII_nonplastic(v.elements[i], TauII, args)
        out = out
    end
end

_compute_εII_nonplastic(v, TauII, args) = first(compute_εII(v, TauII, args))
_compute_εII_nonplastic(v::AbstractPlasticity, TauII, args) = 0.0



"""
    compute_τII_harmonic(v::CompositeRheology, EpsII, args)

Harmonic average of stress of all elements in a `CompositeRheology` structure that are not || elements
"""
@inline @generated function compute_τII_harmonic(
    v::CompositeRheology{T,N}, 
    EpsII::_T, 
    args
) where {T,N, _T}
    quote
        out = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            out += _compute_τII_harmonic_element(v.elements[i], EpsII, args)
        out = 1/out
    end
end

_compute_τII_harmonic_element(v, EpsII, args) = inv(first(compute_τII(v, EpsII, args)))
_compute_τII_harmonic_element(v::AbstractPlasticity, EpsII, args) = 0.0
_compute_τII_harmonic_element(v::Parallel{T, N,  Nplast, is_plastic}, EpsII, args) where {T, N,  Nplast, is_plastic}  = 0.0

"""
    compute_p_harmonic(v::CompositeRheology, EpsVol, args)

Harmonic average of stress of all elements in a `CompositeRheology` structure that are not || elements
"""
@inline @generated function compute_p_harmonic(
    v::CompositeRheology{T,N}, 
    EpsVol::_T, 
    args
) where {T,N, _T}
    quote
        out = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            if isvolumetric(v.elements[i])
                out += _compute_p_harmonic_element(v.elements[i], EpsVol, args)
            end
        out = 1/out
    end
end

_compute_p_harmonic_element(v, EpsVol, args) = inv(compute_p(v, EpsVol, args))
_compute_p_harmonic_element(v::AbstractPlasticity, EpsVol, args) = 0.0
_compute_p_harmonic_element(v::Parallel, EpsVol, args) = 0.0

"""
    compute_εII_harmonic(v::Parallel{T,N}, TauII::_T, args)

Computes the harmonic average of strainrate for a parallel element
"""
@generated function compute_εII_harmonic(
    v::Union{Parallel{T,N},CompositeRheology{T,N}}, 
    TauII::_T, 
    args
) where {T,N, _T}
    quote
        out = zero($_T)
        Base.Cartesian.@nexprs $N i ->
            out += inv(first(compute_εII(v.elements[i], TauII, args)))
        return inv(out)
    end
end

@generated  function compute_εII_harmonic_i(
    v::CompositeRheology{T,N}, 
    TauII::_T, 
    args, I::Int64;
    tol=1e-6, verbose=false
) where {T,_T,N}
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
@inline @generated function compute_εII_elements(
    v::CompositeRheology{T,N}, 
    TauII::_T, 
    args;
    verbose=false
) where {T,N, _T}
    quote
        out = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            out += _compute_εII_nonparallel(v.elements[i], TauII, args)
    end
end

_compute_εII_nonparallel(v, TauII::_T, args) where {_T} = compute_εII(v, TauII, args)
_compute_εII_nonparallel(v::Parallel, TauII::_T, args) where {_T} = zero(_T)
_compute_εII_nonparallel(v::AbstractPlasticity, TauII::_T, args) where {_T} = zero(_T)


"""
    compute_εvol_elements(v::CompositeRheology, TauII, args)

Sums the strainrate of all non-parallel and non-plastic elements in a `CompositeRheology` structure. Mostly internally used for jacobian iterations.
"""
@inline @generated function compute_εvol_elements(
    v::CompositeRheology{T,N}, 
    P::_T, 
    args;
    verbose=false
) where {T,N, _T}
    quote
        out = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            out += _compute_εvol_elements(v.elements[i], P, args)
    end
end

_compute_εvol_elements(v, P::_T, args) where {_T} = compute_εvol(v, P, args)
_compute_εvol_elements(v::Parallel, P::_T, args) where {_T} = zero(_T)
_compute_εvol_elements(v::AbstractPlasticity, P::_T, args) where {_T} = zero(_T)

compute_εvol(v::Any, P::_T; kwargs...) where _T = _T(0) # in case nothing more specific is defined
