# This holds structures and computational routines for compositional rheologies
using StaticArrays
using Setfield

export CompositeRheology, Parallel, create_rheology_string, print_rheology_matrix
export time_τII_0D, compute_εII_harmonic, compute_τII_AD, isplastic, compute_p_τII, local_iterations_εvol, compute_p_harmonic
export computeViscosity_εII, computeViscosity_εII_AD, compute_yieldfunction
export isvolumetricplastic
import Base.getindex

import GeoParams.Units: nondimensionalize, dimensionalize

@inline isCUDA() =  isdefined(Main,:CUDA) 

"""
    Put rheological elements in parallel 
"""
struct Parallel{T, N,  Nplast, is_plastic, Nvol, is_vol} <: AbstractConstitutiveLaw{T}
    elements::T
end

function Parallel(v::T) where T
    v           = tuple(v...)
    n           =   length(v)

    is_plastic = isa.(v,AbstractPlasticity)     # Is one of the elements a plastic element?
    Nplast = count(is_plastic)

    is_vol = isvolumetric.(v);
    Nvol   =   count(is_vol);
   
    return Parallel{typeof(v),n, Nplast, is_plastic, Nvol, is_vol}(v)
end
Parallel(a,b...) = Parallel((a,b...,)) 


@generated function getindex(p::Parallel{T, N}, I::Int64) where {T,N}
    quote
        Base.@_inline_meta
        @assert I ≤ $N
        Base.Cartesian.@nexprs $N i -> I == i && return p.elements[i]
    end
end


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


# define rules to nondimensionalise this 
function nondimensionalize(MatParam::Union{Parallel,CompositeRheology}, g::GeoUnits{TYPE}) where {TYPE}
    field_new = ();
    field = MatParam.elements;
    for i=1:length(MatParam.elements)
        field_nd  = nondimensionalize(field[i], g) 
        field_new =  tuple(field_new..., field_nd)
    end
    MatParam = set(MatParam, Setfield.PropertyLens{:elements}(), field_new)

    return MatParam
end

function dimensionalize(MatParam::Union{Parallel,CompositeRheology}, g::GeoUnits{TYPE}) where {TYPE}
    field_new = ();
    field = MatParam.elements;
    for i=1:length(MatParam.elements)
        field_nd  = dimensionalize(field[i], g) 
        field_new =  tuple(field_new..., field_nd)
    end
    MatParam = set(MatParam, Setfield.PropertyLens{:elements}(), field_new)

    return MatParam
end

# Print info in the REPL
include("CompositeRheologies_print.jl")

function show(io::IO, g::AbstractComposite)
    #println(io,"Composite rheology:   ")

    # Compose a string with rheological elements, so we have an overview in the REPL
    str = print_rheology_matrix(g)
    println.(str)
    

    return nothing
end

function show(io::IO, a::Parallel)
    println(io,"Parallel:   ")  

    # Compose a string with rheological elements, so we have an overview in the REPL
    str = print_rheology_matrix(a)
    println.(str)

    return nothing
end

function isvolumetric(c::NTuple{N,Any}) where {N}
    ntuple(i-> isvolumetric(c[i]), Val(N))
end


# HELPER FUNCTIONS

# determine if 3 element is plastic or not
isplastic(v) = false;
isplastic(v::Parallel{T, N,  0, is_plastic}) where {T,N,is_plastic} = false;
isplastic(v::Parallel{T, N,  Nplast, is_plastic}) where {T,N,Nplast,is_plastic} = true;
isplastic(v::AbstractPlasticity) = true;
isplastic(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic}) where {T, N,  Npar, is_parallel, Nplast, is_plastic} = true;
isplastic(v::CompositeRheology{T, N,  Npar, is_parallel, 0, is_plastic}) where {T, N,  Npar, is_parallel, is_plastic} = false;

isvolumetric(v) = false;
isvolumetric(v::Parallel{T, N,  Nplast, is_plastic, 0, is_vol}) where {T,N, Nplast,is_plastic, is_vol} = false;
isvolumetric(v::Parallel{T, N,  Nplast, is_plastic, Nvol, is_vol}) where {T,N, Nplast,is_plastic, Nvol, is_vol} = true;
isvolumetric(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic, 0, is_vol}) where {T,N,Npar, is_parallel, Nplast,is_plastic, is_vol} = false;
isvolumetric(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic, Nvol, is_vol}) where {T,N,Npar, is_parallel, Nplast,is_plastic, Nvol, is_vol} = true;

isvolumetricplastic(v::CompositeRheology{T, N,  Npar, is_parallel, Nplast, is_plastic, Nvol, is_vol, volumetricplasticity}) where {T,N,Npar, is_parallel, Nplast,is_plastic, Nvol, is_vol, volumetricplasticity} = volumetricplasticity;



# COMPUTE STRAIN RATE
"""
    compute_εII(v::Parallel{T,N}, τII, args; tol=1e-6, verbose=false, n=1)

Computing `εII` as a function of `τII` for a Parallel elements is (usually) a nonlinear problem
"""
function compute_εII(
    v::Parallel{T,N}, 
    τII::_T, 
    args; 
    tol=1e-6, verbose=false, n=1
) where {T,N,_T}
    εII = local_iterations_τII(v, τII, args; tol=tol, verbose=verbose, n=n)
    return εII
end

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


#compute_εII_AD(v, τII, args; tol=1e-6, verbose=false) = compute_εII(v, τII, args, tol=tol, verbose=verbose)

# Here we do need to do iterations
function compute_εII_AD(v::Parallel, τII, args; tol=1e-6, verbose=false)
    return local_iterations_τII_AD(v, τII, args; tol=tol, verbose=verbose)
end

#COMPUTE VOLUMETRIC STRAIN-RATE

"""
    compute_εvol(v::CompositeRheology{T,N}, p, args; tol=1e-6, verbose=false, n=1)

Computing `εvol` as a function of `p` for a composite element is the sum of the individual contributions
"""
@generated  function compute_εvol(
    v::CompositeRheology{T,N}, 
    p::_T, 
    args; 
    tol=1e-6, verbose=false
) where {T,_T,N}
    quote
        Base.@_inline_meta
        εvol = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            if isvolumetric(v.elements[i])
                εvol += compute_εvol(v.elements[i], p, args)
            end    
        return εvol
    end
end


@generated  function compute_εvol(
    v::CompositeRheology{T,N}, 
    p::Quantity, 
    args; 
    tol=1e-6, verbose=false
) where {T,_T,N}
    quote
        Base.@_inline_meta
        εvol = 0/s
        Base.Cartesian.@nexprs $N i ->
            if isvolumetric(v.elements[i])
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
                    0,is_vol, false}, 
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


"""
    τII = compute_τII(v::CompositeRheology{T,N}, εII, args; tol=1e-6, verbose=false)
    
"""
function compute_τII(v::CompositeRheology, εII, args; tol=1e-6, verbose=false, τ_initial=nothing, ε_init=nothing)
    # A composite rheology case with parallel elements
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose, τ_initial=τ_initial, ε_init=ε_init)
    return τII
end

# For a parallel element, τII for a given εII is the sum of each component
@generated  function compute_τII(
    v::Parallel{T,N}, 
    εII::_T, 
    args;
    tol=1e-6, verbose=false
) where {T,_T,N}
    quote
        Base.@_inline_meta
        τII = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            τII += first(compute_τII(v.elements[i], εII, args))
    end
end
compute_τII_AD(v::Parallel{T,N}, εII::_T, args; tol=1e-6, verbose=false) where {T,N,_T} = first(compute_τII(v, εII, args)) 

# make it work for dimensional cases
@generated  function compute_τII(
    v::Parallel{T,N}, 
    εII::Quantity, 
    args;
    tol=1e-6, verbose=false
) where {T,N}
    quote
        Base.@_inline_meta
        τII = 0Pa
        Base.Cartesian.@nexprs $N i ->
            τII += first(compute_τII(v.elements[i], εII, args))
    end
end


function compute_τII_AD(v::CompositeRheology, εII, args; tol=1e-6, verbose=false)
     τII = local_iterations_εII_AD(v, εII, args; tol=tol, verbose=verbose)
     return τII
end

@inline function compute_τII!(
    τII::AbstractArray{T,nDim},
    v::NTuple{N,AbstractConstitutiveLaw},
    εII::AbstractArray{T,nDim},
    args;
    tol=1e-6, verbose=false
    ) where {T,nDim,N}
    for I in eachindex(τII)
        τII[I] = first(compute_τII(v, εII[I], (; zip(keys(args), getindex.(values(args), I))...)))
    end
end

# VISCOSITY COMPUTATIONS
""" 
    η = computeViscosity_εII(v::Union{Parallel{T,N}, CompositeRheology{T,N}, AbstractConstitutiveLaw}, εII::_T, args; tol=1e-6, verbose=false)

This computes the effective viscosity for a given input rheology `v` and strainrate `εII`
"""
function computeViscosity_εII(
    v::Union{Parallel, CompositeRheology}, 
    εII::_T, 
    args;
    tol=1e-6, verbose=false
) where {_T}
    τII, = compute_τII(v, εII, args; tol=tol, verbose=verbose)
    η    = _T(0.5) * τII * inv(εII)
    return η
end

function computeViscosity_εII(v::T, εII::_T, args; tol=1e-6, verbose=false) where {T<:AbstractConstitutiveLaw,_T}
    τII, = compute_τII(v, εII, args)
    η    = 0.5 * τII * inv(εII)
    return η
end

""" 
    η = computeViscosity_εII_AD(v::Union{Parallel{T,N}, CompositeRheology{T,N}, AbstractConstitutiveLaw}, εII::_T, args; tol=1e-6, verbose=false)

This computes the effective viscosity for a given input rheology `v` and strainrate `εII`, while using AD if necessary
"""
function computeViscosity_εII_AD(
    v::Union{Parallel, CompositeRheology, AbstractConstitutiveLaw}, 
    εII::_T, 
    args;
    tol=1e-6, verbose=false
) where {_T}
    τII = compute_τII_AD(v, εII, args; tol=tol, verbose=verbose)
    η   = _T(0.5) * τII * inv(εII)
    return η
end

function computeViscosity_εII_AD(v::T, εII::_T, args; tol=1e-6, verbose=false) where {T<:AbstractConstitutiveLaw,_T}
    return computeViscosity_εII(v, εII, args) 
end

# NONLINEAR ITERATION SCHEMES
"""
    τII =local_iterations_εII(v::CompositeRheology{T,N,0}, εII::_T, args; tol=1e-6, verbose=false)

Performs local iterations versus stress for a given total strain rate for a given `CompositeRheology` element that does NOT include `Parallel` elements
"""

@inline function local_iterations_εII(
    v::CompositeRheology{T,N,
                    0,is_parallel,
                    0,is_plastic,
                    Nvol,is_vol,
                    false},         # no volumetric plasticity
    εII::_T, 
    args; 
    tol=1e-6, verbose=false
) where {N, T, _T, is_parallel, is_plastic, Nvol, is_vol}

    # Initial guess
    τII = compute_τII_harmonic(v, εII, args)
    
    !isCUDA() && verbose && println("initial stress_II = $τII")

    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    τII_prev = τII
    while ϵ > tol
        iter += 1
        #= 
            Newton scheme -> τII = τII - f(τII)/dfdτII. 
            Therefore,
                f(τII) = εII - strain_rate_circuit(v, τII, args) = 0
                dfdτII = - dεII_dτII(v, τII, args) 
                τII -= f / dfdτII
        =#
        τII = muladd(εII - compute_εII(v, τII, args), inv(dεII_dτII(v, τII, args)), τII)

        ϵ = abs(τII - τII_prev) * inv(abs(τII))
        τII_prev = τII

        !isCUDA() && verbose && println(" iter $(iter) $ϵ")
    end
    if !isCUDA() && verbose
        println("final τII = $τII")
        println("---")
    end

    return τII
end

"""
    τII = local_iterations_εII_AD(v::CompositeRheology{T,N}, εII::_T, args; tol=1e-6, verbose=false)

Performs local iterations versus stress for a given strain rate using AD
"""
@inline function local_iterations_εII_AD(
    v::CompositeRheology{T,
            N,
            Npar,is_par,
            Nplast,is_plastic,
            Nvol,is_vol,
            false},
    εII::_T, 
    args; 
    tol=1e-6, verbose=false
) where {N, T, _T, Npar, is_par, Nplast, is_plastic, Nvol, is_vol}
    
    # Initial guess
    τII = compute_τII_harmonic(v, εII, args)
    
    !isCUDA() && verbose && println("initial τII = $τII")

    # extract plastic element if it exists
    v_pl = v[1]
    if Nplast>0
        for i=1:N
            if is_plastic[i]
                v_pl =  v[i]
            end
        end
    end

    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    τII_prev = τII
    ε_pl = 0.0;
    while (ϵ > tol) && (iter<10)
        iter += 1
        #= 
            Newton scheme -> τII = τII - f(τII)/dfdτII. 
            Therefore,
                f(τII) = εII - compute_εII(v, τII, args) = 0
                dfdτII = - dεII_dτII(v, τII, args) 
                τII -= f / dfdτII
        =#
        
        ε_np = compute_εII_nonplastic(v, τII, args)
        dεII_dτII = dεII_dτII_nonplastic_AD(v, τII, args)

        f = εII - ε_np      # non-plastic contributions to residual
        
        if Nplast>0
           
            # in case of plasticity, iterate for ε_pl
            args = merge(args, (ε_np=ε_np,f=f))
            ε_pl += compute_εII(v_pl, τII, args, tol=tol, verbose=verbose)

            # add contributions to dεII_dτII:
            if ε_pl>0.0
                # in fact dε_pl/dτII = d(λ*∂Q∂τII(v_pl, τII))/dτII = 0 for DP
                dεII_dτII += 0
            end
        end
        f -= ε_pl


        τII = muladd(f, inv(dεII_dτII), τII)

        ϵ = abs(τII - τII_prev) * inv(abs(τII))
        τII_prev = τII
        !isCUDA() && verbose && println(" iter $(iter) $ϵ τII=$τII")
    end
    if !isCUDA() && verbose
        println("final τII = $τII")
        println("---")
    end

    return τII
end

"""
    compute_εII(v::AbstractPlasticity, τII::_T, args; tol=1e-6, verbose=true)

Performs local iterations to compute the plastic strainrate. Note that the non-plastic strainrate, ε_np, should be part of `args`
"""
function compute_εII(v::AbstractPlasticity, τII::_T, args; tol=1e-6, verbose=true) where _T

    
    η_np  = (τII - args.τII_old)/(2.0*args.ε_np)
           
    F    = compute_yieldfunction(v, merge(args, (τII=τII,)))

    iter = 0
    λ = 0.0 
    ϵ = 2.0 * tol
    τII_pl = τII
    while (ϵ > tol) && (iter<100) && (F>0.0)
        #   τII_pl = τII -  2*η_np*λ*∂Q∂τII
        #   F(τII_pl)
        #   dF/dλ = (dF/dτII)*(dτII/dλ) = (dF/dτII)*(2*η_np*∂Q∂τII)
        
        iter += 1
        τII_pl = τII -  2*η_np*λ*∂Q∂τII(v, τII_pl, args)       # update stress
        F      = compute_yieldfunction(v, merge(args, (τII=τII_pl,)))
        
        dFdλ = ∂F∂τII(v, τII, args)*(2*η_np*∂Q∂τII(v, τII, args))
      
        λ -= -F / dFdλ

        ϵ = F

        !isCUDA() && verbose && println("    plastic iter $(iter) ϵ=$ϵ λ=$λ, F=$F")
    end

    ε_pl = λ*∂Q∂τII(v, τII_pl, args)
    
    return ε_pl
end


@inline function local_iterations_τII_AD(
    v::Parallel, τII::T, args; tol=1e-6, verbose=false
) where {T}

    # Initial guess
    εII = compute_εII_harmonic(v, τII, args)

    !isCUDA() && verbose && println("initial εII = $εII")

    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    εII_prev = εII
    while ϵ > tol
        iter += 1
        #= 
            Newton scheme -> τII = τII - f(τII)/dfdτII. 
            Therefore,
                f(τII) = εII - strain_rate_circuit(v, τII, args) = 0
                dfdτII = - dεII_dτII(v, τII, args) 
                τII -= f / dfdτII
        =#
        εII = muladd(τII - first(compute_τII(v, εII, args)), inv(dτII_dεII(v, εII, args)), εII)

        ϵ = abs(εII - εII_prev) * inv(εII)
        εII_prev = εII
        !isCUDA() && verbose && println(" iter $(iter) $ϵ")
        
    end
    if !isCUDA() && verbose
        println("final εII = $εII")
        println("---")
    end

    return εII
end

"""
    p =local_iterations_εvol(v::CompositeRheology{T,N,0}, εvol::_T, args; tol=1e-6, verbose=false)

Performs local iterations versus pressure for a given total volumetric strain rate for a given `CompositeRheology` element that does NOT include `Parallel` elements
"""
@inline function local_iterations_εvol(
    v::CompositeRheology{T,N,
                    0,is_parallel,
                    0,is_plastic,
                    Nvol,is_vol}, 
    εvol::_T, 
    args; 
    tol=1e-6, verbose=false
) where {N, T, _T, is_parallel, is_plastic, Nvol, is_vol}

    # Initial guess
    p = compute_p_harmonic(v, εvol, args)
    
    verbose && println("initial p = $p")

    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    p_prev = p
    while ϵ > tol
        iter += 1
        #= 
            Newton scheme -> τII = τII - f(τII)/dfdτII. 
            Therefore,
                f(τII) = εII - strain_rate_circuit(v, τII, args) = 0
                dfdτII = - dεII_dτII(v, τII, args) 
                τII -= f / dfdτII
        =#
        p = muladd(εvol - compute_εvol(v, p, args), inv(dεvol_dp(v, p, args)), p)

        ϵ = abs(p - p_prev) * inv(abs(p))
        p_prev = p

        verbose && println(" iter $(iter) $ϵ")
    end
    if verbose
        println("final p = $p")
        println("---")
    end

    return p
end


"""
Performs local iterations versus strain rate for a given stress
"""
@inline function local_iterations_τII(
    v::Parallel{T,N}, 
    τII::_T, 
    args; 
    tol=1e-6, 
    verbose=false, n=1
) where {T,N, _T}

    # Initial guess (harmonic average of ε of each element)
    εII = compute_εII_harmonic(v, τII, args) # no allocations 

    # Local iterations
    iter = 0
    ϵ = 2 * tol
    εII_prev = εII

    while ϵ > tol
        iter += 1
        f = τII - first(compute_τII(v, εII, args))
        dfdεII = -dτII_dεII(v, εII, args)
        εII -= f / dfdεII

        ϵ = abs(εII - εII_prev) / abs(εII)
        εII_prev = εII
        if !isCUDA() && verbose
            println(" iter $(iter) $ϵ")
        end
    end
    if !isCUDA() && verbose
        println("---")
    end

    return εII
end



"""
    local_iterations_εII(c::CompositeRheology{T,N}, εII_total, args)

This performs nonlinear Newton iterations for `τII` with given `εII_total` for cases where we have both serial and parallel elements.
"""
@inline function local_iterations_εII(
    c::CompositeRheology{T,
    N,
    Npar,is_par,
    0,is_plastic,
    0,is_vol}, 
    εII_total::_T, 
    args; 
    tol=1e-6, 
    verbose=false,
    τ_initial=nothing, ε_init=nothing
) where {T,N,Npar,is_par, _T, is_plastic, is_vol}
    
    # Compute residual
    n = Npar+1;             # total size of unknowns
    x = zero(εII_total)
    
    # Initial guess of stress & strainrate
    if isnothing(τ_initial)
        τ_initial = compute_τII_harmonic(c, εII_total, args)
    end

    !isCUDA() && verbose && println("τII guess = $τ_initial")

    x    = @MVector ones(_T, n)
    x   .= εII_total
    x[1] = τ_initial

    j = 1;
    for i=1:N
        if is_par[i]
            j += 1
            x[j] = compute_εII_harmonic_i(c, τ_initial, args,i)   
        end
    end
    
    r = @MVector zeros(_T,n);
    J = @MMatrix zeros(_T, Npar+1,Npar+1)   # size depends on # of parallel objects (+ likely plastic elements)
    
    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    τII_prev = τ_initial
    τ_parallel = _T(0)
    max_iter = 1000
    while (ϵ > tol) && (iter < max_iter)
        iter += 1

        τ   = x[1]
  
        # Update part of jacobian related to serial elements
        r[1]   = εII_total - compute_εII_elements(c,τ,args)
        J[1,1] = dεII_dτII_elements(c,x[1],args);
        
        # Add contributions from || elements
        fill_J_parallel!(J, r, x, c, τ, args)
      
        # update solution
        dx  = J\r 
        x .+= dx   
        
        ϵ    = sum(abs.(dx)./(abs.(x)))
        !isCUDA() && verbose && println(" iter $(iter) $ϵ")
    end
    !isCUDA() && verbose && println("---")
    
    if (iter == max_iter)
        error("iterations did not converge")
    end

    return (x...,)
end


"""
    local_iterations_εII(c::CompositeRheology{T,N}, εII_total, args)

This performs nonlinear Newton iterations for `τII` with given `εII_total` for cases where we plastic elements
"""
@inline function local_iterations_εII(
    c::CompositeRheology{T,N,
                        Npar,is_par,            # no ||
                        Nplast, is_plastic,     # with plasticity
                        0,is_vol},              # no volumetric
    εII_total::_T, 
    args; 
    tol = 1e-6, 
    verbose = false,
    τ_initial = nothing, 
    ε_init = nothing,
    max_iter = 1000
) where {T,N,Npar,is_par, _T, Nplast, is_plastic, is_vol}
    
    # Compute residual
    n = 1 + Nplast + Npar;             # total size of unknowns
    x = zero(εII_total)

    # Initial guess of stress & strainrate
    if isnothing(τ_initial)
        τ_initial = compute_τII_harmonic(c, εII_total, args)
    end
    
    verbose && println("τII guess = $τ_initial")
    
    x    = @MVector zeros(_T, n)
    x[1] = τ_initial

    j = 1;
    for i=1:N
        if is_plastic[i] && is_par[i]
            # parallel plastic element
            j=j+2
            x[j] = τ_initial    # τ_plastic initial guess     
        end
    end

    r = @MVector zeros(_T,n);
    J = @MMatrix zeros(_T, n,n)   # size depends on # of plastic elements
    
    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    while (ϵ > tol) && (iter < max_iter)
        iter += 1

        τ   = x[1]
        λ   = x[2]

        args = merge(args, (τII=τ, λ=λ))    # update

        # Update part of jacobian related to serial, non-plastic, elements
        r[1]   = εII_total - compute_εII_elements(c,τ,args)     
        J[1,1] = dεII_dτII_elements(c,x[1],args);               
        
        # Add contributions from plastic elements
        fill_J_plastic!(J, r, x, c, args)
        
        # update solution
        dx  = J\r 
        x .+= dx   
       # @show dx x r J
        
        ϵ    = sum(abs.(dx)./(abs.(x .+ 1e-9)))
        verbose && println(" iter $(iter) $ϵ F=$(r[2]) τ=$(x[1]) λ=$(x[2])")
    end
    verbose && println("---")
    if (iter == max_iter)
        error("iterations did not converge")
    end
    

    return (x...,)
end


# Helper functions
@generated function fill_J_parallel!(J, r, x, c::CompositeRheology{T, N, Npar, is_par, Nplast, is_plast}, τ, args) where {T, N, Npar, is_par, Nplast, is_plast}
    quote
        Base.@_inline_meta
        j = 1
        Base.Cartesian.@nexprs $N i -> j = @inbounds _fill_J_parallel!(J, r, x, c.elements[i], τ, args, $(is_par)[i], j)
        return nothing
    end
end

@inline function _fill_J_parallel!(J, r, x, elements, τ, args, is_par, j)
    !is_par && return j

    j += 1
    εII_p = x[j]
    r[1] -= εII_p
    τ_parallel, = compute_τII(elements, εII_p, args)    
    r[j]        =  (τ - τ_parallel) # residual (stress should be equal)
    J[j,j]      = -dτII_dεII(elements, εII_p, args)
    J[j,1]      =  1.0
    J[1,j]      =  1.0
    
    return j
end

@generated function fill_J_plastic!(J, r, x, c::CompositeRheology{T, N, Npar, is_par, Nplast, is_plastic, Nvol, is_vol}, args) where {T, N, Npar, is_par, Nplast, is_plastic, Nvol, is_vol}
    quote
        Base.@_inline_meta
        j = 1
        Base.Cartesian.@nexprs $N i -> j = @inbounds _fill_J_plastic!(J, r, x, c.elements[i], args, static($(is_plastic)[i]), $(is_par)[i], $(is_vol)[i], j)
        return nothing
    end
end

@inline _fill_J_plastic!(J, r, x, element, args, ::False, is_par, is_vol, j) = j

@inline function _fill_J_plastic!(J, r, x, element, args, ::True, is_par, is_vol, j)

    j += 1
    λ̇  = x[j]

    @inline function __fill_J_plastic!(::True, ::False, j, args) # parallel, non-dilatant element 
        τ       = x[1]
        τ_pl    = x[j+1]    # if the plastic element is in || with other elements, need to explicitly solve for this

        args    = merge(args, (τII=τ_pl,))
        F       = compute_yieldfunction(element, args);  # yield function applied to plastic element
    
        ε̇_pl    =  λ̇*∂Q∂τII(element, τ_pl, args)  
        r[1]   -=  ε̇_pl                     #  add plastic strainrate

        if F>=0.0
            J[1,j] = ∂Q∂τII(element, τ_pl, args)     

            J[j,j]     = ∂F∂λ(element.elements[1], τ_pl, args)        # derivative of F vs. λ
            J[j,j+1]   = ∂F∂τII(element.elements[1], τ_pl, args)    
        
            J[j+1,1]   = -1.0;
            J[j+1,2]   = dτII_dεII_nonplastic(element, τ_pl, args)*∂Q∂τII(element, τ_pl, args) ;
            J[j+1,j+1] = 1.0;
            r[j] = -F
            r[j+1] = τ - compute_τII_nonplastic(element, ε̇_pl, args) - τ_pl                
        else
            J[j,j] =  1.0
            
            # In this case set τ_pl=τ
            J[j+1,j+1] = 1.0
            J[j+1,1] = -1.0
            
            r[j] = r[j+1] = 0.0
        end
    end

    @inline function __fill_J_plastic!(::True, ::True, j, args) # parallel, dilatant element 
        τ       = x[1]
        τ_pl    = x[j+1]    # if the plastic element is in || with other elements, need to explicitly solve for this
        P       = P[j+2]    # pressure

        args    = merge(args, (τII=τ_pl,P=P))
        F       = compute_yieldfunction(element,args);  # yield function applied to plastic element
    
        ε̇_pl    =  λ̇*∂Q∂τII(element, τ_pl, args)  
        ε̇vol_pl =  λ̇*∂Q∂P(element, P, args)  
        r[1]   -=  ε̇_pl                     #  contribution of plastic strainrate to residual
        r[j+2] +=  ε̇vol_pl                  #  contribution of vol. plastic strainrate to residual
        
        if F>=0.0
            J[1,j]   =  ∂Q∂τII(element, τ_pl, args)     
            J[j+2,j] = -∂Q∂P(element, P, args)     
            
            J[j,j]     = ∂F∂λ(element.elements[1], τ_pl, args)        # derivative of F vs. λ
            J[j,j+1]   = ∂F∂τII(element.elements[1], τ_pl, args)    
            J[j,j+2]   = ∂F∂P(element, P, args)           # derivative of F vs. P
            
            J[j+1,1]   = -1.0;
            J[j+1,2]   = dτII_dεII_nonplastic(element, τ_pl, args)*∂Q∂τII(element, τ_pl, args) ;
            J[j+1,j+1] = 1.0;

            r[j] = -F
            r[j+1] = τ - compute_τII_nonplastic(element, ε̇_pl, args) - τ_pl                
        else
            J[j,j] = 1.0
            
            # In this case set τ_pl=τ
            J[j+1,j+1] = 1.0
            J[j+1,1] = -1.0
            
            r[j] = r[j+1] = 0.0
        end
    end

    @inline function __fill_J_plastic!(::False, ::False, j, args) #non-parallel, non-dilatant
        τ_pl    = x[1]    # if the plastic element is in || with other elements, need to explicitly solve for this

        args    = merge(args, (τII=τ_pl,))
        F       = compute_yieldfunction(element,args);  # yield function applied to plastic element
    
        ε̇_pl    =  λ̇*∂Q∂τII(element, τ_pl)  
        r[1]   -=  ε̇_pl                     #  add plastic strainrate
        
        if F>=0.0
            J[1,j] = ∂Q∂τII(element, τ_pl, args)     

            # plasticity is not in a parallel element    
            J[j,1] = ∂F∂τII(element, τ_pl, args)  
            J[j,j] = ∂F∂λ(element, τ_pl, args)        # derivative of F vs. λ
            r[j] =  -F      
        else
            J[j,j] = 1.0
            r[j] = 0.0
        end
    end

    @inline function __fill_J_plastic!(::False, ::True, j, args) #non-parallel, dilatant
        τ_pl    = x[1]    # if the plastic element is in || with other elements, need to explicitly solve for this
        P       = x[3]

        args    = merge(args, (τII=τ_pl,P=P))
        F       = compute_yieldfunction(element,args);  # yield function applied to plastic element
    
        ε̇_pl    =  λ̇*∂Q∂τII(element, τ_pl, args)  
        ε̇vol_pl =  λ̇*∂Q∂P(element, P, args)  
        
        r[1]   -=  ε̇_pl                     #  contribution of plastic strainrate to residual
        r[j+1] +=  ε̇vol_pl                  #  contribution of vol. plastic strainrate to residual
        
        if F>=0.0
            J[1,j] = ∂Q∂τII(element, τ_pl, args)     
            J[3,j] = -∂Q∂P(element, P, args)      # minus sign because of P sign convention
            
            # plasticity is not in a parallel element    
            J[j,1] = ∂F∂τII(element, τ_pl, args)      # derivative of F vs. τ
            J[j,2] = ∂F∂λ(element, τ_pl, args)        # derivative of F vs. λ
            J[j,3] = ∂F∂P(element, P, args)           # derivative of F vs. P
            J[j,j] = 0.0
            
            r[j] =  -F                          # residual
        else
            J[j,j] = 1.0
            r[j] = 0.0
        end
    end

    __fill_J_plastic!(static(is_par), static(is_vol), j, args)

    return j
end



"""
    local_iterations_εII(c::CompositeRheology{T,N}, εII_total, εvol_total, args)

This performs nonlinear Newton iterations for `τII` with given `εII_total` for cases where we plastic elements
"""
@inline function local_iterations_εvol_εII(
    c::CompositeRheology{T,N,
                        Npar,is_par,            # no ||
                        Nplast, is_plastic,     # with plasticity
                        Nvol,is_vol,            # with volumetric    
                        true},                  # with volumetric plasticity        
    εII_total::_T, 
    εvol_total::_T, 
    args; 
    tol = 1e-6, 
    verbose = false,
    τ_initial = nothing, 
    p_initial = nothing, 
    ε_init = nothing,
    max_iter = 1000
) where {T,N,Npar,is_par, _T, Nplast, is_plastic, Nvol, is_vol}
   # println("local iterations for εvol_εII")    

    # Compute residual
    n = 1 + Nplast + Npar + 1;             # total size of unknowns (one added for volumetric plasticity)

    # Initial guess of stress & strainrate
    if isnothing(τ_initial)
        τ_initial = compute_τII_harmonic(c, εII_total, args)
    end
    if isnothing(p_initial)
        p_initial = compute_p_harmonic(c, εvol_total, args)
    end    
    verbose && println("τII guess = $τ_initial \n  P guess = $p_initial")
    
    x    = @MVector zeros(_T, n)
    x[1] = τ_initial

    j = 1;
    for i=1:N

        if is_plastic[i] && is_par[i]
            # parallel plastic element
            j=j+1
            x[j] = τ_initial    # τ_plastic initial guess     
        end
        if is_plastic[i] && !is_par[i]
            # parallel plastic element
            j += 1
            x[j] = 0    # λ̇  
            j += 1
            x[j] = p_initial    # initial guess for pressure
        end
            
    end
#=    
        if is_plastic[i] 
           # plastic element 
           j += 2
           x[j] = 0    # λ̇  

            j += 1
            x[j] = p_initial    # initial guess for pressure
        end
        if is_par[i] 
            j += 1
            x[j] = τ_initial    # τ_plastic initial guess  
        end
       # if is_vol[i]
       #     # dilatant element
       #     j += 1
       #     x[j] = p_initial    # initial guess for pressure
       # end
=#
#end

    r = @MVector zeros(_T,n);
    J = @MMatrix zeros(_T, n,n)   # size depends on # of plastic elements
    
    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    while (ϵ > tol) && (iter < max_iter)
        iter += 1

        τ   = x[1]
        λ   = x[2]
        P   = x[n]
        
        args = merge(args, (τII=τ, P=P, λ=λ))    # update (to compute yield function)

        # Update part of jacobian related to serial, non-plastic, elements
        r[1]   = εII_total - compute_εII_elements(c,τ,args)     
        J[1,1] = dεII_dτII_elements(c,τ,args);               
        
        r[n]   = εvol_total - compute_εvol_elements(c,P,args)     
        J[n,n] = dεvol_dp(c,P,args);               
        
        # Add contributions from plastic elements
        fill_J_plastic!(J, r, x, c, args)

        #@show x r J
        #error("stop here")

        # update solution
        dx  = J\r 
        x .+= dx   
        
        ϵ    = sum(abs.(dx)./(abs.(x .+ 1e-9)))
        verbose && println(" iter $(iter) $ϵ F=$(r[2]) τ=$(x[1]) λ=$(x[2]) P=$(x[3])")
    end
    verbose && println("---")
    if (iter == max_iter)
        error("iterations did not converge")
    end

    return (x...,)
end

@generated function ∂Q∂τII(
    v::Parallel{T, N,  Nplast, is_plastic}, τ::_T, args
) where {_T,T, N,  Nplast, is_plastic}
    quote
        Base.@_inline_meta
        Base.Cartesian.@nexprs $N i -> is_plastic[i] == true && return ∂Q∂τII(v[i],τ, args)
    end
end

@generated function compute_yieldfunction(
    v::Parallel{T, N,  Nplast, is_plastic}, args
) where {T, N,  Nplast, is_plastic}
    quote
        Base.@_inline_meta
        Base.Cartesian.@nexprs $N i -> is_plastic[i] == true && return compute_yieldfunction(v[i],args)
    end
end

compute_yieldfunction(v::Parallel{T, N,  0, is_plastic}, args) where {T, N,  is_plastic} = NaN

 

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


function dεII_dτII(
    v::Parallel{T,N}, τII::_T, args
) where {T,N,_T}
    ε  = compute_εII(v, τII, args)
    return inv(dτII_dεII(v, ε, args))
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


"""
    dτII_dεII(v::Parallel{T,N}, TauII::_T, args)

Computes the derivative of `τII` vs `εII` for parallel elements   
"""
@generated function dτII_dεII(
    v::Parallel{T,N}, 
    TauII::_T, 
    args
) where {T,N, _T}
    quote
        dτII_dεII_der = zero($_T)
        Base.Cartesian.@nexprs $N i ->
            dτII_dεII_der += dτII_dεII(v.elements[i], TauII, args)
        return dτII_dεII_der
    end
end


"""
    dτII_dεII_nonplastic(v::Parallel{T,N}, TauII::_T, args)

Computes the derivative of `τII` vs `εII` for parallel elements that are non-plastic  
"""
@generated function dτII_dεII_nonplastic(
    v::Parallel{T,N}, 
    TauII::_T, 
    args
) where {T,N, _T}
    quote
        dτII_dεII_der = zero($_T)
        Base.Cartesian.@nexprs $N i ->
            dτII_dεII_der += dτII_dεII_nonplastic(v.elements[i], TauII, args)
        return dτII_dεII_der
    end
end

dτII_dεII_nonplastic(v, TauII, args)  = dτII_dεII(v, TauII, args)
dτII_dεII_nonplastic(v::AbstractPlasticity, TauII, args)  = 0.0


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
    compute_τII_nonplastic(v::Parallel, EpsII, args)

Harmonic average of stress of all elements in a `CompositeRheology` structure that are not || elements
"""
@inline @generated function compute_τII_nonplastic(
    v::Parallel{T,N}, 
    EpsII::_T, 
    args
) where {T,N, _T}
    quote
        out = zero(_T)
        Base.Cartesian.@nexprs $N i ->
            out += _compute_τII_nonplastic(v.elements[i], EpsII, args)
        out = out
    end
end

_compute_τII_nonplastic(v, EpsII, args) = first(compute_τII(v, EpsII, args))
_compute_τII_nonplastic(v::AbstractPlasticity, EpsII, args) = 0.0


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