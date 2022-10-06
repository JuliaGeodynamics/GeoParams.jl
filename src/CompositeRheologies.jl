# This holds structures and computational routines for compositional rheologies
using StaticArrays
using Setfield

export CompositeRheology, Parallel, create_rheology_string, print_rheology_matrix
export time_τII_0D, compute_εII_harmonic, compute_τII_AD, compute_p_τII, local_iterations_εvol, compute_p_harmonic 
export computeViscosity_εII, computeViscosity_εII_AD
import Base.getindex

import GeoParams.Units: nondimensionalize, dimensionalize

"""
    Put rheological elements in parallel 
"""
struct Parallel{T, N} <: AbstractConstitutiveLaw{T}
    elements::T
end

function Parallel(v::T) where T
    v = tuple(v...)
    return Parallel{typeof(v),length(v)}(v)
end
Parallel(a,b...) = Parallel((a,b...,)) 


@generated function getindex(p::Parallel{T, N}, I::Int64) where {T,N}
    quote
        Base.@_inline_meta
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
                        Nvol, is_vol
                        } <: AbstractComposite
    elements::T
end

# Defines tuples of composite rheologies, while also checking which type of iterations need to be performed
function CompositeRheology(v::T) where {T}
    @show typeof(v)

    # determine if we have parallel elements & if yes: where
    id_parallel =   findall(isa.(v, Parallel));
    Npar        =   length(id_parallel)
    n           =   length(v)
    par         =   zeros(Bool,n)
    if Npar>0
        par[id_parallel] .= 1
    end
    is_parallel =   SVector{n,Bool}(par)

    # determine if we have plastic elements 
    # NOTE: we likely have to expand this to include parallel elements that have plasticity
    id_plastic  =   findall(isa.(v, AbstractPlasticity))
    Nplast      =   length(id_plastic)
    plastic     =   zeros(Bool,n)
    if Nplast>0
        plastic[id_plastic] .= 1
    end
    is_plastic  =   SVector{n,Bool}(plastic)
    
    # determine if we have elements that have volumetric deformation
    # TO BE EXPANDED 

    #Base.Cartesian.@nexprs $N i -> isvolumetric(v.elements[i])
    id_vol      =   findall(isvolumetric(v))
    Nvol        =   length(id_vol)
    volum       =   zeros(Bool,n)
    if Nvol>0
        volum[id_vol] .= 1
    end
    is_vol      =   SVector{n,Bool}(volum)
     
    return CompositeRheology{typeof(v), n, Npar, is_parallel, Nplast, is_plastic, Nvol, is_vol}(v)
end
CompositeRheology(a,b...) = CompositeRheology( (a,b...,)) 
CompositeRheology(a::Parallel) = CompositeRheology( (a,)) 

@generated function getindex(p::CompositeRheology{T, N}, I::Int64) where {T,N}
    quote
        Base.@_inline_meta
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

function compute_p_τII(
    v::CompositeRheology{T,N,
                    0,is_parallel,
                    0,is_plastic,
                    Nvol,is_vol}, 
        εII::_T, 
        εvol,
        args; 
        tol=1e-6, verbose=false
    ) where {T, N, _T, is_parallel, is_plastic, Nvol, is_vol}
    # A composite rheology case with no parallel element; iterations for τII
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose)
    p   = local_iterations_εvol(v, εvol, args; tol=tol, verbose=verbose)
    return p,τII
end

function compute_p_τII(
    v::CompositeRheology{T,N,
                    0,is_parallel,
                    0,is_plastic,
                    0,is_vol}, 
        εII::_T, 
        εvol,
        args;
        tol=1e-6, verbose=false
    ) where {T, N, _T, is_parallel, is_plastic, Nvol, is_vol}
    # A composite rheology case with no parallel element; iterations for τII
    τII = local_iterations_εII(v, εII, args; tol=tol, verbose=verbose)
    p = any(keys(args) .=== :p_old) ? args.p_old : 0
    return p,τII
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
            τII += compute_τII(v.elements[i], εII, args)
    end
end
compute_τII_AD(v::Parallel{T,N}, εII::_T, args; tol=1e-6, verbose=false) where {T,N,_T} = compute_τII(v, εII, args) 

# make it work for dimensional cases
@generated  function compute_τII(
    v::Parallel{T,N}, 
    εII::Quantity, 
    args;
    tol=1e-6, verbose=false
) where {T,_T,N}
    quote
        Base.@_inline_meta
        τII = 0Pa
        Base.Cartesian.@nexprs $N i ->
            τII += compute_τII(v.elements[i], εII, args)
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
        τII[I] = compute_τII(v, εII[I], (; zip(keys(args), getindex.(values(args), I))...))
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
    τII = compute_τII(v, εII, args; tol=tol, verbose=verbose)
    η   = _T(0.5) * τII * inv(εII)
    return η
end

function computeViscosity_εII(v::T, εII::_T, args; tol=1e-6, verbose=false) where {T<:AbstractConstitutiveLaw,_T}
    τII = compute_τII(v, εII, args)
    η   = 0.5 * τII * inv(εII)
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
    v::CompositeRheology{T,N,0,is_par,0,is_plastic,0,is_vol}, 
    εII::_T, 
    args; 
    tol=1e-6, verbose=false
) where {T,N,_T,is_par, is_plastic, is_vol}

    # Initial guess
    τII = compute_τII_harmonic(v, εII, args)
    
    verbose && println("initial τII = $τII")

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

        verbose && println(" iter $(iter) $ϵ")
    end
    if verbose
        println("final τII = $τII")
        println("---")
    end

    return τII
end

@inline function local_iterations_εII(
    v::CompositeRheology{T,N,0,is_par,0,is_plastic, Nvol,is_vol}, 
    εII::_T, 
    args; 
    tol=1e-6, verbose=false
) where {T,N,_T,is_par, is_plastic,Nvol, is_vol}

    # Initial guess
    τII = compute_τII_harmonic(v, εII, args)
    
    verbose && println("initial τII = $τII")

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

        verbose && println(" iter $(iter) $ϵ")
    end
    if verbose
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
    v::CompositeRheology{T,N}, 
    εII::_T, 
    args; 
    tol=1e-6, verbose=false
) where {N, T, _T}

    # Initial guess
    τII = compute_τII_harmonic(v, εII, args)
    
    verbose && println("initial τII = $τII")

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
        τII = muladd(εII - compute_εII(v, τII, args), inv(dεII_dτII_AD(v, τII, args)), τII)

        ϵ = abs(τII - τII_prev) * inv(abs(τII))
        τII_prev = τII

        verbose && println(" iter $(iter) $ϵ")
    end
    if verbose
        println("final τII = $τII")
        println("---")
    end

    return τII
end

@inline function local_iterations_τII_AD(
    v::Parallel, τII::T, args; tol=1e-6, verbose=false
) where {T}
    # Initial guess
    εII = compute_εII_harmonic(v, τII, args)

    verbose && println("initial εII = $εII")

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
        εII = muladd(τII - compute_τII(v, εII, args), inv(dτII_dεII(v, εII, args)), εII)

        ϵ = abs(εII - εII_prev) * inv(εII)
        εII_prev = εII
        verbose && println(" iter $(iter) $ϵ")
        
    end
    if verbose
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
        f = τII - compute_τII(v, εII, args)
        dfdεII = -dτII_dεII(v, εII, args)
        εII -= f / dfdεII

        ϵ = abs(εII - εII_prev) / abs(εII)
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
    local_iterations_εII(c::CompositeRheology{T,N}, εII_total, args)

This performs nonlinear Newton iterations for `τII` with given `εII_total` for cases where we have both serial and parallel elements.
"""
@inline function local_iterations_εII(
    c::CompositeRheology{T,N,Npar,is_par,0,is_plastic,0,is_vol}, 
    εII_total::_T, 
    args; 
    tol=1e-6, 
    verbose=false, τ_initial=nothing, ε_init=nothing
) where {T,N,Npar,is_par, _T, is_plastic, is_vol}
    
    # Compute residual
    n = Npar+1;             # total size of unknowns
    x = zero(εII_total)
    
    # Initial guess of stress & strainrate
    if isnothing(τ_initial)
        τ_initial = compute_τII_harmonic(c, εII_total, args)
    end

    verbose && println("τII guess = $τ_initial")

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
        verbose && println(" iter $(iter) $ϵ")
    end
    verbose && println("---")
    
    if (iter == max_iter)
        error("iterations did not converge")
    end

    τII = x[1]

    return τII
end


# NOTE: we should dispatch on plasticity flag in CompositeRheology struct
"""
    local_iterations_εII(c::CompositeRheology{T,N}, εII_total, args)

This performs nonlinear Newton iterations for `τII` with given `εII_total` for cases where we plastic elements
"""
@inline function local_iterations_εII(
    c::CompositeRheology{T,N,
                        0,is_par,               # no ||
                        Nplast, is_plastic,     # with plasticity
                        0,is_vol},              # no volumetric
    εII_total::_T, 
    args; 
    tol = 1e-6, 
    verbose = false, 
    τ_initial = nothing, 
    ε_init = nothing,
    max_iter = 1000
) where {T,N,is_par, _T, Nplast, is_plastic, is_vol}

    # Compute residual
    n = Nplast+1;             # total size of unknowns
    x = zero(εII_total)

    # Initial guess of stress & strainrate
    if isnothing(τ_initial)
        τ_initial = compute_τII_harmonic(c, εII_total, args)
    end
    
    verbose && println("τII guess = $τ_initial")
    
    x    = @MVector ones(_T, n)
    x   .= εII_total
    x[1] = τ_initial

    j = 1;
    for i=1:N
        #if is_par[i]
        #    j += 1
        #    x[j] = compute_εII_harmonic_i(c, τ_initial, args,i)   
        #end
        if is_plastic[i]
            j += 1
            x[j] = 0        # λ̇  
        end

    end
    
    r = @MVector zeros(_T,n);
    J = @MMatrix zeros(_T, n,n)   # size depends on # of plastic elements
    
    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    τII_prev = τ_initial
    τ_parallel = _T(0)
    while (ϵ > tol) && (iter < max_iter)
        iter += 1

        τ   = x[1]

        args = merge(args, (τII=τ,))    # update

        # Update part of jacobian related to serial, non-plastic, elements
        r[1]   = εII_total - compute_εII_elements(c,τ,args)     
        J[1,1] = dεII_dτII_elements(c,x[1],args);               
        
        # Add contributions from plastic elements
        fill_J_plastic!(J, r, x, c, τ, args)
        
        # update solution
        dx  = J\r 
        x .+= dx   
        
        ϵ    = sum(abs.(dx)./(abs.(x .+ 1e-9)))
        verbose && println(" iter $(iter) $ϵ")
    end
    verbose && println("---")
    if (iter == max_iter)
        error("iterations did not converge")
    end

    τII = x[1]

    return τII
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
    τ_parallel = compute_τII(elements, εII_p, args)    
    r[j]       =  (τ - τ_parallel) # residual (stress should be equal)
    J[j,j]     = -dτII_dεII(elements, εII_p, args)
    J[j,1]     =  1.0
    J[1,j]     =  1.0
    
    return j
end

@generated function fill_J_plastic!(J, r, x, c::CompositeRheology{T, N, Npar, is_par, Nplast, is_plastic}, τ, args) where {T, N, Npar, is_par, Nplast, is_plastic}
    quote
        Base.@_inline_meta
        j = 1
        Base.Cartesian.@nexprs $N i -> j = @inbounds _fill_J_plastic!(J, r, x, c.elements[i], τ, args, $(is_plastic)[i], j)
        return nothing
    end
end

@inline function _fill_J_plastic!(J, r, x, element, τ, args, is_plast, j)
    !is_plast && return j

    j       += 1
    λ̇       = x[j]
    F       = compute_yieldfunction(element,args);
    r[1]   -=  λ̇*∂Q∂τII(element, τ)         # add plastis strainrate
    if F>0
        J[1,j] = ∂Q∂τII(element, τ)     
        J[j,1] = ∂F∂τII(element, τ)    
        J[j,j] = 0 #-F
        r[j] =  -F #-  λ̇*F
    else
        J[j,j] = 1.0
        r[2] = 0.0
    end

    return j
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
dεII_dτII_AD(v::Union{Parallel,CompositeRheology}, τII, args) = ForwardDiff.derivative(x->compute_εII_AD(v, x, args), τII)

# Computes sum of dεII/dτII for all elements that are NOT parallel elements
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
            val += dεvol_dp(v.elements[i], p, args)
        end
        return val
    end
end

"""
    dτII_dεII(v::CompositeRheology, TauII::_T, args)

Computes the derivative of `τII` vs `εII` for `CompositeRheology`   
"""
function dτII_dεII(
    v::CompositeRheology{T,N}, εII::_T, args
) where {T,N,_T}
    τ  = compute_τII(v, εII, args)
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

dτII_dεII_AD(v::Union{Parallel,CompositeRheology}, εII, args) = ForwardDiff.derivative(x->compute_τII_AD(v, x, args), εII)



# HARMONIC AVERAGES (mostly used as initial guesses)

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

_compute_τII_harmonic_element(v, EpsII, args) = inv(compute_τII(v, EpsII, args))
_compute_τII_harmonic_element(v::AbstractPlasticity, EpsII, args) = 0.0

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
            out += inv(compute_εII(v.elements[i], TauII, args))
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


# 0D RHEOLOGY functions

"""
    t_vec, τ_vec = time_τII_0D(v::CompositeRheology, εII::Number, args; t=(0.,100.), τ0=0., nt::Int64=100)

This performs a 0D constant strainrate experiment for a composite rheology structure `v`, and a given, constant, strainrate `εII` and rheology arguments `args`.
The initial stress `τ0`, the time range `t` and the number of timesteps `nt` can be modified 
"""
function time_τII_0D(v::Union{CompositeRheology,Tuple, Parallel}, εII::Number, args; t=(0.,100.), τ0=0., nt::Int64=100, verbose=true)
    t_vec    = range(t[1], t[2], length=nt)
    τ_vec    = zero(t_vec)
    εII_vec  = zero(t_vec) .+ εII
    τ_vec[1] = τ0;

    time_τII_0D!(τ_vec, v, εII_vec, args, t_vec, verbose=verbose);

    return t_vec, τ_vec
end

"""
    time_τII_0D!(τ_vec::Vector{T}, v::CompositeRheology, εII_vec::Vector{T}, args, t_vec::AbstractVector{T}) where {T}

Computes stress-time evolution for a 0D (homogeneous) setup with given strainrate vector (which can vary with time).
"""
function time_τII_0D!(τ_vec::Vector{T}, v::Union{CompositeRheology,Tuple, Parallel}, εII_vec::Vector{T}, args, t_vec::AbstractVector{T}; verbose=false) where {T}

    nt  = length(τ_vec)
    τII = τ_vec[1]

    for i=2:nt  
        dt      = t_vec[i]-t_vec[i-1]
        args    = merge(args, (; τII_old=τII, dt=dt))
        τII     = compute_τII(v, εII_vec[i-1], args, verbose=verbose)
        
        τ_vec[i] = τII
    end


    return nothing
end