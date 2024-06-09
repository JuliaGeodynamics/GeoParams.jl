import GeoParams: isGPU
using LinearAlgebra

# NONLINEAR ITERATION SCHEMES
"""
    τII =local_iterations_εII(v::CompositeRheology{T,N,0}, εII::_T, args; tol=1e-6, verbose=false)

Performs local iterations versus stress for a given total strain rate for a given `CompositeRheology` element that does NOT include `Parallel` elements
"""
function local_iterations_εII(
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
    # @print(verbose, "initial stress_II = $τII")

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
        τII = @muladd τII + (εII - compute_εII(v, τII, args)) * inv(dεII_dτII(v, τII, args)) 
        ϵ = abs(τII - τII_prev) * inv(τII)
        τII_prev = τII
        # @print(verbose, " iter $(iter) $ϵ")
        
        T_check = ϵ isa Union{AbstractFloat,Integer}
        !(T_check) && break
    end
   
    # @print(verbose, "final τII = $τII")
    # @print(verbose, "---")

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
    
    # @print(verbose, "initial stress_II = $τII")

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

        τII = @muladd f * inv(dεII_dτII) + τII

        ϵ = abs(τII - τII_prev) * inv(τII)
        τII_prev = τII
        # @print(verbose, " iter $(iter) $ϵ τII=$τII")
    end
    # @print(verbose, "final τII = $τII")
    # @print(verbose, "---")

    return τII
end

"""
    compute_εII(v::AbstractPlasticity, τII::_T, args; tol=1e-6, verbose=true)

Performs local iterations to compute the plastic strainrate. Note that the non-plastic strainrate, ε_np, should be part of `args`
"""
function compute_εII(v::AbstractPlasticity, τII::_T, args; tol=1e-6, verbose=false) where _T

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

        # @print(verbose, "    plastic iter $(iter) ϵ=$ϵ λ=$λ, F=$F")
    end

    ε_pl = λ*∂Q∂τII(v, τII_pl, args)
    
    return ε_pl
end


@inline function local_iterations_τII_AD(
    v::Parallel, τII::T, args; tol=1e-6, verbose=false
) where {T}

    # Initial guess
    εII = compute_εII_harmonic(v, τII, args)

    # @print(verbose, "initial εII = $εII")

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
        εII = @muladd (τII - first(compute_τII(v, εII, args))) * inv(dτII_dεII(v, εII, args)) + εII

        ϵ = abs(εII - εII_prev) * inv(εII)
        εII_prev = εII
        # @print(verbose," iter $(iter) $ϵ")
        
    end
    
    # @print(verbose,"final εII = $εII")
    # @print(verbose,"---")
    
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
    
    # @print(verbose,"initial p = $p")

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
        p = @muladd (εvol - compute_εvol(v, p, args)) * inv(dεvol_dp(v, p, args)) + p

        ϵ = abs(p - p_prev) * inv(abs(p))
        p_prev = p

        # @print(verbose," iter $(iter) $ϵ")
    end
    
    # @print(verbose,"final p = $p")
    # @print(verbose,"---")
    
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

        ϵ = abs(εII - εII_prev) / εII
        εII_prev = εII
        # @print(verbose," iter $(iter) $ϵ")
    end
    # @print(verbose,"---")
    
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

    # @print(verbose,"τII guess = $τ_initial")

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
        
        ϵ = sum(abs.(dx)./(abs.(x)))
        # @print(verbose," iter $(iter) $ϵ")
    end
    # @print(verbose,"---")
    
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
    
    # @print(verbose,"τII guess = $τ_initial")
    
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
        # @print(verbose," iter $(iter) $ϵ F=$(r[2]) τ=$(x[1]) λ=$(x[2])")
    end
    # @print(verbose,"---")
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
        P       = x[j+2]    # pressure

        args    = merge(args, (τII=τ_pl, P=P))
        F       = compute_yieldfunction(element,args);  # yield function applied to plastic element
    
        ε̇_pl    =  λ̇*∂Q∂τII(element, τ_pl, args)  
        ε̇vol_pl =  λ̇*∂Q∂P(element, P, args)  
        r[1]   -=  ε̇_pl                     #  contribution of plastic strainrate to residual
        r[j+2] -=  ε̇vol_pl                  #  contribution of vol. plastic strainrate to residual
        
        if F>=0.0
            J[1,j]   =  ∂Q∂τII(element, τ_pl, args)     
            J[j+2,j] =  ∂Q∂P(element, P, args)     
            
            J[j,j]     = ∂F∂λ(element.elements[1], τ_pl, args)        # derivative of F vs. λ
            J[j,j+1]   = ∂F∂τII(element.elements[1], τ_pl, args)    
            J[j,j+2]   = ∂F∂P(element.elements[1], P, args)           # derivative of F vs. P
            
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
    
        ε̇_pl    =  λ̇*∂Q∂τII(element, τ_pl, args)  
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

        args    = merge(args, (τII=τ_pl, P=P))
        F       = compute_yieldfunction(element,args);  # yield function applied to plastic element
    
        ε̇_pl    =  λ̇*∂Q∂τII(element, τ_pl, args)  
        ε̇vol_pl =  λ̇*∂Q∂P(element, P, args)  
        
        r[1]   -=  ε̇_pl                     #  contribution of plastic strainrate to residual
        r[j+1] -=  ε̇vol_pl                  #  contribution of vol. plastic strainrate to residual
        if F>=0.0
            J[1,j] = ∂Q∂τII(element, τ_pl, args)     
            J[3,j] = ∂Q∂P(element, P, args)      
            
            # plasticity is not in a parallel element    
            J[j,1] = ∂F∂τII(element, τ_pl, args)      # derivative of F vs. τ
            J[j,j] = ∂F∂λ(element, τ_pl, args)        # derivative of F vs. λ
            J[j,3] = ∂F∂P(element, P, args)           # derivative of F vs. P
            r[j] =  -F                                # residual
        else
            J[j,j] = 1.0
            r[j] = 0.0
        end
    end

    __fill_J_plastic!(static(is_par), static(is_vol), j, args)

    return j
end

"""
    α,rmnorm = linesearch(f!,x,dx, rnorm; αmin=1e-2, lstol=0.95)

Local linesearch algorithm provided that the in-place residual function is `f!(r,x)`
"""
function linesearch(f!,x,dx, rnorm; αmin=1e-2, lstol=0.95)
    r = zero(x)
    xn = zero(x)
    α, rmnorm = 1.0, 1.0
    while α>αmin
        xn = x + α*dx

        f!(r,xn) # inplace residual computation

        rmnorm = norm(r)
        if rmnorm <= lstol*rnorm
            break
        end

        α  /= 2.0 # reduce bisect step length
    end

    return α, rmnorm
end

# compute the local residual
function Res_εvol_εII!(r,J,x,c,εII_total,εvol_total, args, n)

    τ   = x[1]
    P   = x[n]
    
    # Update part of jacobian related to serial, non-plastic, elements
    r[1]   = εII_total - compute_εII_elements(c,τ,args)     
    J[1,1] = dεII_dτII_elements(c,τ,args);               

    r[n]   = εvol_total - compute_εvol_elements(c,P,args)     
    J[n,n] = dεvol_dp(c,P,args);               

    # Add contributions from plastic elements
    fill_J_plastic!(J, r, x, c, args)

    return nothing
end

"""
    local_iterations_εII(c::CompositeRheology{T,N}, εII_total, εvol_total, args;
                                tol = 1e-6, 
                                verbose = false,
                                τ_initial = nothing, 
                                p_initial = nothing, 
                                ε_init  = nothing,
                                max_iter = 1000)

This performs nonlinear Newton iterations for `τII` with given `εII_total` for cases where we have plastic elements
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

    # Compute residual
    n = 1 + Nplast + Npar + 1;             # total size of unknowns (one added for volumetric plasticity)

    # Initial guess of stress & strainrate
    if isnothing(τ_initial)
        τ_initial = compute_τII_harmonic(c, εII_total, args)
    end
    if isnothing(p_initial)
        p_initial = compute_p_harmonic(c, εvol_total, args)
    end    
    if !isGPU
         @print(verbose,"τII guess = $τ_initial,  P guess = $p_initial")
    end
    x    = @MVector zeros(_T, n)
    x[1] = τ_initial

    set_initial_values!(x, c, τ_initial, p_initial)

    r = @MVector zeros(_T,n);
    J = @MMatrix zeros(_T, n,n)   # size depends on # of plastic elements
    
    
    # Local Iterations
    Res_local_εvol_εII!(r,x)  = Res_εvol_εII!(r,J,x,c,εII_total,εvol_total, args, n)

    iter = 0
    ϵ = 2 * tol
    while (ϵ > tol) && (iter < max_iter)
        iter += 1

        τ   = x[1]
        λ   = x[2]
        P   = x[n]
        
        args = merge(args, (τII=τ, P=P, λ=λ))    # update (to compute yield function)
        
        # Update residual and jacobian
        Res_εvol_εII!(r,J,x,c,εII_total,εvol_total, args, n)
        rnorm_old = norm(r) 
      
        # update solution
        dx  = J\r

        # use linesearch to find optimal stepsize
        # remark: we are currently also updating the jacobian while computing the residual; 
        # that is likely an overkill, so we have speedup potential here
        α, rnorm   = linesearch(Res_local_εvol_εII!, x,dx, rnorm_old; αmin=2e-1, lstol=0.99)
        x .+=  α*dx            

        ϵ    = rnorm
        if !isGPU
            @print(verbose," iter $(iter) $ϵ F=$(-r[2]) τ=$(x[1]) λ=$(x[2]) P=$(x[3]) α=$(α)")
        end
    end
    if !isGPU
        @print(verbose,"---")
    end

    if (iter == max_iter)
        error("iterations did not converge")
    end
    return (x...,)
end



# These lines of code are added as they are slightly faster
@generated function set_initial_values!(x, c::CompositeRheology{T, N, Npar, is_par, Nplast, is_plast}, τ_initial, p_initial) where {T, N, Npar, is_par, Nplast, is_plast}
    quote
        Base.@_inline_meta
        j = 1
        Base.Cartesian.@nexprs $N i -> j = @inbounds _set_initial_values!(x, static($(is_plast)[i]),  $(is_par)[i], τ_initial, p_initial, j)
        return nothing
    end
end

@inline _set_initial_values!(x, ::False, is_par, τ_initial, p_initial, j) = j

@inline function _set_initial_values!(x, ::True, is_par, τ_initial, p_initial, j)

    @inline function __set_initial_values!(x, ::True, τ_initial, p_initial, j)
      # parallel plastic element
      j=j+2
      x[j] = τ_initial    # τ_plastic initial guess     
      j += 1
      x[j] = p_initial    # initial guess for pressure

      return nothing
    end

    @inline function __set_initial_values!(x, ::False, τ_initial, p_initial, j)
        # parallel plastic element
        j += 1
        x[j] = 0    # λ̇  
        j += 1
        x[j] = p_initial    # initial guess for pressure

        return nothing
      end

    __set_initial_values!(x, static(is_par),  τ_initial, p_initial, j)

    return j

end