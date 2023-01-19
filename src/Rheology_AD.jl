# i'm lazy so let's define both functions programatically
for tensor in (:εII, :τII)
    fun1 = Symbol("λ_$(tensor)")
    fun2 = Symbol("compute_$(tensor)")
    @eval begin
        # One can't define lambas inside @generated functions (i.e. as needed for FD.derivative)
        # It'd work with one lambda, but for some obscure reason nested lambdas 
        # yield a slightly more performant algorithm
        function $(fun1)(vi, II, args)
            # lambdas party
            @inline λ(vi, II) = begin
                _λ = II -> $(fun2)(vi, II, args)
                derivative_all(_λ, II)
            end
        
            λ(vi, II)
        end
    end
end

# evalute f(x) and f'(x) of second invariant of tensor A w.r.t. second invariant of tensor B
@generated function compute_II_AD_all(v::NTuple{N, AbstractConstitutiveLaw}, λ::F, BII, args) where {N,F}
    quote
        Base.@_inline_meta
        AII, dεII_dBII = 0.0, 0.0
        Base.@nexprs $N i -> @inbounds ( 
            V = λ(v[i], BII, args);
            AII += V[1];
            dεII_dBII += V[2];
        )
        return AII, dεII_dBII 
    end
end

# thin wrappers to compute the second invariant of the strain rate using forward mode AD
compute_εII_AD_all(v, τII, args) = compute_II_AD_all(v, λ_εII, τII, args)
compute_εII_AD_all(v::Union{Parallel, CompositeRheology}, τII, args) = compute_εII_AD_all(v.elements, τII, args)

# thin wrappers to compute the second invariant of the deviatoric stress using forward mode AD
compute_τII_AD_all(v, εII, args) = compute_II_AD_all(v, λ_τII, εII, args)
compute_τII_AD_all(v::Union{Parallel, CompositeRheology}, εII, args) = compute_τII_AD_all(v.elements, εII, args)

# generic function to compute the second invariant of tensor B 
# given the second invariant of tensor AII
@inline function compute_II_AD(
    v::CompositeRheology{T,N,
        0,is_parallel,
        0,is_plastic,
        Nvol,is_vol,
        false},
    f_harmonic::F1,
    f_tensor_AD::F2,
    AII, 
    args, 
    tol
) where {N, T, is_parallel, is_plastic, Nvol, is_vol, F1, F2}

    # Initial guess
    BII = f_harmonic(v, AII, args)

    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    BII_prev = BII
    while ϵ > tol
        iter += 1
        # Newton scheme -> τII = τII - f(τII)/dfdτII. Therefore,
        #   f(τII) = εII - strain_rate_circuit(v, τII, args) = 0
        #   dfdτII = - dεII_dτII(v, τII, args) 
        #   τII -= f / dfdτII 
        AII_n, dAII_dBII_n = f_tensor_AD(v, BII, args)
        BII = muladd(AII - AII_n, inv(dAII_dBII_n), BII)

        ϵ = abs(BII - BII_prev) * inv(BII)
        BII_prev = BII
    end

    return BII
end

compute_τII_AD(v, εII, args; tol=1e-6) = compute_II_AD(v, compute_τII_harmonic, compute_εII_AD_all, εII, args, tol)
compute_εII_AD(v, τII, args; tol=1e-6) = compute_II_AD(v, compute_εII_harmonic, compute_τII_AD_all, τII, args, tol)

# generic function to compute the second invariant of tensor B 
# given the second invariant of tensor AII
@inline function local_iterations_II_AD(v::Parallel, f_harmonic::F1, f_tensor_AD::F2, AII, args, tol) where {F1, F2}
    # Initial guess
    BII = f_harmonic(v, AII, args)

    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    BII_prev = BII
    while ϵ > tol
        iter += 1
        # Newton scheme -> τII = τII - f(τII)/dfdτII. Therefore,
        #   f(τII) = εII - strain_rate_circuit(v, τII, args) = 0
        #   dfdτII = - dεII_dτII(v, τII, args) 
        #   τII -= f / dfdτII 
        AII_n, dAII_dBII_n = f_tensor_AD(v, BII, args)
        BII = muladd(AII - AII_n, inv(dAII_dBII_n), BII)

        ϵ = abs(BII - BII_prev) * inv(BII)
        BII_prev = BII
    end
    
    return BII
end

local_iterations_τII_AD(v, εII, args; tol=1e-6) = local_iterations_II_AD(v, compute_τII_harmonic, compute_εII_AD_all, εII, args, tol)
local_iterations_εII_AD(v, τII, args; tol=1e-6) = local_iterations_II_AD(v, compute_εII_harmonic, compute_τII_AD_all, τII, args, tol)


# PIRACY TIME FROM CompositeRheology.jl ------------------------------------------------------------

# COMPUTE DEVIATORIC STRESS
"""
    τII = compute_τII_AD(v::CompositeRheology{T,N}, εII, args; tol=1e-6, verbose=false) 
"""
function compute_τII_AD(v::CompositeRheology, εII, args; tol=1e-6)
    # A composite rheology case with parallel elements
    τII = local_iterations_εII_AD(v, εII, args; tol=tol)
    return τII
end

compute_τII(v::CompositeRheology{T,N,0}, εII, args; tol=1e-6) where {T,N} = compute_τII_AD(v, εII, args; tol=tol)

"""
    p,τII = compute_p_τII(v::CompositeRheology, εII, εvol, args; tol=1e-6, verbose=false) 

This updates pressure `p` and deviatoric stress invariant `τII` in case the composite rheology structure has volumetric components, but does not contain plastic or parallel elements.
The 'old' pressure should be stored in `args` as `args.P_old`   
"""
# TODO
function compute_p_τII_AD(
        v::CompositeRheology{T,N,
                    Npar,is_parallel,
                    Nplastic,is_plastic,
                    Nvol,is_vol,
                    false}, 
        εII,
        εvol,
        args; 
        tol=1e-6, verbose=false
    ) where {T, N, Npar, is_parallel, Nplastic, is_plastic, Nvol, is_vol}

    # A composite rheology case that may have volumetric elements, but the are not 
    # tightly coupled, so we do NOT perform coupled iterations.
    τII = local_iterations_εII_AD(v, εII, args; tol=tol)
    P   = local_iterations_εvol_AD(v, εvol, args; tol=tol) # TODO

    return P,τII
end

"""
    p,τII = compute_p_τII(v::CompositeRheology, εII, εvol, args; tol=1e-6, verbose=false) 

This updates pressure `p` and deviatoric stress invariant `τII` in case the composite rheology structure has no volumetric elemnts, but may contain plastic or parallel elements. 
In that case, pressure is not updated (`args.P` is used instead).     
"""
function compute_p_τII_AD(
        v::CompositeRheology{T,N,
                    Npar,is_parallel,
                    Nplast,is_plastic,
                    0,is_vol, 
                    false}, 
        εII, 
        εvol,
        args;
        tol=1e-6
    ) where {T, N, _T, Npar, is_parallel, Nplast, is_plastic, is_vol}
    
    # A composite rheology case with no parallel element; iterations for τII
    τII, = local_iterations_εII_AD(v, εII, args; tol=tol)
    
    P = any(keys(args) .=== :P_old) ? args.P_old : 0.0

    return P, τII
end

"""
    p,τII = compute_p_τII(v::CompositeRheology, εII, εvol, args; tol=1e-6, verbose=false) 

This updates pressure `p` and deviatoric stress invariant `τII` in case the composite rheology structure has volumetric components and has volumetric plasticity
The 'old' pressure should be stored in `args` as `args.P_old`   
"""
# TODO
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
    out = local_iterations_εvol_εII_AD(v, εII, εvol, args; tol=tol, verbose=verbose) #TODO

    τII = out[1]
    P = out[end]

    return P, τII, out[2:N]
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
    tol=1e-6
) where {N, T, _T, Npar, is_par, Nplast, is_plastic, Nvol, is_vol}
    
    # Initial guess
    τII = compute_τII_harmonic(v, εII, args)
    
    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    τII_prev = τII
    ε_pl = 0.0

    while (ϵ > tol) && (iter<10)
        iter += 1
        #  Newton scheme -> τII = τII - f(τII)/dfdτII. 
        #  Therefore,
        #      f(τII) = εII - compute_εII(v, τII, args) = 0
        #      dfdτII = - dεII_dτII(v, τII, args) 
        #      τII -= f / dfdτII
        
        ε_np, dεII_dτII = compute_εII_AD_nonplastic_all(v, τII, args)
        f = εII - ε_np      # non-plastic contributions to residual
        
        if Nplast>0      
            # in case of plasticity, iterate for ε_pl
            args = merge(args, (ε_np=ε_np,f=f))            
            ε_pl += compute_εII_plastic(v, τII, args, tol)

            # add contributions to dεII_dτII:
            if ε_pl > 0.0
                # in fact dε_pl/dτII = d(λ*∂Q∂τII(v_pl, τII))/dτII = 0 for DP
                dεII_dτII += 0
            end
        end
        f -= ε_pl

        τII = muladd(f, inv(dεII_dτII), τII)

        ϵ = abs(τII - τII_prev) * inv(τII)
        τII_prev = τII
    end

    return τII
end

function local_iterations_τII_AD(v::Parallel, τII::T, args; tol=1e-6)
    # Initial guess
    εII = compute_εII_harmonic(v, τII, args)

    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    εII_prev = εII
    while ϵ > tol
        iter += 1
        εII_n, dτII_dεII_n = compute_εII_AD_all(v, εII, args)
        εII = muladd(τII - εII_n, inv(dτII_dεII_n), εII)

        ϵ = abs(εII - εII_prev) * inv(εII)
        εII_prev = εII
        
    end
    
    return εII
end

# NOTE; not working for general composites, only for one single parallel element
function compute_τII_par(
    c::CompositeRheology{
        T,    N,
        Npar, is_par,
        0,    is_plastic,
        0,    is_vol
    },
    εII, 
    args; 
    tol = 1e-6
) where {T, N, Npar, is_par, is_plastic, is_vol}
    
    # number of unknowns  
    n = 1 + Npar; 

    # initial guesses for stress and strain
    εII_total = εII
    εp_n = εII_total 
    τ_n = compute_τII_harmonic(c, εII_total, args)

    # x[1:end-1] = strain rate of parallel(s) element, x[end] = total stress invariant
    x = SVector{n,Float64}(i != n ? εp_n : τ_n for i in 1:n)

    # define closures to be differentiated
    f1(x, args) = εII_total - compute_εII_elements(c, x[2], args) - x[1]
    f2(x, args) = x[2] - compute_τII_parallel(c, x[1], args)
    funs = f1, f2

    # funs = ntuple(Val(n)) do i
    #     f = if i ==1
    #         (x, args) -> εII_total - compute_εII_elements(c, x[end], args) - sum(x[i] for i in 1:Npar)
    #     else
    #         (x, args) -> x[end] - compute_τII_parallel(c.elements, x[i-1], args)
    #     end
    # end

    # Local Iterations
    iter = 0
    ϵ = 2 * tol
    max_iter = 20
    while (ϵ > tol) && (iter < max_iter)
        iter += 1

        x_n = x
        args_n = merge(args, (τII = x[end], )) # we define a new arguments tuple to maintain type stability 
        r, J = jacobian(input -> SInput(funs, input, args_n), x)

        # update solution
        dx  = J\r 
        x  -= dx   
        
        ϵ = norm(@. abs(x-x_n))
    end
    
    return x[end]

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
    εvol,
    args; 
    tol=1e-6
) where {N, T, is_parallel, is_plastic, Nvol, is_vol}

    # Initial guess
    p = compute_p_harmonic(v, εvol, args)
    
    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    p_prev = p
    while ϵ > tol
        iter += 1
        εvol_n, dεvol_dp_n = derivative_all(p -> compute_εvol(v, p, args), p)
        p = muladd(εvol - εvol_n, inv(dεvol_dp_n), p)

        ϵ = abs(p - p_prev) * inv(abs(p))
        p_prev = p

    end
    
    return p
end