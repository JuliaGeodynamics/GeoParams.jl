using GeoParams, ForwardDiff, StaticArrays, LinearAlgebra
import GeoParams.MaterialParameters.ConstitutiveRelationships: compute_τII_parallel, compute_τII_harmonic, compute_εII_elements, compute_τII_elements

# Define a range of rheological components
v1 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
v2 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
v3 = LinearViscous()
v4 = LinearViscous(; η=1e22Pa * s)
e1 = ConstantElasticity()           # elasticity
e2 = SetConstantElasticity(; G=5e10, Kb=1e11)
#pl1= DruckerPrager(C=1e6)                # plasticity
pl1 = DruckerPrager(; C=1e6 / cosd(30))        # plasticity which ends up with the same yield stress as pl3
pl2 = DruckerPrager(; C=1e6, ϕ=0, Ψ=10)      # plasticity
pl3 = DruckerPrager(; C=1e6, ϕ=0)            # plasticity

# Parallel elements
p1 = Parallel(v3, v4)                # linear elements
p2 = Parallel(v1, v2)                # includes nonlinear viscous elements
p3 = Parallel(v1, v2, v3)             # includes nonlinear viscous elements
p4 = Parallel(pl1, LinearViscous(; η=1e20Pa * s)) # viscoplastic regularisation

# CompositeRheologies
c1 = CompositeRheology(v1, v2)
c2 = CompositeRheology(v3, v4)       # two linear rheologies
c3 = CompositeRheology(v1, v2, e1)   # with elasticity
c4 = CompositeRheology(v1, v3, p1)   # with linear || element
c5 = CompositeRheology(v1, v4, p2)   # with nonlinear || element
c6 = CompositeRheology(v1, v4, p1, p1) # with 2 || elements
c7 = CompositeRheology(v4, e1)       # viscoelastic with linear viscosity
c8 = CompositeRheology(v4, e1, pl1)   # with plastic element
c9 = CompositeRheology(v4, e1, p4)    # with visco-plastic parallel element
c10 = CompositeRheology(e1, pl3)      # elastoplastic
c11 = CompositeRheology(e1, Parallel(pl3, LinearViscous(; η=1e19Pa * s)))      # elasto-viscoplastic
c14 = CompositeRheology(SetConstantElasticity(G=1e10, Kb=1e10), LinearViscous(η=1e20), DruckerPrager(C=3e5, Ψ=1))   # case A


# Check derivatives 
vec1 = [c1 c2 c3 c4 c5 p1 p2 p3]
args = (T=900.0, d=100e-6, τII_old=1e6, dt=1e8)
εII, τII = 2e-15, 2e6

v = c4   # problem with c3 (elasticity) and c4 (that has || elements) 

compute_εII(c6, τII, args)
compute_τII(c6, εII, args)

ProfileCanvas.@profview for i in 1:1000000
    compute_εII(v, τII, args)
end

@inline function derivative_all(f::F, x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    res = f(ForwardDiff.Dual{T}(x, one(x)))
    f, ∂f∂x = if res != 0.0
        res.value, res.partials.values[1]
    else
        0.0, 0.0
    end
    return f, ∂f∂x
end

# im lazy so we define both functions here
for tensor in (:εII, :τII)
    fun1 = Symbol("λ_$(tensor)")
    fun2 = Symbol("compute_$(tensor)")
    @eval begin
        function $(fun1)(vi, II, args)
            # trick with nested lambdas so that it works with @generated
            # It'd work with one lambda, but for some reason this is
            # slightly more optimal
            @inline λ(vi, II) = begin
                _λ = II -> $(fun2)(vi, II, args)
                derivative_all(_λ, II)
            end
        
            λ(vi, II)
        end
    end
end

function λ_εII_nonplastic_AD(vi, τII, args)
    # trick with nested lambdas so that it works with @generated
    # It'd work with one lambda, but for some reason this is
    # slightly more optimal
    @inline λ(vi, τII) = begin
        _λ = τII -> _compute_εII_nonplastic(vi, τII, args)
        derivative_all(_λ, τII)
    end

    λ(vi, τII)
end

@generated function compute_II_AD_all(v::NTuple{N, AbstractConstitutiveLaw},
    λ::F,
    II,
    args
) where {N,F}
    quote
        Base.@_inline_meta
        εII, dεII_dτII = 0.0, 0.0
        Base.@nexprs $N i -> (
            V = λ(v[i], II, args);
            εII += V[1];
            dεII_dτII += V[2];
        )
        return εII, dεII_dτII 
    end
end

compute_εII_AD_all(v, τII, args) = compute_II_AD_all(v, λ_εII, τII, args)
compute_εII_AD_all(v::Union{Parallel, CompositeRheology}, τII, args) = compute_εII_AD_all(v.elements, τII, args)

compute_εII_AD_nonplastic_all(v, τII, args) = compute_II_AD_all(v, λ_εII_nonplastic_AD, τII, args)
compute_εII_AD_nonplastic_all(v::Union{Parallel, CompositeRheology}, τII, args) = compute_εII_AD_nonplastic_all(v.elements, τII, args)

compute_τII_AD_all(v, εII, args) = compute_II_AD_all(v, λ_τII, εII, args)
compute_τII_AD_all(v::Union{Parallel, CompositeRheology}, εII, args) = compute_τII_AD_all(v.elements, εII, args)

function compute_τII_AD(
    v::CompositeRheology{T,N,
        0,is_parallel,
        0,is_plastic,
        Nvol,is_vol,
        false},
    εII, 
    args; 
    tol=1e-6
) where {N, T, is_parallel, is_plastic, Nvol, is_vol}

    # Initial guess
    τII = compute_τII_harmonic(v, εII, args)

    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    τII_prev = τII
    while ϵ > tol
        iter += 1
        # Newton scheme -> τII = τII - f(τII)/dfdτII. Therefore,
        #   f(τII) = εII - strain_rate_circuit(v, τII, args) = 0
        #   dfdτII = - dεII_dτII(v, τII, args) 
        #   τII -= f / dfdτII 
        εII_n, dεII_dτII_n = compute_εII_AD_all(v, τII, args)
        τII = muladd(εII - εII_n, inv(dεII_dτII_n), τII)

        ϵ = abs(τII - τII_prev) * inv(τII)
        τII_prev = τII
    end

    return τII
end

@inline function local_iterations_τII_AD_all(
    v::Parallel, τII, args; tol=1e-6
)
    # Initial guess
    εII = compute_εII_harmonic(v, τII, args)

    # Local Iterations
    iter = 0
    ϵ = 2.0 * tol
    εII_prev = εII
    while ϵ > tol
        iter += 1
        # Newton scheme -> τII = τII - f(τII)/dfdτII. Therefore,
        #   f(τII) = εII - strain_rate_circuit(v, τII, args) = 0
        #   dfdτII = - dεII_dτII(v, τII, args) 
        #   τII -= f / dfdτII 
        τII_n, dτII_dεII_n = compute_τII_AD_all(v, εII, args)
        εII = muladd(τII - τII_n, inv(dτII_dεII_n), εII)

        ϵ = abs(εII - εII_prev) * inv(εII)
        εII_prev = εII
    end
    
    return εII
end

## PROTOTYPES

# # create a Static Vector of a vector of evaluated functions
# @generated function SInput(funs::NTuple{N1, Function}, inputs::SVector{N2, T}) where {N1, N2, T}
#     quote
#         Base.Cartesian.@nexprs $N1 i -> f_i = funs[i](inputs...)
#         Base.Cartesian.@ncall $N1 SVector{$N1, $T} f
#     end
# end


# # create a Mutable Vector of a vector of evaluated functions
# @generated function MInput(funs::NTuple{N1, Function}, inputs::MVector{N2, T}) where {N1, N2, T}
#     quote
#         Base.Cartesian.@nexprs $N1 i -> f_i = funs[i](inputs...)
#         Base.Cartesian.@ncall  $N1 MVector{$N1, $T} f
#     end
# end

SInput(funs::NTuple{N1, Function}, x::SVector{N2, T}) where {N1, N2, T} = SVector{N1,T}(funs[i](x) for i in 1:N1)
# SInput(funs::NTuple{N1, Function}, inputs::SVector{N2, T}, args) where {N1, N2, T} = SVector{N1,T}(funs[i](inputs, args) for i in 1:N1)

function SInput(funs::NTuple{N1, Function}, x::SVector{N2, T}, args) where {N1, N2, T} 
    SVector{N1,T}(
        funs[i](
            i == 1 ?  SVector{2, T}(sum(x[i] for i in 1:N2-1), x[end]) : SVector{2, T}(x[1], x[end]), 
            args
        ) 
        for i in 1:N1
    )
end

MInput(funs::NTuple{N1, Function}, x::MVector{N2, T}) where {N1, N2, T} = MVector{N1,T}(funs[i](x) for i in 1:N1)
MInput(funs::NTuple{N1, Function}, x::MVector{N2, T}, args) where {N1, N2, T} = MVector{N1,T}(funs[i](x, args) for i in 1:N1)


@inline function jacobian(f, x::StaticArray)
    T = typeof(ForwardDiff.Tag(f, eltype(x)))
    result = ForwardDiff.static_dual_eval(T, f, x)
    J = ForwardDiff.extract_jacobian(T, result, x)
    f = extract_value(result)
    return f, J
end

extract_value(result::SVector{N, ForwardDiff.Dual{Tag, T, N}}) where {N,T,Tag} = SVector{N,T}(result[i].value for i in 1:N)

########################################################################################################################################## 
##########################################################################################################################################
c=c4

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

    # x[1:end-1] = strain rate of parallel element, x[end] = total stress invariant
    x = SVector{n,Float64}(i == 1 ? εp_n : τ_n for i in 1:n)

    # define closures that will be differentiated
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
    return x[2], x[1]
end
compute_τII_par(c4, εII, args)

@btime compute_τII_par($c, $εII, $args)

compute_τII_par(c, εII, args)

compute_εII(c, τII, args)
compute_τII(c, εII, args)

@btime compute_τII($v, $εII, $args)


fp = ntuple()

ntuple(Val(2)) do i
    if i == 1
        (x, args) -> εII_total - compute_εII_elements(c, x[2], args) - sum(x[i] for i in 1:(2-1))
    else        
        (x, args) -> x[i+1] - compute_τII_parallel(c, x[1], args)
    end
end

macro foos(x, args)

    esc(
        quote
            # Tuple(
                (x,args) -> x[1] 
            #     for i in length(x)
            # )
        end
    )
        
end
