using GeoParams, ForwardDiff
import GeoParams.MaterialParameters.ConstitutiveRelationships: compute_τII_harmonic

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
c6 = CompositeRheology(v1, v4, p1, p2) # with 2 || elements
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

v = c3   # problem with c3 (elasticity) and c4 (that has || elements) 

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


# Define a range of rheological components
v1 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
v2 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
e1 = ConstantElasticity()           # elasticity
c = CompositeRheology(v1, v2, e1)   # with elasticity

compute_τII_AD(c, εII, args) == compute_τII(c, εII, args)

@btime compute_τII_AD($c, $εII, $args)
@btime compute_τII($c, $εII, $args)

@btime compute_τII_AD($c7, $εII, $args)
@btime compute_τII($c7, $εII, $args)

ProfileCanvas.@profview for i in 1:1000000
    local_iterations_τII_AD_all(p, εII, args)
end
@edit compute_τII(c4, εII, args)

p=c4.elements[3]
@btime local_iterations_τII_AD_all($p, $εII, $args)
@btime        local_iterations_τII($p, $εII, $args)


p=rand()
derivative_all(p -> compute_εvol(c14, p, args), p)

c = CompositeRheology((v1,  e1, Parallel(pl1, v1, v2)))

ProfileCanvas.@profview for i in 1:100000
    compute_τII(c, εII, args)
end

"""
    local_iterations_εII(c::CompositeRheology{T,N}, εII_total, args)

This performs nonlinear Newton iterations for `τII` with given `εII_total` for cases where we have both serial and parallel elements.
"""
@inline function foo_AD(
    c::CompositeRheology{T,N,
                        Npar,is_par,            # no ||
                        Nplast, is_plastic,     # with plasticity
                        0,is_vol},              # no volumetric
    εII_total::_T, 
    args; 
    tol = 1e-6, 
    τ_initial = nothing, 
    max_iter = 1000
) where {T,N,Npar,is_par, _T, Nplast, is_plastic, is_vol}
    
    # Compute residual
    n = 1 + Nplast + Npar;             # total size of unknowns
    x = zero(_T)

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
    
end

@btime foo_AD($c, $εII, $args) # original non-AD version 402.00 ns

function fa(w::AbstractVector{T}) where {T}
    t = norm(w)
    if t == 0
        (@SVector [1.0, 0.0, 0.0]) * one(T)
    else
        (@SVector [w[1], -w[2], w[3]]) / t
    end
end
pp = (@SVector [0.2, -0.1, 0.1])
ForwardDiff.jacobian(fa, pp)

f1(x, y) = 2.1*x + 1.2*y
f2(x, y) = x + y - 3.3

function f(x)
    @SVector [
        f1(x[1], x[2]),
        f2(x[1], x[2]),
    ]
end
@btime ForwardDiff.jacobian($f, $input)

# create a Static Vector of a vector of evaluated functions
@generated function SInput(funs::NTuple{N1, Function}, inputs::SVector{N2, T}) where {N1, N2, T}
    quote
        Base.Cartesian.@nexprs $N1 i -> f_i = funs[i](inputs...)
        Base.Cartesian.@ncall $N1 SVector{$N1, $T} f
    end
end

# create a Mutable Vector of a vector of evaluated functions
@generated function MInput(funs::NTuple{N1, Function}, inputs::MVector{N2, T}) where {N1, N2, T}
    quote
        Base.Cartesian.@nexprs $N1 i -> f_i = funs[i](inputs...)
        Base.Cartesian.@ncall  $N1 MVector{$N1, $T} f
    end
end

input = @SVector [1.0, 2.0]
function _foo(input::T) where T
    fs = (f1, f2)
    # result = DiffResults.JacobianResult(input)
    # result = ForwardDiff.jacobian!(result, x-> SInput(fs, x), input)
    myjacobian(x-> SInput(fs, x), input)

    # DiffResults.value(result), DiffResults.jacobian(result)
end

foo() = _foo( @SVector [1.0, 2.0])
val1, J1 = _foo(input)

f1(x, y) = sin(exp(x)) + cos(y)/log10(x+y)
f2(x, y) = cos(x + y) * tanh(x*y)

df1dx(x, y) = exp(x)*cos(exp(x)) - 0.43429448190325176((x + y)^-1)*(cos(y) / (log10(x + y)^2))
df1dy(x, y) = (-sin(y)) / log10(x + y) - 0.43429448190325176((x + y)^-1)*(cos(y) / (log10(x + y)^2))
df2dx(x, y) = y*(1 - (tanh(x*y)^2))*cos(x + y) - sin(x + y)*tanh(x*y)
df2dy(x, y) = x*(1 - (tanh(x*y)^2))*cos(x + y) - sin(x + y)*tanh(x*y)

function bar(x)
    f = @SVector [
        f1(x[1], x[2]),
        f2(x[1], x[2]),
    ]

    df = @SMatrix [
        df1dx(x[1], x[2]) df1dy(x[1], x[2])
        df2dx(x[1], x[2]) df2dy(x[1], x[2])
    ]
    f, df
end
@btime bar($input)
val2, J2 = bar(input)

@variables x y
J=Symbolics.jacobian([f1(x, y), f2(x, y)],[x, y])


@edit ForwardDiff.jacobian( x-> SInput(fs, x), input)


@inline myjacobian(f, x::StaticArray) = myvector_mode_jacobian(f, x)

@inline function myvector_mode_jacobian(f, x::StaticArray)
    T = typeof(ForwardDiff.Tag(f, eltype(x)))
    result = ForwardDiff.static_dual_eval(T, f, x)
    J = ForwardDiff.extract_jacobian(T, result, x)
    f = myvalue(result)
    return f, J
end

res = myjacobian( x-> SInput(fs, x), input)

ForwardDiff.extract_jacobian(T, ForwardDiff.static_dual_eval(T, f, x), x)

SVector{N, ForwardDiff.Dual{Tag, T, N}}

@generated function myvalue(result::SVector{N, ForwardDiff.Dual{Tag, T, N}}) where {N,T,Tag}
    quote
        Base.Cartesian.@nexprs $N i -> result_i = result[i].value
        Base.Cartesian.@ncall  $N SVector{$N, $T} result
    end
end


myvalue(res)