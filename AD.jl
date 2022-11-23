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

# Check derivatives 
vec1 = [c1 c2 c3 c4 c5 p1 p2 p3]
args = (T=900.0, d=100e-6, τII_old=1e6, dt=1e8)
εII, τII = 2e-15, 2e6

v = c3   # problem with c3 (elasticity) and c4 (that has || elements) 

@inline function derivative_all(f::F, x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    res = f(ForwardDiff.Dual{T}(x, one(x)))
    return res.value, res.partials.values[1]
end

# function λ_εII(vi, τII, args)
#     # trick with nested lambdas so that it works with @generated
#     # It'd work with one lambda, but for some reason this is
#     # slightly more optimal
#     @inline λ(vi, τII) = begin
#         _λ = τII -> compute_εII(vi, τII, args)
#         derivative_all(_λ, τII)
#     end

#     λ(vi, τII)
# end

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
compute_εII_AD_all(v::CompositeRheology, τII, args) = compute_εII_AD_all(v.elements, τII, args)

compute_τII_AD_all(v, εII, args) = compute_II_AD_all(v, εII, λ_τII, args)
compute_τII_AD_all(v::CompositeRheology, εII, args) = compute_εII_AD_all(v.elements, εII, args)


function compute_τII_AD(
    v::CompositeRheology{T,N,
        0,is_parallel,
        0,is_plastic,
        Nvol,is_vol,
        false},
    εII::_T, 
    args; 
    tol=1e-6
) where {N, T, _T, is_parallel, is_plastic, Nvol, is_vol}

    # Initial guess
    τII = compute_τII_harmonic(v, εII, args)

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
        εII_n, dεII_dτII_n = compute_εII_AD_all(v, τII, args)

        τII = muladd(εII - εII_n, inv(dεII_dτII_n), τII)

        ϵ = abs(τII - τII_prev) * inv(τII)
        τII_prev = τII

    end

    return τII
end


# Define a range of rheological components
v1 = SetDiffusionCreep("Dry Anorthite | Rybacki et al. (2006)")
v2 = SetDislocationCreep("Dry Anorthite | Rybacki et al. (2006)")
e1 = ConstantElasticity()           # elasticity
c = CompositeRheology(v1, v2, e1)   # with elasticity

compute_τII_AD(c, εII, args) == compute_τII(c, εII, args)

@btime compute_τII_AD($c, $εII, $args)
@btime compute_τII($c, $εII, $args)


# ProfileCanvas.@profview for i in 1:1000000
#     compute_τII_AD(v, εII, args)
# end

# @btime compute_τII_harmonic($v, $εII, $args)