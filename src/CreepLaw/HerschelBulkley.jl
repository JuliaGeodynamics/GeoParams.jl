export HerschelBulkley,
    compute_εII!,
    compute_εII,
    compute_τII!,
    compute_τII,
    compute_hb_viscosity_εII,
    compute_hb_viscosity_τII

struct HerschelBulkley{T, U1, U2, U3} <: AbstractCreepLaw{T}
    n::T  # shear thinning exponent
    η0::GeoUnit{T, U1} # "rigid" viscosity
    τ0::GeoUnit{T, U2} # critical stress
    ηr::GeoUnit{T, U1} # reference viscosity at the critical strain rate, which is given by 0.5*τ0/η0 and the critical temperature
    Q::GeoUnit{T, U3} # temperature dependence of ηr, activation energy divided by R, unit is K
    Tr::GeoUnit{T, U3} # reference temperature
    function HerschelBulkley(;
            n = 3.0,
            η0 = 1.0e24Pa * s,
            τ0 = 100.0e6Pa,
            ηr = 1.0e20Pa * s,
            Q = 0.0K,
            Tr = 1273K,
        )
        # Convert to GeoUnits
        η0U = convert(GeoUnit, η0)
        τ0U = convert(GeoUnit, τ0)
        ηrU = convert(GeoUnit, ηr)
        QU = convert(GeoUnit, Q)
        Tr = convert(GeoUnit, Tr)
        # Extract struct types
        T = typeof(η0U).types[1]
        U1 = typeof(η0U).types[2]
        U2 = typeof(τ0U).types[2]
        U3 = typeof(Tr).types[2]
        # Create struct
        return new{T, U1, U2, U3}(
            n, η0U, τ0U, ηrU, QU, Tr
        )
    end

    function HerschelBulkley(n, η0, τ0, ηr, Q, Tr)
        return HerschelBulkley(; n = n, η0 = η0, τ0 = τ0, ηr = ηr, Q = Q, Tr = Tr)
    end
end

function compute_εII(a::HerschelBulkley, TauII; T = one(precision(a)), kwargs...)
    η = compute_hb_viscosity_τII(a, TauII; T = T)
    EpsII = TauII / (2 * η)
    return EpsII
end

function compute_εII(a::HerschelBulkley, TauII::Quantity; T = 1K, kwargs...)
    η = compute_hb_viscosity_τII(a, TauII; T = T)
    EpsII = TauII / (2 * η)
    return EpsII
end

"""
    compute_εII!(EpsII::AbstractArray{_T, N}, a::HerschelBulkley, TauII::AbstractArray; T = one(_T), kwargs...)

In-place function for the second invariant of the strain rate for Herschel-Bulkley rheology.

`T` may be a scalar, applied to every element, or an array indexed alongside `TauII`.
"""
function compute_εII!(
        EpsII::AbstractArray{_T, N},
        a::HerschelBulkley,
        TauII::AbstractArray;
        T = one(_T),
        kwargs...,
    ) where {_T, N}
    for i in each_argument_index(EpsII, TauII, T)
        EpsII[i] = compute_εII(a, convert_precision(_T, TauII[i]); T = argument_at(T, i))
    end

    return nothing
end

"""
    compute_τII(a::HerschelBulkley, EpsII; T = one(precision(a)), kwargs...)

"""
function compute_τII(a::HerschelBulkley, EpsII; T = one(precision(a)), kwargs...)
    η = compute_hb_viscosity_εII(a, EpsII; T = T)
    TauII = 2 * EpsII * η
    return TauII
end

function compute_τII(a::HerschelBulkley, EpsII::Quantity; T = 1K, kwargs...)
    η = compute_hb_viscosity_εII(a, EpsII; T = T)
    TauII = 2 * EpsII * η
    return TauII
end

"""
    compute_τII!(TauII::AbstractArray{_T, N}, a::HerschelBulkley, EpsII::AbstractArray; T = one(_T), kwargs...)

In-place function for the second invariant of the stress for Herschel-Bulkley rheology.

`T` may be a scalar, applied to every element, or an array indexed alongside `EpsII`.
"""
function compute_τII!(
        TauII::AbstractArray{_T, N},
        a::HerschelBulkley,
        EpsII::AbstractArray;
        T = one(_T),
        kwargs...,
    ) where {_T, N}
    for i in each_argument_index(TauII, EpsII, T)
        TauII[i] = compute_τII(a, convert_precision(_T, EpsII[i]); T = argument_at(T, i))
    end

    return nothing
end


"""
    compute_hb_viscosity_εII(a::HerschelBulkley, EpsII; T = one(precision(a)), kwargs...)

function to compute the viscosity if EpsII is given
"""
@inline function compute_hb_viscosity_εII(v::HerschelBulkley, εII; T = 1.0, kwargs...)
    Tc = precision_of(εII)
    T = convert_precision(Tc, T)
    η0, τ0, ηr, Q, Tr = if εII isa Quantity
        @unpack_units Tc η0, τ0, ηr, Q, Tr = v
        η0, τ0, ηr, Q, Tr
    else
        @unpack_val Tc η0, τ0, ηr, Q, Tr = v
        η0, τ0, ηr, Q, Tr
    end
    n = convert_precision(Tc, v.n)

    ηT = ηr * exp(Q * (1 / T - 1 / Tr)) # temperature dependence
    εr = τ0 / (2 * η0) # strain rate at which the Bingham yield stress is reached, this is defined as the reference strain rate
    η = @pow (1 - exp(-2 * η0 * εII / τ0)) * (τ0 / (2 * εII) + ηT * (εII / εr)^(one(n) / n - 1))
    return η
end


"""
compute_hb_viscosity_τII(a::HerschelBulkley, EpsII; T = one(precision(a)), kwargs...)

function to compute the viscosity if TauII is given
"""

@inline function compute_hb_viscosity_τII(v::HerschelBulkley, τII; T = one(precision(v)), kwargs...)
    Tc = precision_of(τII)
    T = convert_precision(Tc, T)

    η0, τ0, ηr, Q, Tr = if τII isa Quantity
        @unpack_units Tc η0, τ0, ηr, Q, Tr = v
        η0, τ0, ηr, Q, Tr
    else
        @unpack_val Tc η0, τ0, ηr, Q, Tr = v
        T = ustrip(T)
        η0, τ0, ηr, Q, Tr
    end
    n = convert_precision(Tc, v.n)

    ηT = ηr * exp(Q * (1 / T - 1 / Tr))
    εr = τ0 / (2 * η0)

    # strip ALL quantities to plain floats before the Newton iteration, so that
    # ForwardDiff never sees Quantity{Dual} types, and so that the initial guess
    # can raise them to a fractional power: Unitful represents the dimensions of
    # such a power only approximately, and the result no longer reduces to Pa s
    τII_s = ustrip(τII)
    η0_s = ustrip(η0)
    τ0_s = ustrip(τ0)
    ηT_s = ustrip(ηT)
    εr_s = ustrip(εr)

    # initial guess
    η_s = if τII_s < τ0_s
        η0_s
    elseif τII_s == τ0_s
        (1 - exp(-one(η0_s))) * (η0_s + ηT_s)
    else
        # raised as a whole: for laboratory parameters the separate factors reach
        # 1e60 and 1e-49, outside Float32 range, where the product they form is an
        # ordinary viscosity
        @pow (ηT_s * (τII_s / (2 * εr_s))^(inv(n) - 1) / (1 - τ0_s / τII_s))^n
    end

    εII_s = τII_s / (2 * η_s)
    εII_unit = τII isa Quantity ? unit(τII / η0) : one(εII_s)

    # the residual 2ηε - τII cancels to the last bits of τII near the root, so the
    # step stalls a few eps above zero; √eps is the step whose quadratic
    # convergence puts the answer at that level
    tol = sqrt(eps(Tc))
    it_max = 100
    res = one(tol)

    for _ in 1:it_max
        f, dfdε = value_and_partial(
            ε -> fres_hb(ε, τII_s, η0_s, τ0_s, ηT_s, n, εr_s),  # all plain floats
            εII_s
        )
        Δε = f / dfdε
        εII_s -= Δε
        res = abs(Δε) / (abs(εII_s) + eps(typeof(εII_s)))
        res < tol && break
    end
    res < tol || error(
        "compute_hb_viscosity_τII: iterations did not converge for τII=$τII: relative step $res after $it_max iterations, tolerance $tol"
    )

    εII = εII_s * εII_unit
    η = @pow (1 - exp(-2 * η0 * εII / τ0)) * (τ0 / (2 * εII) + ηT * (εII / εr)^(one(n) / n - 1))

    return η
end

# non-dimensional residual
function fres_hb(εII::Number, τII::Number, η0::Number, τ0::Number, ηT::Number, n::Number, εr::Number)
    η = @pow (1 - exp(-2 * η0 * εII / τ0)) * (τ0 / (2 * εII) + ηT * (εII / εr)^(one(n) / n - 1))
    return 2 * η * εII - τII
end


# print info
function show(io::IO, g::HerschelBulkley)
    return print(
        io,
        "Hershel Bulkley viscosity: η0=$(Value(g.η0)), τ0=$(Value(g.τ0)), ηr=$(Value(g.ηr)), n=$(g.n), Q=$(Value(g.Q)), Tr=$(Value(g.Tr))",
    )
end
