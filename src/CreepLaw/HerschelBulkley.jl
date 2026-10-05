export HerschelBulkley,
    compute_εII!,
    compute_εII,
    compute_τII!,
    compute_τII,
    compute_hb_viscosity_εII,
    compute_hb_viscosity_τII

"""
    HerschelBulkley(; n=3.0, η0=1e24Pa*s, τ0=100e6Pa, ηr=1e20Pa*s, Q, Tr)

A Herschel-Bulkley creep law: a temperature-dependent viscoplastic fluid with a yield stress `τ0`,
"rigid" viscosity `η0`, shear-thinning exponent `n`, and reference viscosity `ηr` at reference
temperature `Tr`, with `Q` the activation energy divided by the gas constant.
"""
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
    return compute_hb_εII(a, TauII; T)
end

function compute_εII(a::HerschelBulkley, TauII::Quantity; T = 1K, kwargs...)
    return compute_hb_εII(a, TauII; T)
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
    # in x = εII / εr: the form τ0 / (2εII) has the derivative τ0 / (2εII²), which overflows Float32
    x = εII / εr
    η = @pow (-expm1(-x)) * (η0 / x + ηT * x^(inv(n) - 1))
    return η
end


"""
    compute_hb_viscosity_τII(a::HerschelBulkley, TauII; T = one(precision(a)), kwargs...)

function to compute the viscosity if TauII is given
"""
@inline function compute_hb_viscosity_τII(v::HerschelBulkley, τII; T = one(precision(v)), kwargs...)
    return compute_hb_viscosity_εII(v, compute_hb_εII(v, τII; T); T)
end

"""
    compute_hb_εII(a::HerschelBulkley, TauII; T = one(precision(a)), kwargs...)

Strain rate for a given stress, by Newton iteration.
"""
@inline function compute_hb_εII(v::HerschelBulkley, τII; T = one(precision(v)), kwargs...)
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

    # Solved for x = εII / εr, in which every term is O(1): in εII the Newton derivatives
    # reach 1e42 for laboratory parameters and overflow Float32. The ratios are unitless,
    # so ForwardDiff never sees Quantity{Dual} types.
    τ̃ = ustrip(τII / τ0)
    ηratio = ustrip(ηT / η0)

    # residual 2ηεII / τ0 - τII / τ0, with η = (1 - exp(-x)) (η0 / x + ηT x^(1/n - 1))
    fres(x) = (-expm1(-x)) * (1 + ηratio * x^inv(n)) - τ̃

    # initial guess: below yield η ≈ η0, above it the power-law branch 1 + ηratio x^(1/n) ≈ τ̃,
    # offset by τ̃ so that the guess at the yield stress is not x = 0, where the residual is singular
    x = τ̃ < 1 ? τ̃ * one(ηratio) : τ̃ + ((τ̃ - 1) / ηratio)^n

    # the residual cancels to the last bits of τ̃ near the root, so the step stalls a few
    # eps above zero; √eps is the step whose quadratic convergence puts the answer at that level
    tol = sqrt(eps(Tc))
    it_max = 100
    for _ in 1:it_max
        f, dfdx = value_and_partial(fres, x)
        Δx = f / dfdx
        x -= Δx
        if abs(Δx) < tol * abs(x)
            # one more step: the derivatives carried by a Dual x converge one iteration
            # behind its value
            f, dfdx = value_and_partial(fres, x)
            return (x - f / dfdx) * εr
        end
    end
    return error("compute_hb_εII: iterations did not converge for τII=$τII after $it_max iterations, tolerance $tol")
end

# print info
function show(io::IO, g::HerschelBulkley)
    return print(
        io,
        "Hershel Bulkley viscosity: η0=$(Value(g.η0)), τ0=$(Value(g.τ0)), ηr=$(Value(g.ηr)), n=$(g.n), Q=$(Value(g.Q)), Tr=$(Value(g.Tr))",
    )
end
