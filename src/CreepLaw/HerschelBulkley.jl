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
    EpsII = 0.5 * TauII / η
    return EpsII
end

function compute_εII(a::HerschelBulkley, TauII::Quantity; T = 1K, kwargs...)
    η = compute_hb_viscosity_τII(a, TauII; T = T)
    EpsII = 0.5 * TauII / η
    return EpsII
end

"""
    compute_εII!(EpsII::AbstractArray{T, N}, a::HerschelBulkley, TauII::AbstractArray{T, N}; T = ones(size(TauII)), kwargs...)

In-place function for the second invariant of the strain rate for Herschel-Bulkley rheology.
"""
function compute_εII!(
        EpsII::AbstractArray{_T, N},
        a::HerschelBulkley,
        TauII::AbstractArray{_T, N};
        T = ones(size(TauII))::AbstractArray{_T, N},
        kwargs...,
    ) where {_T, N}
    @inbounds for i in eachindex(EpsII)
        EpsII[i] = compute_εII(a, TauII[i]; T = T[i])
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
    compute_τII!(TauII::AbstractArray{T, N}, a::HerschelBulkley, EpsII::AbstractArray{T, N}; T = ones(size(EpsII)), kwargs...)

In-place function for the second invariant of the stress for Herschel-Bulkley rheology.
"""
function compute_τII!(
        TauII::AbstractArray{_T, N},
        a::HerschelBulkley,
        EpsII::AbstractArray{_T, N};
        T = ones(size(EpsII)),
        kwargs...,
    ) where {_T, N}
    @inbounds for i in eachindex(TauII)
        TauII[i] = compute_τII(a, EpsII[i]; T = T[i])
    end

    return nothing
end


"""
    compute_hb_viscosity_εII(a::HerschelBulkley, EpsII; T = one(precision(a)), kwargs...)

function to compute the viscosity if EpsII is given
"""
@inline function compute_hb_viscosity_εII(v::HerschelBulkley, εII; T = 1.0, kwargs...)
    η0, τ0, ηr, Q, Tr = if εII isa Quantity
        @unpack_units η0, τ0, ηr, Q, Tr = v
        η0, τ0, ηr, Q, Tr
    else
        @unpack_val η0, τ0, ηr, Q, Tr = v
        η0, τ0, ηr, Q, Tr
    end
    (; n) = v

    ηT = ηr * exp(Q * (1 / T - 1 / Tr)) # temperature dependence
    εr = 0.5 * τ0 / η0 # strain rate at which the Bingham yield stress is reached, this is defined as the reference strain rate
    η = @pow (1.0 - exp(-2.0 * η0 * εII / τ0)) * (0.5 * τ0 / εII + ηT * (εII / εr)^(one(n) / n - 1))
    return η
end


"""
compute_hb_viscosity_τII(a::HerschelBulkley, EpsII; T = one(precision(a)), kwargs...)

function to compute the viscosity if TauII is given
"""

@inline function compute_hb_viscosity_τII(v::HerschelBulkley, τII; T = one(precision(v)), kwargs...)

    η0, τ0, ηr, Q, Tr = if τII isa Quantity
        @unpack_units η0, τ0, ηr, Q, Tr = v
        η0, τ0, ηr, Q, Tr
    else
        @unpack_val η0, τ0, ηr, Q, Tr = v
        T = ustrip(T)
        η0, τ0, ηr, Q, Tr
    end
    (; n) = v

    ηT = ηr * exp(Q * (1 / T - 1 / Tr))
    εr = 0.5 * τ0 / η0

    # initial guess
    η = if τII < τ0
        η0
    elseif τII == τ0
        (1 - exp(-one(η0))) * (η0 + ηT)
    else
        @pow (ηT / (1 - τ0 / τII))^n * (τII / (2 * εr))^(1 - n)
    end

    εII = 0.5 * τII / η
    εII_unit = τII isa Quantity ? unit(εII) : one(εII)

    # strip ALL quantities to plain floats before Newton iteration
    # so that ForwardDiff never sees Quantity{Dual} types
    τII_s = ustrip(τII)
    η0_s = ustrip(η0)
    τ0_s = ustrip(τ0)
    ηT_s = ustrip(ηT)
    εr_s = ustrip(εr)
    εII_s = ustrip(εII)

    tol = 1.0e-10
    it_max = 100

    for _ in 1:it_max
        f, dfdε = value_and_partial(
            ε -> fres_hb(ε, τII_s, η0_s, τ0_s, ηT_s, n, εr_s),  # all plain floats
            εII_s
        )
        Δε = f / dfdε
        εII_s -= Δε
        abs(Δε) / (abs(εII_s) + eps(typeof(εII_s))) < tol && break
    end

    εII = εII_s * εII_unit
    η = @pow (1.0 - exp(-2.0 * η0 * εII / τ0)) * (0.5 * τ0 / εII + ηT * (εII / εr)^(one(n) / n - 1))

    return η
end

# non-dimensional residual
function fres_hb(εII::Number, τII::Number, η0::Number, τ0::Number, ηT::Number, n::Number, εr::Number)
    η = @pow (1.0 - exp(-2.0 * η0 * εII / τ0)) * (0.5 * τ0 / εII + ηT * (εII / εr)^(one(n) / n - 1))
    return 2.0 * η * εII - τII
end


# print info
function show(io::IO, g::HerschelBulkley)
    return print(
        io,
        "Hershel Bulkley viscosity: η0=$(Value(g.η0)), τ0=$(Value(g.τ0)), ηr=$(Value(g.ηr)), n=$(g.n), Q=$(Value(g.Q)), Tr=$(Value(g.Tr))",
    )
end
