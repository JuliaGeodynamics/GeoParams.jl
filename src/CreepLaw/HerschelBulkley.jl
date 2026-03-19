

struct HerschelBulkley{T, U1, U2, U3} <: AbstractCreepLaw{T}
    n::T  # shear thinning exponent
    η0::GeoUnit{T, U1} # "rigid" viscosity
    τ0::GeoUnit{T, U2} # critical stress
    ηr::GeoUnit{T, U1} # reference viscosity at the critical strain rate, which is given by 0.5*τ0/η0 and the critical temperature
    Q::GeoUnit{T, U3} # temperature dependence of ηr, activation energy divided by R, unit is K
    Tr::GeoUnit{T,U3} # reference temperature
    function HerschelBulkley(;
        n = 3.0,
        η0 = 1e24Pa*s,
        τ0 = 100e6Pa,
        ηr = 1e20Pa*s,
        Q  = 0.0K,
        Tr = 1273K,
        )
        # Convert to GeoUnits
        η0U = convert(GeoUnit, η0)
        τ0U = convert(GeoUnit, τ0)
        ηrU = convert(GeoUnit, ηr)
        QU  = convert(GeoUnit,Q)
        Tr  = convert(GeoUnit,Tr)
        # Extract struct types
        T = typeof(η0U).types[1]
        U1 = typeof(η0U).types[2]
        U2 = typeof(τ0U).types[2]
        U3 = typeof(Tr).types[2]
        # Create struct
        return new{T, U1, U2, U3}(
            n, η0U,τ0U,ηrU,QU, Tr
        )
    end

    function HerschelBulkley(n, η0,τ0,ηr,Q,Tr)
        return HerschelBulkley(;n = n, η0 = η0,τ0 = τ0,ηr = ηr,Q = Q,Tr = Tr)
    end
end

function compute_εII(a::HershelBulkley, TauII; T = one(precision(a)), kwargs...)
    η = compute_viscosity_τII(a, TauII; T = T)
    EpsII = 0.5 * TauII / η
    return EpsII
end

function compute_εII(a::HershelBulkley, TauII; T = 1K, kwargs...)
    η = compute_viscosity_τII(a, TauII; T = T)
    EpsII = 0.5 * TauII / η
    return EpsII
end


function compute_τII(a::HershelBulkley, EpsII; T = one(precision(a)), kwargs...)
    η = compute_viscosity_εII(a, EpsII; T = T)
    TauII = 2 * EpsII * η
    return TauII
end

function compute_τII(a::HershelBulkley, EpsII; T = 1K, kwargs...)
    η = compute_viscosity_εII(a, EpsII; T = T)
    TauII = 2 * EpsII * η
    return TauII
end

# print info
function show(io::IO, g::HershelBulkley)
    return print(
        io,
        "Hershel Bulkley viscosity: η0=$(Value(g.η0)), τ0=$(Value(g.τ0)), ηr=$(Value(g.ηr)), n=$(Value(g.n)), Q=$(Value(g.Q)), Tr=$(Value(g.Tr))",
    )
end