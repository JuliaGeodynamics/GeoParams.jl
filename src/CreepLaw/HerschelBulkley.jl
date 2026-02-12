

struct HerschelBulkley{T, U1, U2, U3} <: AbstractCreepLaw{T}
    n::T  # shear thinning exponent
    η0::GeoUnit{T, U1} # "rigid" viscosity
    τ0::GeoUnit{T, U2} # critical stress
    ηr::GeoUnit{T, U1} # reference viscosity at the critical strain rate, which is given by 0.5*τ0/η0 and the critical temperature
    Q::GeoUnit{T, U3} # temperature dependence of ηr, activation energy divided by R, unit is K
    Tr::GeoUnit{T,U3} # reference temperature
    function HerschelBulkley(;
        n = 3,
        η0 = 1e24Pa*s,
        τ0 = 100e6Pa,
        ηr = 1e20Pa*s,
        Q  = 0K,
        Tr = 1600K,
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
        U2 = typeof(σ0U).types[2]
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


# function when eII is given 
function HB_rheology(εII,η0,τ0,εr,ηr,n,ηl)
    # compute viscosity using the Souza-Mendes approach
     η = (1.0 - exp(-2.0*η0*εII/τ0)) * (0.5*τ0/εII  + ηr*(εII./εr)^(1/n - 1)) + ηl
     return η
end



# @inline function HerschelBulkley(
#     a::DislocationCreep,
#     EpsII;
#     T = one(precision(a)),
#     args...,
# )
#     n, η0, σ0, Kv, A, B= if EpsII isa Quantity
#         @unpack_units n, η0, σ0, Kv, A, B = a
#         n, η0, σ0, Kv, A, B
#     else
#         @unpack_val n, η0, σ0, Kv, A, B = a
#         n, η0, σ0, Kv, A, B
#     end

#     KvT = Kv * A * exp(B * (T - 273.15))
#     τ = @pow (1 - exp(-η0 * EpsII / σ0)) * (σ0 + KvT * (εII^n))
   
# end

# @inline function compute_viscosity(
#     a::HerschelBulkley; T = 0.0, kwargs...
#     )
#     n, η0, σ0, Kv, A, B= if EpsII isa Quantity
#         @unpack_units n, η0, σ0, Kv, A, B = a
#         n, η0, σ0, Kv, A, B
#     else
#         @unpack_val n, η0, σ0, Kv, A, B = a
#         n, η0, σ0, Kv, A, B
#     end

#     KvT = Kv * A * exp(B * (T - 273.15))

#     η = (1 - exp(-η0 *εII / σ0)) * (σ0 + KvT * (εII^n))
#     return  η
# end