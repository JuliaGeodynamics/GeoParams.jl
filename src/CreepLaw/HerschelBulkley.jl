
struct HerschelBulkley{T, U1, U2, U3, U4} <: AbstractCreepLaw{T}
    n::T  # shear thinning exponent
    η0::GeoUnit{T, U1} # cutoff viscosity
    σ0::GeoUnit{T, U2} # critical stress
    Kv::GeoUnit{T, U3} # consistency
    A::T
    B::GeoUnit{T, U4}

    function HerschelBulkley(;
            n = (29//50),
            η0 = 1e8Pa*s,
            σ0 = 0.0100Pa,
            Kv = 0.140Pa*s^n,
            A = 1.6927,
            B = -0.0257/K,
        )
        # Convert to GeoUnits
        # nU = convert(GeoUnit, rat2float(n))
        η0U = convert(GeoUnit, η0)
        σ0U = convert(GeoUnit, σ0)
        KvU = convert(GeoUnit, Kv)
        # AU = convert(GeoUnit, A)
        BU = convert(GeoUnit, B)
        # Extract struct types
        T = typeof(η0U).types[1]
        U1 = typeof(η0U).types[2]
        U2 = typeof(σ0U).types[2]
        U3 = typeof(KvU).types[2]
        U4 = typeof(BU).types[2]
        # Create struct
        return new{T, U1, U2, U3, U4}(
            n, η0U, σ0U, KvU, A, BU
        )
    end

    function HerschelBulkley(n, η0, σ0, Kv, A, B)
        return HerschelBulkley(;
            n = n, η0 = η0, σ0 = σ0, Kv = Kv, A = A, B = B
        )
    end
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
