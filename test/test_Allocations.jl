# this tests allocations for a range of scenarious
using Test
using GeoParams

@testset "Allocations" begin

    #issue 129, part 2 (was allocating on juli 1.9 but not on 1.10)
    function main1()
        # Unit system
        CharDim    = SI_units(length=1000m, temperature=1000C, stress=1e7Pa, viscosity=1e20Pas)
        # Numerical parameters
        Ncy        = 100
        # Allocate arrays
        ε0         =    1e-4
        Tv         =    rand(Ncy+1)
        ε̇ii        = ε0*ones(Ncy+1) 
        η          =   zeros(Ncy+1)
        
        # Configure viscosity model
        flow_nd0 = DislocationCreep(;
            Name = "Diabase | Caristan (1982)",
            n = 3.05NoUnits,
            A = 6.0e-2MPa^(-61 // 20) / s,
            E = 276kJ / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
        flow_nd  = Transform_DislocationCreep(flow_nd0, CharDim)
    
        # Setup up viscosity model
        a = @allocated begin
            for i in eachindex(ε̇ii)
                η[i] = compute_viscosity_εII(flow_nd, ε̇ii[i], (;T=Tv[i]))
            end
        end

        return a
    end

    a = main1()
    @test a==0


"""
    entry = extract_database_entry(Name::String, database::NTuple{N,AbstractCreepLaw})

Extracts an entry from a creeplaw database
"""
function extract_database_entry(Name::String, database::NTuple{N,AbstractCreepLaw}) where {N}
    names = extract_database_names(database)
    found = false
    name_pad = rpad(Name, length(names[1]))
    entry = database[1]
    for i = 1:N
        if (names[i] .==  name_pad)
            entry = database[i]
            found = true
        end
    end
    if !found; error("Unknown database entry: $Name"); end

   return entry
end

"""
    names = extract_database_names(database::Tuple)
Returns a vector with all `names` in the `database`
"""
function extract_database_names(database::Tuple)
    return [String(collect(f.Name)) for f in database]
end

include("./src/CreepLaw/Data_deprecated/DislocationCreep.jl")

function main()

    flowlaws = Tuple([i[2][1] for i in DislocationCreep_info]);

    flowlaws=NTuple()
    for i in DislocationCreep_info
        flowlaws = merge(flowlaws..., i[2][1])
    end
    #=
    flowlaws = ( 
                DislocationCreep(;
                    Name = "Dry Olivine | Hirth & Kohlstedt (2003)",
                    n = 3.5NoUnits,
                    r = 0.0NoUnits,
                    A = 1.1e5MPa^(-7 // 2) / s,
                    E = 530.0kJ / mol,
                    V = 14e-6m^3 / mol,
                    Apparatus = AxialCompression,
                ),
                DislocationCreep(;
                    Name = "2. Wet Olivine | Hirth & Kohlstedt (2003)",
                    n = 3.0NoUnits,
                    A = 1600MPa^(-3) / s,
                    E = 520.0kJ / mol,
                    V = 22.0m^3 / mol,
                    r = 1.2NoUnits,
                    Apparatus = AxialCompression,
                ),
                DislocationCreep(;
                    Name = "Diabase | Caristan (1982)",
                    n = 3.05NoUnits,
                    A = 6.0e-2MPa^(-61 // 20) / s,
                    E = 276kJ / mol,
                    V = 0m^3 / mol,
                    r = 0NoUnits,
                    Apparatus = AxialCompression,
                )
    )
=#

    # Unit system
    CharDim    = SI_units(length=1000m, temperature=1000C, stress=1e7Pa, viscosity=1e20Pas)
    # Numerical parameters
    Ncy        = 100
    # Allocate arrays
    ε0         =    1e-4
    Tv         =    rand(Ncy+1)
    ε̇ii        = ε0*ones(Ncy+1) 
    η          =   zeros(Ncy+1)

    # Configure viscosity model
    flow_nd0 = extract_database_entry("Diabase | Caristan (1982)", flowlaws)
    flow_nd  = Transform_DislocationCreep(flow_nd0, CharDim)

    # Setup up viscosity model
    a = @allocated begin
        for i in eachindex(ε̇ii)
            η[i] = compute_viscosity_εII(flow_nd, ε̇ii[i], (;T=Tv[i]))
        end
    end
    return a
end

end