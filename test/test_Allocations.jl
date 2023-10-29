# this tests allocations for a range of scenarious
using Test
using GeoParams

@testset "Allocations" begin

    #issue 129
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


end