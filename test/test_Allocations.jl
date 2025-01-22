# this tests allocations for a range of scenarios
using Test
using GeoParams
import GeoParams.Dislocation

@testset "Allocations" begin

    @static if VERSION ≥ v"1.8.0"
        # Numerical parameters
        Ncy = 100
        # Allocate arrays
        ε0 = 1.0e-4
        Tv = rand(Ncy + 1)
        ε̇ii = fill(ε0, Ncy + 1)
        η = zeros(Ncy + 1)

        #issue 129, part 2 (was allocating on Julia 1.9 but not on 1.10)
        function main1(Tv, ε̇ii, η)
            # Unit system
            CharDim = SI_units(length = 1000m, temperature = 1000C, stress = 1.0e7Pa, viscosity = 1.0e20Pas)

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
            flow_nd = Transform_DislocationCreep(flow_nd0, CharDim)

            # Setup up viscosity model
            a = @allocated begin
                for i in eachindex(ε̇ii)
                    η[i] = compute_viscosity_εII(flow_nd, ε̇ii[i], (; T = Tv[i]))
                end
            end

            return a
        end

        @test main1(Tv, ε̇ii, η) == 0

        function main2(Tv, ε̇ii, η)
            # Unit system
            CharDim = SI_units(length = 1000m, temperature = 1000C, stress = 1.0e7Pa, viscosity = 1.0e20Pas)
            # Configure viscosity model
            flow_nd = SetDislocationCreep(Dislocation.diabase_Caristan_1982, CharDim)

            # Setup up viscosity model
            a = @allocated begin
                for i in eachindex(ε̇ii)
                    η[i] = compute_viscosity_εII(flow_nd, ε̇ii[i], (; T = Tv[i]))
                end
            end

            return a
        end

        @test main2(Tv, ε̇ii, η) == 0
    end

end
