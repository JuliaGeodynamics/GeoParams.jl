using Test

files = readdir(@__DIR__)
test_files = filter(startswith("test_"), files)

for f in test_files
    include(f)
end

# include("test_GeoUnits.jl")
# include("test_MaterialParameters.jl")
# include("test_CreepLaw.jl")
# include("test_CompositeRheologies.jl")
# include("test_Density.jl")
# include("test_Energy.jl")
# include("test_MeltingParam.jl")
# include("test_TAS_classification.jl")
# include("test_ZirconAge.jl")
# include("test_SeismicVelocity.jl")
# include("test_DiffusionCreep.jl")
# include("test_DislocationCreep.jl")
# include("test_Plasticity.jl")
