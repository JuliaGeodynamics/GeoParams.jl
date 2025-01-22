# Package extension for adding GLMakie-based features to GeoParams.jl
module GeoParamsGLMakieExt

using GeoParams

# We do not check `isdefined(Base, :get_extension)` as recommended since
# Julia v1.9.0 does not load package extensions when their dependency is
# loaded from the main environment.
if VERSION >= v"1.9.1"
    using GLMakie
else
    using ..GLMakie
end

include("../src/Plotting/Plotting.jl")
include("../src/Plotting/StrengthEnvelope.jl")

end # module
