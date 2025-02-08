# Package extension for adding Makie-based features to GeoParams.jl
module GeoParamsMakieExt

using GeoParams

# We do not check `isdefined(Base, :get_extension)` as recommended since
# Julia v1.9.0 does not load package extensions when their dependency is
# loaded from the main environment.
if VERSION >= v"1.9.1"
    using Makie
else
    using ..Makie
end

include("../src/Plotting/Plotting.jl")
include("../src/Plotting/StrengthEnvelope.jl")

end # module
