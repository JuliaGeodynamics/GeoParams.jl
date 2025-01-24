module Garnet

using GeoParams

include("La.jl")
include("Ce.jl")
include("Pr.jl")
include("Nd.jl")
include("Sm.jl")
include("Gd.jl")
include("Tb.jl")
include("Dy.jl")
include("Ho.jl")
include("Er.jl")
include("Yb.jl")
include("Lu.jl")

"""
    chemical_diffusion_list()

List all available chemical diffusion data for garnet.
"""
function chemical_diffusion_list()
    m = @__MODULE__
    all = string.(names(m; all = true, imported = true))
    Grt = Symbol.(filter(x -> startswith(x, "Grt"), all))
    return [getfield(m, fun) for fun in Grt]
end


end
