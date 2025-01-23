module Garnet

using GeoParams

include("Fe.jl")
include("Mg.jl")
include("Mn.jl")
include("REE.jl")

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
