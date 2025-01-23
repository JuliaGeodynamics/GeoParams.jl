module Rutile

using GeoParams

include("Hf.jl")
include("Zr.jl")

"""
    chemical_diffusion_list()

List all available chemical diffusion data for rutile.
"""
function chemical_diffusion_list()
    m = @__MODULE__
    all = string.(names(m; all = true, imported = true))
    Rt = Symbol.(filter(x -> startswith(x, "Rt"), all))
    return [getfield(m, fun) for fun in Rt]
end


end
