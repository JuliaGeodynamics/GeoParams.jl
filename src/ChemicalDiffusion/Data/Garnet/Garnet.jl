module Garnet

using GeoParams

include("Elements/Fe.jl")
include("Elements/Mg.jl")
include("Elements/Mn.jl")
include("Elements/REE.jl")

"""
    chemical_diffusion_list(search::String="")

List all available chemical diffusion data for melt.
Includes an argument to search for a specific term, i.e. an element ("La") or an author.
"""
function chemical_diffusion_list(search::String="")
    m = @__MODULE__
    all = string.(names(m; all = true, imported = true))
    Grt = Symbol.(filter(x -> startswith(x, "Grt") && contains(x, search), all))
    return [getfield(m, fun) for fun in Grt]
end


end
