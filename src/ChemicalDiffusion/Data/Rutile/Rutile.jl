module Rutile

using GeoParams

include("Elements/Hf.jl")
include("Elements/Zr.jl")

"""
    chemical_diffusion_list(search::String="")

List all available chemical diffusion data for melt.
Includes an argument to search for a specific term, i.e. an element ("La") or an author.
"""
function chemical_diffusion_list(search::String="")
    m = @__MODULE__
    all = string.(names(m; all = true, imported = true))
    Rt = Symbol.(filter(x -> startswith(x, "Rt") && contains(x, search), all))
    return [getfield(m, fun) for fun in Rt]
end


end
