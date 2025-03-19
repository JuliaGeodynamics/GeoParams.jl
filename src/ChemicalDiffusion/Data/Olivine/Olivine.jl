module Olivine

using GeoParams
using Unitful

include("Elements/Mg.jl")

"""
    chemical_diffusion_list(search::String="")

List all available chemical diffusion data for rutile.
Includes an argument to search for a specific term, i.e. an element ("La") or an author.
"""
function chemical_diffusion_list(search::String = "")
    m = @__MODULE__
    all = string.(names(m; all = true, imported = true))
    Ol = Symbol.(filter(x -> startswith(x, "Ol") && contains(x, search), all))
    return [getfield(m, fun) for fun in Ol]
end


end
