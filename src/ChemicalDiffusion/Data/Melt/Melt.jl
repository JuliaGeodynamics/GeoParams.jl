module Melt

using GeoParams


include("Elements/Sc.jl")
include("Elements/V.jl")
include("Elements/Y.jl")
include("Elements/Zr.jl")
include("Elements/Hf.jl")
include("Elements/Nb.jl")
include("Elements/La.jl")
include("Elements/Ce.jl")
include("Elements/Pr.jl")
include("Elements/Nd.jl")
include("Elements/Sm.jl")
include("Elements/Eu.jl")
include("Elements/Gd.jl")
include("Elements/Tb.jl")
include("Elements/Dy.jl")
include("Elements/Ho.jl")
include("Elements/Er.jl")
include("Elements/Yb.jl")
include("Elements/Lu.jl")
include("Elements/Th.jl")
include("Elements/U.jl")

"""
    chemical_diffusion_list(search::String="")

List all available chemical diffusion data for melt.
Includes an argument to search for a specific term, i.e. an element ("La") or an author.
"""
function chemical_diffusion_list(search::String="")
    m = @__MODULE__
    all = string.(names(m; all=true, imported=true))
    Melt = Symbol.(filter(x -> startswith(x, "Melt") && contains(x, search), all))
    return [getfield(m, fun) for fun in Melt]
end


end
