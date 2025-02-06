module Rutile

using GeoParams
using Unitful

include("Elements/H.jl")
include("Elements/He.jl")
include("Elements/Li.jl")
include("Elements/O.jl")
include("Elements/Al.jl")
include("Elements/Si.jl")
include("Elements/Ti.jl")
include("Elements/Zr.jl")
include("Elements/Nb.jl")
include("Elements/Hf.jl")
include("Elements/Ta.jl")
include("Elements/Pb.jl")

"""
    chemical_diffusion_list(search::String="")

List all available chemical diffusion data for melt.
Includes an argument to search for a specific term, i.e. an element ("La") or an author.
"""
function chemical_diffusion_list(search::String = "")
    m = @__MODULE__
    all = string.(names(m; all = true, imported = true))
    Rt = Symbol.(filter(x -> startswith(x, "Rt") && contains(x, search), all))
    return [getfield(m, fun) for fun in Rt]
end


end
