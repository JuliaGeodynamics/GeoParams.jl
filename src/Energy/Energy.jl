include("HeatCapacity.jl")
include("Conductivity.jl")
include("LatentHeat.jl")
include("RadioactiveHeat.jl")
include("Shearheating.jl")

# export reduce_energy_source_terms

# # Single phase methods

# @inline function reduce_energy_source_terms(f::NTuple{N, Any}, rheology::MaterialParams, specific_args) where N
#     sum(_reduce_energy_source_terms(fᵢ, rheology, specific_argsᵢ) for (fᵢ, specific_argsᵢ) in zip(f, specific_args))
# end

# @inline function reduce_energy_source_terms(f::NTuple{N, Any}, rheology::MaterialParams, specific_args, source_terms::Vararg{NH, T}) where {N, NH, T}
#     H = reduce_energy_source_terms(f, rheology, specific_args)
#     return sum(source_terms; init = H)
# end

# @inline function _reduce_energy_source_terms(f::F, rheology::MaterialParams, specific_args) where F
#     f(rheology, specific_args...)
# end

# @inline function _reduce_energy_source_terms(f::F, rheology::MaterialParams, specific_args::Union{NamedTuple, Tuple{}}) where F
#     f(rheology, specific_args)
# end
