# Package extension for running GeoParams laws inside Metal kernels.
module GeoParamsMetalExt

using GeoParams, Metal

# Metal has no Float64, so on the device `retry_wider` does not widen Float32:
# a Float32 result that overflows or underflows is returned as Inf or 0.
# Nondimensionalized laws stay within the Float32 range.
Base.Experimental.@overlay Metal.method_table GeoParams._wider(::Type{Float32}) = Float32

function __init__()
    ccall(:jl_generating_output, Cint, ()) == 1 && return
    @info """GeoParams on Metal: Metal has no Float64, and SI creep-law parameters such as 1e-55 lie outside the Float32 range.
    Nondimensionalize the parameters, then convert them once on the host before passing them to a kernel:
        phases = convert_precision(Float32, nondimensionalize(phases, GEO_units()))"""
    return nothing
end

end
