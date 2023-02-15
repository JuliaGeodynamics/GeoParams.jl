# Add a list of pre-defined grain boundary sliding (GBS) values
export GrainBoundarySliding_info

"""
    SetGrainBoundarySliding["Name of GBS"]
This is a dictionary with pre-defined creep laws    
"""
SetGrainBoundarySliding(name::String; kwargs...) = Transform_GrainBoundarySliding(name; kwargs)

# predefined grain boundary sliding laws are to be added in the dictionary as it is done for dislocation creep laws (see 'DislocationCreep.jl')!
const GrainBoundarySliding_info = Dict([])