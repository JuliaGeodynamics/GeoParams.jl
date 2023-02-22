# Add a list of pre-defined peierls creep law values
export PeierlsCreep_info

"""
    SetPeierlsCreep["Name of peierls creep law"]
This is a dictionary with pre-defined creep laws    
"""
SetPeierlsCreep(name::String; kwargs...) = Transform_PeierlsCreep(name; kwargs)

# predefined peierls creep laws are to be added in the dictionary as it is done for dislocation creep laws (see 'DislocationCreep.jl')!
const GrainBoundarySliding_info = Dict([


])