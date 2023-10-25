# Add a list of pre-defined grain boundary sliding (GBS) values
export GrainBoundarySliding_data

"""
    SetGrainBoundarySliding["Name of GBS"]
This is a dictionary with pre-defined creep laws    
"""
function SetGrainBoundarySliding(name::String)
    return Transform_GrainBoundarySliding(name)
end

function SetGrainBoundarySliding(
    name::String,
    CharDim::GeoUnits{T}
) where {T<:Union{GEO,SI}}
    nondimensionalize(Transform_GrainBoundarySliding(name), CharDim)
end

function GrainBoundarySliding_data(name::String)
    if name === "Dry Olivine < 1523K | Hirth and Kohlstedt (2003)"
        return GrainBoundarySliding(;
            Name = "Dry Olivine < 1523K | Hirth and Kohlstedt (2003)",
            n = 3.5NoUnits,                           # power-law exponent
            p = -2.0NoUnits,                          # grain size exponent
            A = 6500.0MPa^(-7 // 2) * μm^(2) / s, # material specific rheological parameter
            E = 400.0kJ / mol,                        # activation energy
            V = 18.0e-6m^3 / mol,                     # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Test GBS"
        return GrainBoundarySliding(;
            Name = "Test GBS",
            n = 3.5NoUnits,                         # power-law exponent
            p = -2.0NoUnits,                        # grain size exponent
            A = 1.506190693026593e2MPa^(-7 // 2) * m^(2) / s, # material specific rheological parameter
            E = 600.0kJ / mol,                      # activation energy
            V = 18.0e-6m^3 / mol,                   # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Olivine >= 1523K | Hirth and Kohlstedt (2003)"
        return GrainBoundarySliding(;
            Name = "Dry Olivine >= 1523K | Hirth and Kohlstedt (2003)",
            n = 3.5NoUnits,                           # power-law exponent
            p = -2.0NoUnits,                          # grain size exponent
            A = 4.7e10MPa^(-7 // 2) * μm^(2) / s, # material specific rheological parameter
            E = 600.0kJ / mol,                        # activation energy
            V = 18e-6m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )

    end

    return GrainBoundarySliding()
end
