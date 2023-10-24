# Add a list of pre-defined peierls creep law values
export PeierlsCreep_data

"""
    SetPeierlsCreep["Name of peierls creep law"]
This is a dictionary with pre-defined creep laws    
"""
SetPeierlsCreep(name::String) = Transform_PeierlsCreep(name)

function SetPeierlsCreep(name::String, CharDim::GeoUnits{GEO})
    return nondimensionalize(Transform_PeierlsCreep(name), CharDim)
end

function PeierlsCreep_data(name::String)
    if name === "Dry Olivine | Goetze and Evans (1979)"
        return PeierlsCreep(;
            Name = "Dry Olivine | Goetze and Evans (1979)",
            n = 1.0NoUnits,     # power-law exponent
            o = 2.0NoUnits,     # exponent of water-fugacity
            q = 1.0NoUnits,     # grain size exponent
            TauP = 8.5e9Pa,     # Peierls stress
            A = (5.7e11) / s, # material specific rheological parameter
            E = 536.0kJ / mol,  # activation energy
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Olivine | Demouchy (2013)"
        return PeierlsCreep(;
            Name = "Dry Olivine | Demouchy (2013)",
            n = 1.0NoUnits,    # power-law exponent
            o = 2.0NoUnits,    # exponent of water-fugacity
            q = 0.5NoUnits,    # grain size exponent
            TauP = 15.0e9Pa,   # Peierls stress
            A = (1.0e6) / s, # material specific rheological parameter
            E = 450.0kJ / mol, # activation energy
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Olivine | Idrissei (2016)"
        return PeierlsCreep(;
            Name = "Dry Olivine | Idrissei (2016)",
            n = 1.0NoUnits,    # power-law exponent
            o = 2.0NoUnits,    # exponent of water-fugacity
            q = 0.5NoUnits,    # grain size exponent
            TauP = 3.8e9Pa,    # Peierls stress
            A = (1e6) / s,   # material specific rheological parameter
            E = 566.0kJ / mol, # activation energy
            Apparatus = AxialCompression,
        )
    end
    return PeierlsCreep()
end
