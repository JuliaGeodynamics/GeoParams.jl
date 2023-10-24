# Add a list of pre-defined non linear peierls creep law values
export NonLinearPeierlsCreep_data

"""
    SetNonLinearPeierlsCreep["Name of non linear peierls creep law"]
This is a dictionary with pre-defined creep laws    
"""
function SetNonLinearPeierlsCreep(name::String)
    return Transform_NonLinearPeierlsCreep(name)
end

function SetNonLinearPeierlsCreep(
    name::String,
    CharDim::GeoUnits{T}
) where {T<:Union{GEO,SI}}
    return nondimensionalize(Transform_NonLinearPeierlsCreep(name), CharDim)
end

function NonLinearPeierlsCreep_data(name::String)
    if name === "Dry Olivine | Mei et al. (2010)"
        return NonLinearPeierlsCreep(;
            Name = "Dry Olivine | Mei et al. (2010)",
            n = 2.0NoUnits,                         # power-law exponent
            q = 1.0NoUnits,                         # exponent of water-fugacity
            o = 0.5NoUnits,                         # grain size exponent
            TauP = 5.9e9Pa,                         # Peierls stress
            A = (1.4e-7) * (2.0^((1.0 + (3.0 / 2.0)) / 2.0))MPa^(-2) / s,    # material specific rheological parameter
            E = 320.0kJ / mol,                      # activation energy
            Apparatus = AxialCompression,
        )
    elseif name === "Test Peierls"
        return NonLinearPeierlsCreep(;
            Name = "Test Peierls",
            n = 2.0NoUnits,                           # power-law exponent
            q = 1.0NoUnits,                           # exponent of water-fugacity
            o = 0.5NoUnits,                           # grain size exponent
            TauP = 5.9e9Pa,                           # Peierls stress
            A = 4.55657940893437e-9MPa^(-2) / s, # material specific rheological parameter
            E = 320.0kJ / mol,                        # activation energy
            Apparatus = AxialCompression,
        )
    end
    return NonLinearPeierlsCreep()
end
