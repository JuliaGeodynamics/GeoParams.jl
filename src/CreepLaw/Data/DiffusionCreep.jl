# Add a list of pre-defined diffusion creep law values
export DiffusionCreep_data

"""
    SetDiffusionCreep["Name of Diffusion Creep"]
This is a dictionary with pre-defined creep laws    
"""
@inline SetDiffusionCreep(name::String) = Transform_DiffusionCreep(name)

function SetDiffusionCreep(name::String, CharDim::GeoUnits{T}) where {T<:Union{GEO,SI}}
    return nondimensionalize(Transform_DiffusionCreep(name), CharDim)
end

function DiffusionCreep_data(name::String)
    if name === "Test Diff"
        return DiffusionCreep(;
            Name = "Test Diff",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = 2.070729911135297e-7/ MPa * m^3 / s,    # material specific rheological parameter
            E = 375.0kJ / mol,                      # activation energy
            V = 6.0e-6m^3 / mol,                    # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Anorthite | Rybacki et al. (2006)"
        return DiffusionCreep(;
            Name = "Dry Anorthite | Rybacki et al. (2006)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^12.1)/ MPa * μm^3 / s,  # material specific rheological parameter
            E = 460.0kJ / mol,                      # activation energy
            V = 24e-6m^3 / mol,                     # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Olivine | Mei & Kohlstedt (2000a)"
        return DiffusionCreep(;
            Name = "Wet Olivine | Mei & Kohlstedt (2000a)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 1.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^4.7)/ MPa * μm^3 / s, # material specific rheological parameter
            E = 295.0kJ / mol,                      # activation energy
            V = 20e-6m^3 / mol,                     # activation Volume
            Apparatus = AxialCompression,
        )

    elseif name === "Wet Olivine | Hirth & Kohlstedt (2003)"
        return DiffusionCreep(;
            Name = "Wet Olivine | Hirth & Kohlstedt (2003)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 1.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^7.4)/ MPa * μm^3 / s,    # material specific rheological parameter
            E = 375.0kJ / mol,                        # activation energy
            V = 20e-6m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Olivine | Hirth & Kohlstedt (2003)"
        return DiffusionCreep(;
            Name = "Dry Olivine | Hirth & Kohlstedt (2003)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^9.2)/ MPa * μm^3 / s,    # material specific rheological parameter
            E = 375.0kJ / mol,                        # activation energy
            V = 10e-6m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Olivine | Faul & Jackson (2006)"
        return DiffusionCreep(;
            Name = "Dry Olivine | Faul & Jackson (2006)",
            n = 1.4NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^10.3)MPa^(-7 // 5) * μm^3 / s,    # material specific rheological parameter
            E = 484.0kJ / mol,                        # activation energy
            V = 0.0m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Clinopyroxene | Hier-Majumder et al. (2005)"
        return DiffusionCreep(;
            Name = "Dry Clinopyroxene | Hier-Majumder et al. (2005)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = ((10^25.3) / 64.9e3)/ MPa * μm^3 / s,    # material specific rheological parameter
            E = 760.0kJ / mol,                        # activation energy
            V = 0.0m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Clinopyroxene | Hier-Majumder et al. (2005)"
        return DiffusionCreep(;
            Name = "Wet Clinopyroxene | Hier-Majumder et al. (2005)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 1.4NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^7.9 / 64.9e3)/ MPa * μm^3 / s,    # material specific rheological parameter
            E = 340.0kJ / mol,                        # activation energy
            V = 14e-6m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Clinopyroxene | Bystricky & Mackwell (2001)"
        return DiffusionCreep(;
            Name = "Dry Clinopyroxene | Bystricky & Mackwell (2001)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^15.1)/ MPa * μm^3 / s,    # material specific rheological parameter
            E = 560.0kJ / mol,                        # activation energy
            V = 0.0m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Diopside | Dimanov & Dresen (2005)"
        return DiffusionCreep(;
            Name = "Dry Diopside | Dimanov & Dresen (2005)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^14)/ MPa * μm^3 / s,    # material specific rheological parameter
            E = 528.0kJ / mol,                        # activation energy
            V = 0.0m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Diopside | Dimanov & Dresen (2005)"
        return DiffusionCreep(;
            Name = "Wet Diopside | Dimanov & Dresen (2005)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^8.1)/ MPa * μm^3 / s,    # material specific rheological parameter
            E = 528.0kJ / mol,                        # activation energy
            V = 0.0m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Clinopyroxene | Chen et al. (2006)"
        return DiffusionCreep(;
            Name = "Wet Clinopyroxene | Chen et al. (2006)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^8.1)/ MPa * μm^3 / s,    # material specific rheological parameter
            E = 528.0kJ / mol,                        # activation energy
            V = 0.0m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Anorthite | Rybacki & Dresen (2000)"
        return DiffusionCreep(;
            Name = "Dry Anorthite | Rybacki & Dresen (2000)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^12.1)/ MPa * μm^3 / s,    # material specific rheological parameter
            E = 467.0kJ / mol,                        # activation energy
            V = 0.0m^3 / mol,                       # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Anorthite | Rybacki & Dresen (2000)"
        return DiffusionCreep(;
            Name = "Wet Anorthite | Rybacki & Dresen (2000)",
            n = 1.0NoUnits,                # power-law exponent
            r = 0.0NoUnits,                # exponent of water-fugacity
            p = -3.0NoUnits,               # grain size exponent
            A = (10^1.7)/ MPa * μm^3 / s,  # material specific rheological parameter
            E = 170.0kJ / mol,             # activation energy
            V = 0.0m^3 / mol,              # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Anorthite | Rybacki et al. (2006)"
        return DiffusionCreep(;
            Name = "Wet Anorthite | Rybacki et al. (2006)",
            n = 1.0NoUnits,                # power-law exponent
            r = 1.0NoUnits,                # exponent of water-fugacity
            p = -3.0NoUnits,               # grain size exponent
            A = (10^-0.7)/ MPa * μm^3 / s, # material specific rheological parameter
            E = 159.0kJ / mol,             # activation energy
            V = 38.0m^3 / mol,             # activation Volume
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Quartzite | Rutter & Brodie (2004)"
        return DiffusionCreep(;
            Name = "Wet Quartzite | Rutter & Brodie (2004)",
            n = 1.0NoUnits,                # power-law exponent
            r = 0.0NoUnits,                # exponent of water-fugacity
            p = -2.0NoUnits,               # grain size exponent
            A = 0.39810717055349726/ MPa * μm^2 / s, # material specific rheological parameter
            # A = (10^-0.4)/ MPa * μm^3 / s, # material specific rheological parameter
            E = 220.0kJ / mol,             # activation energy
            V = 0.0m^3 / mol,              # activation Volume
            Apparatus = AxialCompression,
        )
    end

    return DiffusionCreep()
end
