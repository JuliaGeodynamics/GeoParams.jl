export SetDislocationCreep, DislocationCreep_data
# This contains predefined dislocation creep values - Feel free to expand
"""
    SetDislocationCreep["Name of Dislocation Creep"]

Sets predefined dislocation creep data from a dictionary

"""
SetDislocationCreep(name::String) = Transform_DislocationCreep(name)

function SetDislocationCreep(name::String, CharDim::GeoUnits{T}) where {T<:Union{GEO,SI}}
    return nondimensionalize(Transform_DislocationCreep(name), CharDim)
end

function DislocationCreep_data(name::String)
    if name === "Dry Olivine | Hirth & Kohlstedt (2003)"
        return DislocationCreep(;
            Name = "Dry Olivine | Hirth & Kohlstedt (2003)",
            n = 3.5NoUnits,
            r = 0.0NoUnits,
            A = 1.1e5MPa^(-7 // 2) / s,
            E = 530.0kJ / mol,
            V = 14e-6m^3 / mol,
            Apparatus = AxialCompression,
        )
    elseif name === "Test Disl"
        return DislocationCreep(;
            Name = "Test Disl",
            n = 3.5NoUnits,
            r = 0.0NoUnits,
            A = 7.440446357925734e7MPa^(-7 // 2) / s,
            E = 530.0kJ / mol,
            V = 14e-6m^3 / mol,
            Apparatus = AxialCompression,
        )
    elseif name === "1. Wet Olivine | Hirth & Kohlstedt (2003)"
        return DislocationCreep(;
            Name = "1. Wet Olivine | Hirth & Kohlstedt (2003)",
            n = 3.5NoUnits,
            A = 90MPa^(-7 // 2) / s,
            E = 480kJ / mol,
            V = 11e-6m^3 / mol,
            r = 1.2NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "2. Wet Olivine | Hirth & Kohlstedt (2003)"
        return DislocationCreep(;
            Name = "2. Wet Olivine | Hirth & Kohlstedt (2003)",
            n = 3.0NoUnits,
            A = 1600MPa^(-3) / s,
            E = 520.0kJ / mol,
            V = 22.0m^3 / mol,
            r = 1.2NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Quartz Diorite | Hansen & Carter (1982)"
        return DislocationCreep(;
            Name = "Quartz Diorite | Hansen & Carter (1982)",
            n = 2.25NoUnits,
            A = 3.5e-2MPa^(-9 // 4) / s,
            E = 212kJ / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Diabase | Caristan (1982)"
        return DislocationCreep(;
            Name = "Diabase | Caristan (1982)",
            n = 3.05NoUnits,
            A = 6.0e-2MPa^(-61 // 20) / s,
            E = 276kJ / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Tumut Pond Serpentinite | Raleigh and Paterson (1965)"
        return DislocationCreep(;
            Name = "Tumut Pond Serpentinite | Raleigh and Paterson (1965)",
            n = 2.8NoUnits,
            A = 6.3e-7MPa^(-14 // 5) / s,
            E = 66kJ / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Maryland strong diabase | Mackwell et al. (1998)"
        return DislocationCreep(;
            Name = "Maryland strong diabase | Mackwell et al. (1998)",
            n = 4.7NoUnits,
            A = 8MPa^(-47 // 10) / s,
            E = 485kJ / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Mafic Granulite | Wilks and Carter (1990)"
        return DislocationCreep(;
            Name = "Mafic Granulite | Wilks and Carter (1990)",
            n = 4.2NoUnits,
            A = 1.4e4MPa^(-21 // 5) / s,
            E = 445kJ / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Quartzite | Ueda et al. (2008)"
        return DislocationCreep(;
            Name = "Wet Quartzite | Ueda et al. (2008)",
            n = 2.3NoUnits,
            A = 1 * exp10(-3.5)MPa^(-23 // 10) / s,
            E = 154kJ / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Granite | Carter and Tsenn (1987)"
        return DislocationCreep(;
            Name = "Granite | Carter and Tsenn (1987)",
            n = 3.3NoUnits,
            A = 1.0 * exp10(-5.7)MPa^(-33 // 10) / s,
            E = 186.5kJ / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Plagioclase An75 | Ji and Zhao (1993)"
        return DislocationCreep(;
            Name = "Plagioclase An75 | Ji and Zhao (1993)",
            n = 3.2NoUnits,
            A = 3.27e-4MPa^(-16 // 5) / s,
            E = 238kJ / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Anorthite | Rybacki et al. (2006)"
        return DislocationCreep(;
            Name = "Dry Anorthite | Rybacki et al. (2006)",
            n = 3.0NoUnits,
            A = exp10(12.7)MPa^(-3) / s,
            E = 641kJ / mol,
            V = 24e-6m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Anorthite | Rybacki et al. (2006)"
        return DislocationCreep(;
            Name = "Wet Anorthite | Rybacki et al. (2006)",
            n = 3.0NoUnits,
            A = (10^0.2)MPa^(-3) / s,
            E = 345kJ / mol,
            V = 38m^3 / mol,
            r = 1NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Olivine | Hirth and Kohlstedt (2003)"
        return DislocationCreep(;
            Name = "Wet Olivine | Hirth and Kohlstedt (2003)",
            n = 3.5NoUnits,
            A = 1600.0MPa^(-7 // 2) / s,
            E = 520kJ / mol,
            V = 11.0e-6m^3 / mol,
            r = 1.2NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Quarzite | Kirby (1983)"
        return DislocationCreep(;
            Name = "Wet Quarzite | Kirby (1983)",
            n = 2.3NoUnits,
            A = uconvert(Pa^(-23 // 10) / s, 3.2e-4MPa^(-23 // 10) / s),
            E = 154e3J / mol,
            V = 0m^3 / mol,
            r = 0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Upper Mantle Olivine | Afonso and Ranalli (2004)"
        return DislocationCreep(;
            Name = "Wet Upper Mantle Olivine | Afonso and Ranalli (2004)",
            n = 4.0NoUnits,
            A = 2.0e3MPa^(-4) / s,
            E = 471kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Granite | Tirel et al. (2008)"
        return DislocationCreep(;
            Name = "Granite | Tirel et al. (2008)",
            n = 3.2NoUnits,
            A = 1.25e-9MPa^(-16 // 5) / s,
            E = 123kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Olivine | Gerya (2019)"
        return DislocationCreep(;
            Name = "Dry Olivine | Gerya (2019)",
            n = 3.5NoUnits,
            A = uconvert(MPa^(-7 // 2) / s, 2.5e-17Pa^(-7 // 2) / s),
            E = 532.0kJ / mol,
            V = 0.0m^3 / mol,
            Apparatus = AxialCompression,   # used in book (according to matlab script)
            r = 0.0NoUnits,
        )
    elseif name === "Salado rock salt | Li & Urai (2016)"
        return DislocationCreep(;
            Name = "Salado rock salt | Li & Urai (2016)",
            n = 5.0NoUnits,
            A = 7.26e-6MPa^(-5) / s,
            E = 53.92kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Rock salt | Li & Urai (2016)"
        return DislocationCreep(;
            Name = "Rock salt | Li & Urai (2016)",
            n = 5.0NoUnits,
            A = 7.26e-6MPa^(-5) / s,
            E = 53.92kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Olivine | Mei & Kohlstedt (2000b)"
        return DislocationCreep(;
            Name = "Wet Olivine | Mei & Kohlstedt (2000b)",
            n = 3.0NoUnits,
            A = (10^3.2)MPa^(-3) / s,
            E = 470.0kJ / mol,
            V = 20.0m^3 / mol,
            r = 0.98NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Olivine | Karato & Jung (2003)"
        return DislocationCreep(;
            Name = "Dry Olivine | Karato & Jung (2003)",
            n = 3.0NoUnits,
            A = (10^6.1)MPa^(-3) / s,
            E = 510.0kJ / mol,
            V = 14.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Olivine | Karato & Jung (2003)"
        return DislocationCreep(;
            Name = "Wet Olivine | Karato & Jung (2003)",
            n = 3.0NoUnits,
            A = (10^2.9)MPa^(-3) / s,
            E = 510.0kJ / mol,
            V = 24.0m^3 / mol,
            r = 1.2NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Clinopyroxene | Chen et al. (2006)"
        return DislocationCreep(;
            Name = "Wet Clinopyroxene | Chen et al. (2006)",
            n = 2.7NoUnits,
            A = (10^6.7)MPa^(-27 // 10) / s,
            E = 670.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 3.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Clinopyroxene | Bystricky & Mackwell (2001)"
        return DislocationCreep(;
            Name = "Dry Clinopyroxene | Bystricky & Mackwell (2001)",
            n = 4.7NoUnits,
            A = (10^9.8)MPa^(-47 // 10) / s,
            E = 760.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Clinopyroxene | Bystricky & Mackwell (2001)"
        return DislocationCreep(;
            Name = "Dry Clinopyroxene | Bystricky & Mackwell (2001)",
            n = 4.7NoUnits,
            A = (10^10.8)MPa^(-47 // 10) / s,
            E = 760.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Diopside | Dimanov & Dresen (2005)"
        return DislocationCreep(;
            Name = "Dry Diopside | Dimanov & Dresen (2005)",
            n = 5.5NoUnits,
            A = uconvert(MPa^(-55 // 10) / s, 3.01e-28Pa^(-55 // 10) / s),
            E = 691.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Diopside | Dimanov & Dresen (2005)"
        return DislocationCreep(;
            Name = "Wet Diopside | Dimanov & Dresen (2005)",
            n = 5.5NoUnits,
            A = uconvert(MPa^(-55 // 10) / s, 5.16e-33Pa^(-55 // 10) / s),
            E = 534.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Omphacite | Zhang et al. (2006)"
        return DislocationCreep(;
            Name = "Wet Omphacite | Zhang et al. (2006)",
            n = 3.5NoUnits,
            A = (10^-2)MPa^(-7 // 2) / s,
            E = 310.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Jadeit | Orzol et al. (2006)"
        return DislocationCreep(;
            Name = "Wet Jadeit | Orzol et al. (2006)",
            n = 3.7NoUnits,
            A = (10^-3.3)MPa^(-37 // 10) / s,
            E = 326.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Anorthite | Rybacki & Dresen (2000)"
        return DislocationCreep(;
            Name = "Dry Anorthite | Rybacki & Dresen (2000)",
            n = 3.0NoUnits,
            A = (10^12.7)MPa^(-3) / s,
            E = 648.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Anorthite | Rybacki & Dresen (2000)"
        return DislocationCreep(;
            Name = "Wet Anorthite | Rybacki & Dresen (2000)",
            n = 3.0NoUnits,
            A = (10^0.2)MPa^(-3) / s,
            E = 356.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Quartzite | Rutter & Brodie (2004)"
        return DislocationCreep(;
            Name = "Wet Quartzite | Rutter & Brodie (2004)",
            n = 3.0NoUnits,
            A = (10^-4.9)MPa^(-3) / s,
            E = 242.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 1.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Quartzite | Hirth et al. (2001)"
        return DislocationCreep(;
            Name = "Wet Quartzite | Hirth et al. (2001)",
            n = 4.0NoUnits,
            A = (10^-11.2)MPa^(-4) / s,
            E = 135.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 1.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Dry Quartzite | Jaoul et al. (1984)"
        return DislocationCreep(;
            Name = "Dry Quartzite | Jaoul et al. (1984)",
            n = 2.8NoUnits,
            A = (10^-5.415)MPa^(-14 // 5) / s,
            E = 184.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Quartzite | Jaoul et al. (1984)"
        return DislocationCreep(;
            Name = "Wet Quartzite | Jaoul et al. (1984)",
            n = 2.8NoUnits,
            A = (10^-5.045)MPa^(-14 // 5) / s,
            E = 163.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 0.0NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Quartzite | Tokle et al. (2019)"
        return DislocationCreep(;
            Name = "Wet Quartzite | Tokle et al. (2019)",
            n = 3.0NoUnits,
            A = (10^-11.959)MPa^(-3) / s,
            E = 115.0kJ / mol,
            V = 0.0m^3 / mol,
            r = 1.2NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Quartzite | Lu and Jiang (2019)"
        return DislocationCreep(;
            Name = "Wet Quartzite | Lu and Jiang (2019)",
            n = 3.0NoUnits,
            A = (10^-14.2218)MPa^(-3) / s,
            E = 132.0kJ / mol,
            V = 35.3e-6m^3 / mol,
            r = 2.7NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "low pressure wet Quartzite | Lusk et al. (2021)"
        return DislocationCreep(;
            Name = "low pressure wet Quartzite | Lusk et al. (2021)",
            n = 3.5NoUnits,
            A = (10^-9.3)MPa^(-7 // 2) / s,
            E = 118.0kJ / mol,
            V = 2.59e-6m^3 / mol,
            r = 0.49NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "high pressure wet Quartzite | Lusk et al. (2021)"
        return DislocationCreep(;
            Name = "high pressure wet Quartzite | Lusk et al. (2021)",
            n = 2.1NoUnits,
            A = (10^-6.36)MPa^(-21 // 10) / s,
            E = 94.0kJ / mol,
            V = 1.44e-6m^3 / mol,
            r = 0.2NoUnits,
            Apparatus = AxialCompression,
        )
    elseif name === "Wet Quartzite | Lusk et al. (2021)"
        return DislocationCreep(;
            Name = "Wet Quartzite | Lusk et al. (2021)",
            n = 2.0NoUnits,
            A = (10^-7.9)MPa^(-2) / s,
            E = 77.0kJ / mol,
            V = 2.59e-6m^3 / mol,
            r = 0.49NoUnits,
            Apparatus = AxialCompression,
        )
    end

    return DislocationCreep()
end
