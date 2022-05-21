# Add a list of pre-defined diffusion creep law values
"""
    SetDiffusionCreep["Name of Diffusion Creep"]
This is a dictionary with pre-defined creep laws    
"""
SetDiffusionCreep(name::String) = DiffusionCreep_info[name][1]

# predefined diffusion creep laws are to be added in the dictionary as it is done for dislocation creep laws (see 'DislocationCreep.jl')!
DiffusionCreep_info = Dict([

    # Dry Plagioclase rheology 
    ("Dry Anorthite | Bürgmann & Dresen (2008)", 
        (DiffusionCreep(
            Name = "Dry Anorthite | Bürgmann & Dresen (2008)",
            n = 1.0NoUnits,                         # power-law exponent
            r = 0.0NoUnits,                         # exponent of water-fugacity
            p = -3.0NoUnits,                        # grain size exponent
            A = (10^12.1)MPa^(-1)*μm^3.0*s^(-1),    # material specific rheological parameter
            E = 460.0kJ/mol,                        # activation energy
            V = 24e-6m^3/mol,                       # activation Volume
            Apparatus = AxialCompression),
            
            MaterialParamsInfo(Comment = "Taken from Bürgmann & Dresen (2008), supplementary table 1 .",
            
            BibTex_Reference = parse_bibtex("""
                @article{Bürgmann_Dresen_2008,
                address = {Washington, D. C.},
                title = {Rheology of the Lower Crust and Upper Mantle: Evidence from Rock Mechanics, Geodesy, and Field Observations},
                volume = {36},
                ISSN={0084-6597, 1545-4495},
                journal={Annual Review of Earth and Planetary Sciences}, 
                author={Bürgmann, Roland and Dresen, Georg}, 
                year={2008}, 
                month={May}, 
                pages={531–567},
                }
            """))
        )
    )




])  # end of list
