# Add a list of pre-defined diffusion creep law values
export DiffusionCreep_info

"""
    SetDiffusionCreep["Name of Diffusion Creep"]
This is a dictionary with pre-defined creep laws    
"""
SetDiffusionCreep(name::String) = Transform_DiffusionCreep(name)

# predefined diffusion creep laws are to be added in the dictionary as it is done for dislocation creep laws (see 'DislocationCreep.jl')!
const DiffusionCreep_info = Dict([

# Dry Plagioclase rheology 
(
    "Dry Anorthite | Rybacki et al. (2006)",
    (
        DiffusionCreep(;
            Name="Dry Anorthite | Rybacki et al. (2006)",
            r=0.0NoUnits,                         # exponent of water-fugacity
            p=-3.0NoUnits,                        # grain size exponent
            A=(10^12.1)MPa^(-1) * μm^3.0 * s^(-1),    # material specific rheological parameter
            E=460.0kJ / mol,                        # activation energy
            V=24e-6m^3 / mol,                       # activation Volume
            Apparatus=AxialCompression,
        ),
        MaterialParamsInfo(;
            Comment="Checked values; not yet plots (BK)",
            BibTex_Reference="
                @article{Rybacki_Gottschalk_Wirth_Dresen_2006, 
                title={Influence of water fugacity and activation volume on the flow properties of fine-grained anorthite aggregates}, 
                volume={111}, 
                DOI={10.1029/2005JB003663}, 
                number={B3}, 
                journal={Journal of Geophysical Research: Solid Earth}, 
                author={Rybacki, E. and Gottschalk, M. and Wirth, R. and Dresen, G.}, 
                year={2006}, 
                month={Mar}
                }
            ",
        ),
    ),
)])  # end of list
