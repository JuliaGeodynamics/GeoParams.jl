# Add a list of pre-defined peierls creep law values
export PeierlsCreep_info

"""
    SetPeierlsCreep["Name of peierls creep law"]
This is a dictionary with pre-defined creep laws    
"""
SetPeierlsCreep(name::String; kwargs...) = Transform_PeierlsCreep(name; kwargs)

# predefined peierls creep laws are to be added in the dictionary as it is done for dislocation creep laws (see 'DislocationCreep.jl')!
const PeierlsCreep_info = Dict([

# Dry Olivine rheology 
(
    "Dry Olivine | Goetze and Evans (1979)",
    (
        PeierlsCreep(;
            Name="Dry Olivine | Goetze and Evans (1979)",
            n=1.0NoUnits,                         # power-law exponent
            o=2.0NoUnits,                         # exponent of water-fugacity
            q=1.0NoUnits,                        # grain size exponent
            TauP=8.5e9Pa,                         # Peierls stress
            A=(5.7e11)s^(-1),    # material specific rheological parameter
            E=536.0kJ / mol,                        # activation energy
            Apparatus=AxialCompression,
        ),
        MaterialParamsInfo(;
            Comment="Checked values; not yet plots (NM)",
            BibTex_Reference="
                @article{goetze1979stress,
                title={Stress and temperature in the bending lithosphere as constrained by experimental rock mechanics},
                author={Goetze, Christopher and Evans, Brian},
                journal={Geophysical Journal International},
                volume={59},
                number={3},
                pages={463--478},
                year={1979},
                publisher={Blackwell Publishing Ltd Oxford, UK}
                }
        "),
    ),
)

# Dry Olivine rheology 
(
    "Dry Olivine | Demouchy (2013)",
    (
        PeierlsCreep(;
            Name="Dry Olivine | Demouchy (2013)",
            n=1.0NoUnits,                         # power-law exponent
            o=2.0NoUnits,                         # exponent of water-fugacity
            q=0.5NoUnits,                        # grain size exponent
            TauP=15.0e9Pa,                         # Peierls stress
            A=(1.0e6)s^(-1),    # material specific rheological parameter
            E=450.0kJ / mol,                        # activation energy
            Apparatus=AxialCompression,
        ),
        MaterialParamsInfo(;
            Comment="Paper not found; plots not checked (NM)",
            BibTex_Reference="
                @article{demouchy2013low,
                title={Low strength of Earth’s uppermost mantle inferred from tri-axial deformation experiments on dry olivine crystals},
                author={Demouchy, Sylvie and Tommasi, Andr{\'e}a and Ballaran, Tiziana Boffa and Cordier, Patrick},
                journal={Physics of the Earth and Planetary Interiors},
                volume={220},
                pages={37--49},
                year={2013},
                publisher={Elsevier}
                }
        "),
    ),
)

# Dry Olivine rheology 
(
    "Dry Olivine | Idrissei (2016)",
    (
        PeierlsCreep(;
            Name="Dry Olivine | Idrissei (2016)",
            n=1.0NoUnits,                         # power-law exponent
            o=2.0NoUnits,                         # exponent of water-fugacity
            q=0.5NoUnits,                        # grain size exponent
            TauP=3.8e9Pa,                         # Peierls stress
            A=(1e6)s^(-1),    # material specific rheological parameter
            E=566.0kJ / mol,                        # activation energy
            Apparatus=AxialCompression,
        ),
        MaterialParamsInfo(;
            Comment="Paper not found; plots not checked (NM)",
            BibTex_Reference=""),
    ),
)

])