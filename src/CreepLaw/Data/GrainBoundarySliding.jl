# Add a list of pre-defined grain boundary sliding (GBS) values
export GrainBoundarySliding_info

"""
    SetGrainBoundarySliding["Name of GBS"]
This is a dictionary with pre-defined creep laws    
"""
SetGrainBoundarySliding(name::String; kwargs...) = Transform_GrainBoundarySliding(name; kwargs)

# predefined grain boundary sliding laws are to be added in the dictionary as it is done for dislocation creep laws (see 'DislocationCreep.jl')!
const GrainBoundarySliding_info = Dict([
# Dry Olivine rheology, T < 1523 K
(
    "Dry Olivine < 1523K | Hirth and Kohlstedt (2003)",
    (
        DiffusionCreep(;
            Name="Dry Olivine < 1523K | Hirth and Kohlstedt (2003)",
            n=3.5NoUnits,                         # power-law exponent
            p=-2.0NoUnits,                        # grain size exponent
            A=6500MPa^(-3.5) * μm^3.0 * s^(-1),    # material specific rheological parameter
            E=400.0kJ / mol,                        # activation energy
            V=18e-6m^3 / mol,                       # activation Volume
            Apparatus=AxialCompression,
        ),
        MaterialParamsInfo(;
            Comment="Checked values; not yet plots (NM)",
            BibTex_Reference="
                @article{hirth2004rheology,
                title={Rheology of the Upper Mantle and the Mantle Wedge: A View from the Experimentalists},
                author={Hirth, Greg and Kohlstedt, David},
                journal={Inside the Subduction Factory},
                volume={138},
                pages={83--105},
                year={2004},
                publisher={Wiley Online Library}
                }
        "),
    ),
)
])