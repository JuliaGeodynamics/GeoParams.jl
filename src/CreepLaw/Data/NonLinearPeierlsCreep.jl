# Add a list of pre-defined non linear peierls creep law values
export NonLinearPeierlsCreep_info

"""
    SetNonLinearPeierlsCreep["Name of non linear peierls creep law"]
This is a dictionary with pre-defined creep laws    
"""
SetNonLinearPeierlsCreep(name::String; kwargs...) = Transform_NonLinearPeierlsCreep(name; kwargs)

# predefined peierls creep laws are to be added in the dictionary as it is done for dislocation creep laws (see 'DislocationCreep.jl')!
const NonLinearPeierlsCreep_info = Dict([

# Wet Olivine rheology of Mei 2010
(
    "Dry Olivine | Mei et al. (2010)",
    (
        NonLinearPeierlsCreep(;
            Name="Dry Olivine | Mei et al. (2010)",
            n=2.0NoUnits,                         # power-law exponent
            q=1.0NoUnits,                         # exponent of water-fugacity
            o=0.5NoUnits,                        # grain size exponent
            TauP=5.9e9Pa,                         # Peierls stress
            A=(1.4e-7)*(2.0^((1.0+(3.0/2.0))/2.0))MPa^(-2.0) * s^(-1.0),    # material specific rheological parameter
            E=320.0kJ / mol,                        # activation energy
            Apparatus=AxialCompression,
        ),
        MaterialParamsInfo(;
            Comment="Checked values; not yet plots (NM)",
            BibTex_Reference="
                @article{mei2010experimental,
                title={Experimental constraints on the strength of the lithospheric mantle},
                author={Mei, S and Suzuki, AM and Kohlstedt, DL and Dixon, NA and Durham, WB},
                journal={Journal of Geophysical Research: Solid Earth},
                volume={115},
                number={B8},
                year={2010},
                publisher={Wiley Online Library}
                }
        "),
    ),
)

# Wet Olivine rheology of Mei 2010 TEST
(
    "TEST PEIERLS",
    (
        NonLinearPeierlsCreep(;
            Name="TEST PEIERLS",
            n=2.0NoUnits,                         # power-law exponent
            q=1.0NoUnits,                         # exponent of water-fugacity
            o=0.5NoUnits,                        # grain size exponent
            TauP=5.9e9Pa,                         # Peierls stress
            A=4.55657940893437e-9MPa^(-2.0) * s^(-1.0),    # material specific rheological parameter
            E=320.0kJ / mol,                        # activation energy
            Apparatus=AxialCompression,
        ),
        MaterialParamsInfo(;
            Comment="Law to check Thorsten Beckers book rheology and its conversion of the A factor",
            BibTex_Reference=""),
    ),
)

# Wet Olivine rheology of Mei 2010 TEST
(
    "THORSTEN PEIERLS",
    (
        NonLinearPeierlsCreep(;
            Name="THORSTEN PEIERLS",
            n=2.0NoUnits,                         # power-law exponent
            q=1.0NoUnits,                         # exponent of water-fugacity
            o=0.5NoUnits,                        # grain size exponent
            TauP=5.9e9Pa,                         # Peierls stress
            A=2^((1.0+(3.0/3.5))/2.0) * 1.4e-7MPa^(-2.0) * s^(-1.0),    # material specific rheological parameter
            E=320.0kJ / mol,                        # activation energy
            Apparatus=AxialCompression,
        ),
        MaterialParamsInfo(;
            Comment="Law to check Thorsten Beckers book rheology and its conversion of the A factor",
            BibTex_Reference=""),
    ),
)
])