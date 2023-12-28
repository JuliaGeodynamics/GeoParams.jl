module NonLinearPeierls

using GeoParams
export nonlinear_peierls_database, nonlinear_peierls_database_info

function dry_olivine_Mei_2010()
    data = NonLinearPeierlsCreep(;
        Name="Dry Olivine | Mei et al. (2010)",
        n=2.0NoUnits,                         # power-law exponent
        q=1.0NoUnits,                         # exponent of water-fugacity
        o=0.5NoUnits,                        # grain size exponent
        TauP=5.9e9Pa,                         # Peierls stress
        A=(1.4e-7) * (2.0^((1.0 + (3.0 / 2.0)) / 2.0))MPa^(-2) * s^(-1),    # material specific rheological parameter
        E=320.0kJ / mol,                        # activation energy
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
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
    ",
    )
    return data, info
end

function test_Peierls()
    data = NonLinearPeierlsCreep(;
        Name="Test Peierls",
        n=2.0NoUnits,                         # power-law exponent
        q=1.0NoUnits,                         # exponent of water-fugacity
        o=0.5NoUnits,                        # grain size exponent
        TauP=5.9e9Pa,                         # Peierls stress
        A=4.55657940893437e-9MPa^(-2) * s^(-1),    # material specific rheological parameter
        E=320.0kJ / mol,                        # activation energy
        Apparatus=AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment="Law to check Thorsten Beckers book rheology and its conversion of the A factor",
        BibTex_Reference="",
    )
    return data, info
end

@inbounds nonlinear_peierls_database(f::F) where {F} = first(f())
@inbounds nonlinear_peierls_database_info(f::F) where {F} = last(f())

end