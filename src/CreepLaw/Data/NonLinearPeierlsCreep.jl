module NonLinearPeierls

using GeoParams
export nonlinear_peierls_database, nonlinear_peierls_database_info

function dry_olivine_Mei_2010()
    data = NonLinearPeierlsCreep(;
        Name = "Dry Olivine | Mei et al. (2010)",
        n = 2.0NoUnits,                         # power-law exponent
        q = 1.0NoUnits,                         # outer exponent q in (1 - (σ/σP)^p)^q
        o = 0.5NoUnits,                         # inner stress exponent p
        TauP = 5.9e9Pa,                         # Peierls stress
        # A_P of eq. (3) as published, for axial strain rate and differential stress. The former extra factor
        # 2^((1+3/n)/2) belongs to the stress-scale convention B = 2^((1+3/n)/2) A^(-1/n) of other codes;
        # the conversion to invariants is done here by Apparatus = AxialCompression.
        A = 1.4e-7MPa^(-2) * s^(-1),            # material specific rheological parameter
        E = 320.0kJ / mol,                        # activation energy
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Mei et al. (2010), eq. (3) and paragraph [23]: A_P = 1.4e-7 MPa^-2/s, σP = 5.9 GPa, E_k(0) = 320 kJ/mol, p = 1/2, q = 1 (axial strain rate, differential stress). Reproduces the paper's data: log10 misfit -0.09 ± 0.18 (N = 11, Table 1).",
        BibTex_Reference = "
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

@inbounds nonlinear_peierls_database(f::F) where {F} = first(f())
@inbounds nonlinear_peierls_database_info(f::F) where {F} = last(f())

end
