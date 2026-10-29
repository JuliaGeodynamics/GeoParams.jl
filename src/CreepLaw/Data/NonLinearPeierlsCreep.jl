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

function antigorite_Burdette_2022()
    data = NonLinearPeierlsCreep(;
        Name = "Antigorite | Burdette & Hirth (2022)",
        n = 2.0NoUnits,                                     # power-law exponent of the (σ/μ)² term
        q = 1.18NoUnits,                                    # outer exponent q
        o = 1.0NoUnits,                                     # inner stress exponent p
        TauP = 2.42e9Pa,                                    # Peierls stress
        A = exp(-0.624) / (35.0e3)^2 * MPa^(-2) * s^(-1),   # A/μ² with μ = 35 GPa
        E = 86.3kJ / mol,                                   # activation energy
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Burdette & Hirth (2022), Eq. 1 and Table 2: ε̇ = A (σ/μ)² exp(-(ΔF/RT)(1-(σ/τ0)^p)^q), ln A = -0.624 (A in s^-1), μ = 35 GPa, ΔF = 86.3 kJ/mol, τ0 = 2.42 GPa, p = 1 (imposed), q = 1.18; σ is differential stress. Intact antigorite cores, 1 GPa, 75-550 °C, Griggs apparatus. Reproduces the paper's data: log10 misfit -0.01 ± 0.10 (N = 57, digitized from Fig. 5b; tabulated data are in Supporting Information S1), and the Sect. 4.5 viscosities at 50 MPa shear stress within 0.1 dex.",
        BibTex_Reference = "
            @article{burdette2022creep,
            title={Creep Rheology of Antigorite: Experiments at Subduction Zone Conditions},
            author={Burdette, Eric and Hirth, Greg},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={127},
            number={7},
            year={2022},
            doi={10.1029/2022JB024260}
            }
    ",
    )
    return data, info
end

function glaucophane_glide_Hufford_2026()
    data = NonLinearPeierlsCreep(;
        Name = "Glaucophane glide | Hufford et al. (2026)",
        n = 2.0NoUnits,                                     # power-law exponent of the (σ/μ)² term
        q = 2.0NoUnits,                                     # outer exponent q
        o = 1.0NoUnits,                                     # inner stress exponent p
        TauP = 22.0e9Pa,                                    # Peierls stress
        # A/μ² with μ = 60.6 GPa; × √3 because the paper uses ε̇ = γ̇/√3 and SimpleShear uses ε̇ = γ̇.
        A = sqrt(3) * 10^38.0181 / (60.6e3)^2 * MPa^(-2) * s^(-1),
        E = 800.0kJ / mol,                                  # activation enthalpy
        Apparatus = SimpleShear,
    )
    info = MaterialParamsInfo(;
        Comment = "Hufford et al. (2026), glide law (p. 11): ε̇ = A (σ/μ)² exp(-(H/RT)(1-σ/σP)^2), A = 10^38.0181 s^-1, μ = 60.6 GPa, σP = 22 GPa, H = 800 kJ/mol, p = 1, q = 2, for σ_equiv = 2τ and ε̇_equiv = γ̇/√3; stored as SimpleShear with A/μ² × √3. Glaucophane aggregates, 1 GPa, 600-700 °C, Griggs general shear. Reproduces the paper's data: log10 misfit -0.02 ± 0.21 (N = 12 high-stress steps, Table 1).",
        BibTex_Reference = "
            @article{hufford2026blueschist,
            title={Blueschist Dislocation Creep and Glide in Subduction Systems: Constraints From Glaucophane Experiments},
            author={Hufford, Lonnie J. and Tokle, Leif and Behr, Whitney M. and Morales, Luiz F. G. and Madonna, Claudio},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={131},
            number={3},
            year={2026},
            doi={10.1029/2025JB033622}
            }
    ",
    )
    return data, info
end

@inbounds nonlinear_peierls_database(f::F) where {F} = first(f())
@inbounds nonlinear_peierls_database_info(f::F) where {F} = last(f())

end
