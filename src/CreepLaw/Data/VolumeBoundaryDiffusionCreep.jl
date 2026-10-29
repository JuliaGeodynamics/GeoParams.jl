module VolumeBoundaryDiffusion

using GeoParams
export volume_boundary_diffusion_database, volume_boundary_diffusion_database_info

function ice_Goldsby_2001()
    data = VolumeBoundaryDiffusionCreep(;
        Name = "Ice | Goldsby and Kohlstedt (2001)",
        D0v = 9.1e-4m^2 / s,      # pre-exponential factor of volume diffusion
        Qv = 59.4kJ / mol,        # activation energy of volume diffusion
        D0b = 6.0e-4m^2 / s,      # pre-exponential factor of grain-boundary diffusion
        Qb = 49.0kJ / mol,        # activation energy of grain-boundary diffusion
        δ = 9.04e-10m,            # grain-boundary width
        Vm = 1.97e-5m^3 / mol,    # molar volume
        Apparatus = AxialCompression,
    )
    info = MaterialParamsInfo(;
        Comment = "Goldsby & Kohlstedt (2001), eq. (4) and Table 6 (p. 11,025): Vm, D0v, Qv (Ramseier, 1967b), δ = 2b = 9.04e-10 m, Qb = 49 kJ/mol. D0b is not in Table 6; section 5.5 bounds it by 8.4e-4 (208 K) and 5.8e-4 m^2/s (193 K), and section 5.6 (p. 11,026) uses D0b = 6e-4 m^2/s, giving Db(273 K) = 2.53e-13 m^2/s (paper: ~2.5e-13). σ is the differential stress and ε̇ the axial strain rate of uniaxial compression experiments. Eq. (4) is not fitted to diffusion-creep data: GK2001 observed no diffusional flow. Eq. (4) with these values gives 1.8e-8 s^-1 (log10 misfit -0.71) for the section 5.5 bound of 9e-8 s^-1 (208 K, 2 MPa, d = 4 μm, Dv = 0, D0b = 8.4e-4 m^2/s), and at d = 1 mm, 233 K it exceeds the superplastic rate of Table 5 below about 0.01 MPa, whereas section 5.9 predicts no diffusional flow there.",
        BibTex_Reference = "
            @article{goldsby2001superplastic,
            title={Superplastic deformation of ice: Experimental observations},
            author={Goldsby, D. L. and Kohlstedt, D. L.},
            journal={Journal of Geophysical Research: Solid Earth},
            volume={106},
            number={B6},
            pages={11017--11030},
            year={2001},
            doi={10.1029/2000JB900336}
            }
    ",
    )
    return data, info
end

volume_boundary_diffusion_database(f::F) where {F} = first(f())
volume_boundary_diffusion_database_info(f::F) where {F} = last(f())

end
