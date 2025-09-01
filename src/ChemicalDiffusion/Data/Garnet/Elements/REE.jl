"""
    Grt_REE_Bloch2020_fast()

Diffusion data of REE in garnet for fast diffusion mechanism. Calibrated with experiments conducted between 950-1050°C at 1 atm with the QFM buffer and dry condition plus data from the literature.
From Bloch et al. (2020) (https://doi.org/10.1093/petrology/egaa055) combined with data from Van Orman et al. (2002), Tirone et al. (2005) and Bloch et al. (2015).
"""
function Grt_REE_Bloch2020_fast()
    data = DiffusionData(
        Name = "REE diffusion in Garnet (QFM) | Bloch et al. (2020)",
        Phase = "Garnet",  # name of the mineral
        Formula = "X3Y2(SiO4)3",  # chemical formula of the mineral
        Species = "REE",  # element or species being diffused
        Orientation = "Isotropic",  # Crystal orientation from the diffusion experiment
        Crystallography = "Isometric",  # Crystallographic system of the mineral
        Buffer = "QFM",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = (10^(-10.24))u"m^2 / s",  # pre-exponential factor (log10(-10.24) in the original paper)
        log_D0_1σ = (0.21 * log(10))NoUnits,
        Ea = (221057)u"J/mol",  # activation energy
        Ea_1σ = (4284)u"J/mol",  # uncertainty at 1σ of the activation energy
        Charge = 3,  # charge of the cation
        T_range_min = 950C,  # temperature min of the experiment
        T_range_max = 1050C,  # temperature max of the experiment
        P0 = 1.0u"atm"  # pressure of calibration
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (31.08.25). Values originally reported in log10 but converted here to ln.",
        BibTex_Reference = "
          @article{bloch2020multispecies,
          title={Multispecies diffusion of yttrium, rare earth elements and hafnium in garnet},
          author={Bloch, EM and Jollands, MC and Devoir, A and Bouvier, A-S and Ibañez-Mejia, M and Baumgartner, LP},
          journal={Journal of Petrology},
          volume={61},
          number={7},
          pages={egaa055},
          year={2020},
          publisher={Oxford University Press}
          }
          ",
    )

    return data, info
end


"""
    Grt_REE_Bloch2020_slow()

Diffusion data of REE in garnet for slow diffusion mechanism. Calibrated with experiments conducted between 950-1050°C at 1 atm with the QFM buffer and dry condition plus data from the literature.
From Bloch et al. (2020) (https://doi.org/10.1093/petrology/egaa055) combined with data from Van Orman et al. (2002), Tirone et al. (2005) and Bloch et al. (2015).
"""
function Grt_REE_Bloch2020_slow()
    data = DiffusionData(
        Name = "REE diffusion in Garnet (QFM) | Bloch et al. (2020)",
        Phase = "Garnet",  # name of the mineral
        Formula = "X3Y2(SiO4)3",  # chemical formula of the mineral
        Species = "REE",  # element or species being diffused
        Orientation = "Isotropic",  # Crystal orientation from the diffusion experiment
        Crystallography = "Isometric",  # Crystallographic system of the mineral
        Buffer = "QFM",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = (10^(-9.28))m^2 / s,  # pre-exponential factor (log10(-9.28) in the original paper)
        log_D0_1σ = 0.65NoUnits,
        Ea = (265200)u"J/mol",  # activation energy (fitted in log10 in the original paper)
        Ea_1σ = (38540)u"J/mol",  # uncertainty at 1σ of the activation energy
        ΔV = (10800.0e-9)u"m^3 / mol",  # activation volume (with Pressure in GPa in the original paper)
        ΔV_1σ = (2600.0e-9)u"m^3 / mol",  # uncertainty at 1σ of the activation volume
        T_range_min = 950C,  # temperature min of the experiment
        T_range_max = 1050C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (31.08.25). Values originally reported in log10 but converted here to ln.",
        BibTex_Reference = "
          @article{bloch2020multispecies,
          title={Multispecies diffusion of yttrium, rare earth elements and hafnium in garnet},
          author={Bloch, EM and Jollands, MC and Devoir, A and Bouvier, A-S and Ibañez-Mejia, M and Baumgartner, LP},
          journal={Journal of Petrology},
          volume={61},
          number={7},
          pages={egaa055},
          year={2020},
          publisher={Oxford University Press}
          }
          ",
    )

    return data, info
end
