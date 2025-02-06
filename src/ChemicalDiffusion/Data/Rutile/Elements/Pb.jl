"""
    Rt_Pb_Cherniak2000_unor

Diffusion data of Pb in rutile. With anhydrous conditions and unoriented. 
Calibrated between 700-1100°C. From Cherniak et al. (2000) 
(https://doi.org/10.1007/PL00007671).
"""
function Rt_Pb_Cherniak2000_unor()
    data = DiffusionData(
        Name = "Pb diffusion in Rutile (unoriented and anhydrous conditions) | Cherniak et al. (2000)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Pb",  # element or species being diffused
        Orientation = "unoriented",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "QFM, NNO, self-buffered",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 3.9e-10m^2 / s,  # pre-exponential factor
        log_D0_1σ = 1.04NoUnits,
        Ea = 250kJ / mol,  # activation energy
        Ea_1σ = 12kJ / mol,  # uncertainty at 1σ of the activation energy
        Charge = 2,  # charge of the cation
        T_range_min = 700C,  # temperature min of the experiment
        T_range_max = 1100C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (15.01.25)",
        BibTex_Reference = "
          @article{Cherniak2000,
            abstract = {Diffusion of Pb was measured in natural and synthetic rutile under dry, 1 atmosphere conditions, using mixtures of Pb titanate or Pb sulfide and TiO2 as the sources of diffusant. Pb depth profiles were then measured with Rutherford Backscattering Spectrometry (RBS). Over the temperature range 700-1100 °C, the following Arrhenius relation was obtained for the synthetic rutile: D = 3.9 x 10-10exp(-250 ± 12 kJ mol-1/RT) m2s-1. Results for diffusion in natural and synthetic rutile were quite similar, despite significant differences in trace element compositions. Mean closure temperatures calculated from the diffusion parameters are around 600 °C for rutile grains of ~100 μm size. This is about 100 °C higher than rutile closure temperature determinations from past field-based studies, suggesting that rutile is more resistant to Pb loss through volume diffusion than previously thought.},
            author = {D. J. Cherniak},
            doi = {10.1007/PL00007671},
            issn = {00107999},
            issue = {2},
            journal = {Contributions to Mineralogy and Petrology},
            pages = {198-207},
            title = {Pb diffusion in rutile},
            volume = {139},
            year = {2000},
            }
          ",
    )

    return data, info
end
