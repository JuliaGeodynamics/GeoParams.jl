"""
    Rt_Al_Cherniak2019_para_c

Diffusion data of Al in rutile. With anhydrous conditions and parallel to c-axis. 
Calibrated between 1100-1400°C. From Cherniak et al. (2019) (https://doi.org/10.2138/am-2019-7030).
"""
function Rt_Al_Cherniak2019_para_c()
    data = DiffusionData(
        Name = "Al diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Cherniak et al. (2019)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Al",  # element or species being diffused
        Orientation = "Ξc",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "NNO",  # Buffer condition (e.g., NNO) during the experiment
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 1.21e-2m^2 / s,  # pre-exponential factor
        log_D0_1σ = 2.1NoUnits,
        Ea = 531kJ / mol,  # activation energy
        Ea_1σ = 27kJ / mol,  # uncertainty at 1σ of the activation energy
        Charge = 3,  # charge of the cation
        T_range_min = 1100C,  # temperature min of the experiment
        T_range_max = 1400C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (15.01.25)",
        BibTex_Reference = "
          @article{Cherniak2019,
            abstract = {Diffusion of Al and Si has been measured in synthetic and natural rutile under anhydrous conditions. Experiments used Al2O3 or Al2O3-TiO2 powder mixtures for Al diffusant sources, and SiO2-TiO2 powder mixtures or quartz-rutile diffusion couples for Si. Experiments were run in air in crimped Pt capsules, or in sealed silica glass ampoules with solid buffers (to buffer at NNO or IW). Al profiles were measured with Nuclear Reaction Analysis (NRA) using the reaction 27Al(p,g)28Si. Rutherford Backscattering spectrometry (RBS) was used to measure Si diffusion profiles, with RBS also used in measurements of Al to complement NRA profiles. We determine the following Arrhenius relations from these measurements: For Al diffusion parallel to c, for experiments buffered at NNO, over the temperature range 1100-1400 °C: DAl = 1.21 × 10-2 exp(-531 ± 27 kJ/mol-1/RT) m2s-1. For Si diffusion parallel to c, for both unbuffered and NNO-buffered experiments, over the temperature range 1100-1450 °C: DSi = 8.53 × 10-13 exp(-254 ± 31 kJ/mol-1/RT) m2s-1. Diffusion normal to (100) is similar to diffusion normal to (001) for both Al and Si, indicating little diffusional anisotropy for these elements. Diffusivities measured for synthetic and natural rutile are in good agreement, indicating that these diffusion parameters can be applied in evaluating diffusivities in rutile in natural systems Diffusivities of Al and Si for experiments buffered at IW are faster (by a half to three-quarters of a log unit) than those buffered at NNO. Si and Al are among the slowest-diffusing species in rutile measured thus far. Diffusivities of Al and Si are significantly slower than the diffusion of Pb and slower than the diffusion of tetravalent Zr and Hf and pentavalent Nb and Ta. These data indicate that Al compositional information will be strongly retained in rutile, providing evidence for the robustness of the recently developed Al in rutile thermobarometer. For example, at 900 °C, Al compositional information would be preserved over ~3 Gyr in the center of 250 mm radius rutile grains, but Zr compositional information would be preserved for only about 300 000 yr at this temperature. Al-in-rutile compositions will also be much better preserved during subsolidus thermal events subsequent to crystallization than those for Ti-in-quartz and Zr-in-titanite crystallization thermometers.},
            author = {Daniele J. Cherniak and E. Bruce Watson},
            doi = {10.2138/am-2019-7030},
            issn = {19453027},
            issue = {11},
            journal = {American Mineralogist},
            keywords = {Rutherford backscattering,Rutile,aluminum,diffusion,geobarometry,geothermometry,nuclear reaction analysis,silicon},
            pages = {1638-1649},
            title = {Al and Si diffusion in rutile},
            volume = {104},
            year = {2019},
            }
          ",
    )

    return data, info
end
