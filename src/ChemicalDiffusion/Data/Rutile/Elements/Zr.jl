# ------------------------------------- Zr -------------------------------------

"""
    Rt_Zr_Cherniak2007_para_c

Diffusion data of Zr in rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 750-1100°C.
From Cherniak et al. (2007) (https://doi.org/10.1016/j.epsl.2007.06.027).
"""
function Rt_Zr_Cherniak2007_para_c()
    data = DiffusionData(
        Name = "Zr diffusion in Rutile (Ξc and anhydrous conditions) | Cherniak et al. (2007)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "Zr",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Buffer = "QFM, NNO, Air",
        Fluid = "Anhydrous",
        D0 = 9.8e-15m^2 / s,
        log_D0_1σ = 1.25NoUnits,
        Ea = 170kJ / mol,
        Ea_1σ = 30kJ / mol,
        Charge = 4,  # charge of the cation
        T_range_min = 750C,
        T_range_max = 1100C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (20.01.25)",
        BibTex_Reference = "
          @article{cherniak2007zr,
          title={Zr and Hf diffusion in rutile},
          author={Cherniak, DJ and Manchester, J and Watson, EB},
          journal={Earth and Planetary Science Letters},
          volume={261},
          number={1-2},
          pages={267--279},
          year={2007},
          publisher={Elsevier}
          }
          ",
    )

    return data, info
end

# """
#     Rt_Zr_Dohmen2019_para_c

# Diffusion data of Zr in rutile. Parallel c-axis. Calibrated between 800-1100°C.
# From Dohmen et al. (2019) (10.1007/s00269-018-1005-7).
# """
# function Rt_Zr_Dohmen2019_para_c()
#     data = DiffusionData(
#         Name = "Zr diffusion in Rutile (parallel c and anhydrous conditions) | Dohmen et al. (2019)",
#         Phase = "Rutile",
#         Formula = "TiO2",
#         Species = "Zr",
#         Orientation = "Ξc",
#         Crystallography = "Tetragonal",
#         Buffer = "QFM, fO2=-7",
#         Fluid = "Anhydrous",
#         D0 = (10^(-0.253))m^2 / s,
#         log_D0_1σ = (0.019*log(10))NoUnits,
#         Ea = (414/log(10))kJ / mol,
#         Ea_1σ = (11/log(10))kJ / mol,
#         Charge = 4,  # charge of the cation
#         bfO2 = 1e-7NoUnits,
#         T_range_min = 800C,
#         T_range_max = 1100C
#     )
#     info = MaterialParamsInfo(;
#         Comment = "Checked values by HD (20.01.25)",
#         BibTex_Reference = "
#           @article{Dohmen2019,
#             abstract = {We performed experiments with thin film diffusion couples to simultaneously measure diffusion coefficients of Zr, Hf, Nb and Ta parallel to the a- and c-axes of synthetic rutile in a gas mixing furnace at controlled oxygen fugacity at temperatures between 800 and 1100∘C. Depth profiles of the diffusion couples were measured using secondary-ion mass spectrometry. Some of the diffusion profiles show a concentration dependence, which indicates different diffusion mechanisms above and below a particular trace-element concentration level (∼1000μg/g). The diffusion coefficients for the mechanism dominant at high-concentration levels are approximately two orders of magnitude smaller than for the low-concentration mechanism. Below the critical concentration the diffusion coefficient is constant, as consistently shown in all of the experiments. For this diffusion coefficient we have found that D Zr ∼ D Nb > D Hf > > D Ta , and diffusion is isotropic for the four elements at all investigated T and fO 2 conditions. At 1000∘C for log fO 2 < FMQ+1, the diffusion coefficients decrease with increasing oxygen fugacity where D is proportional to fO2n with exponents n≈ - 0.25 for Zr and Hf and n≈ - 0.30 for Nb and Ta. Diffusivites of Nb and Ta strongly differ from each other at all investigated conditions, thus providing the potential to fractionate these geochemical twins, as suggested earlier. The present data and literature data for Zr and Ti self diffusion are interpreted and predicted based on published quantitative point defect models. Two end-member diffusion mechanisms were identified for impurity diffusion of Zr: (i) an interstitialcy mechanism involving Ti 3 + on interstitial sites, which is dominant at approximately log fO 2 < FMQ+2; (ii) a vacancy mechanism involving Ti vacancies, which is dominant at approximately log fO 2 > FMQ+2. The point defect calculations also explain the observed effects of heterovalent substitutions, such as Nb 5 + for Ti 4 + at high concentration levels for changes in the diffusion mechanism and hence diffusion rates. In the case of rutile, this concentration effect becomes much more sensitive to the substitution level at lower temperature. In natural rutile penta- and hexavalent cations may largely be charge balanced by mono-, di- and trivalent cations, such that the doping effect on diffusion may be reduced or may even be reversed. The Arrhenius relationships established here may therefore not be directly applicable to natural rutile. We obtained the following Arrhenius relationships (with diffusion coefficients D in m 2 / s , fO 2 in Pascal and T in Kelvin), which are only applicable for log fO 2 < FMQ+2: logDZr=(-0.40±0.47)+(-0.253±0.019)logfO210-7-414±11kJ/molRTln10logDHf=(-0.08±0.63)+(-0.266±0.023)logfO210-7-428±15kJ/molRTln10logDNb=(-0.19±0.36)+(-0.294±0.014)logfO210-7-421±9kJ/molRTln10logDTa=(0.45±0.73)+(-0.304±0.015)logfO210-7-463±18kJ/molRTln10.},
#             author = {Ralf Dohmen and Horst R. Marschall and Thomas Ludwig and Joana Polednia},
#             doi = {10.1007/s00269-018-1005-7},
#             isbn = {0123456789},
#             issn = {14322021},
#             issue = {3},
#             journal = {Physics and Chemistry of Minerals},
#             keywords = {Diffusion,Experiments,HFSE,Impurities,TiO 2},
#             pages = {311-332},
#             publisher = {Springer Berlin Heidelberg},
#             title = {Diffusion of Zr, Hf, Nb and Ta in rutile: effects of temperature, oxygen fugacity, and doping level, and relation to rutile point defect chemistry},
#             volume = {46},
#             url = {http://dx.doi.org/10.1007/s00269-018-1005-7},
#             year = {2019},
#             }

#           ",
#     )

#     return data, info
# end
