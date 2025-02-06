
"""
    Rt_Ti_Hoshino1985_para_c

Diffusion data of Ti self-diffusion in rutile. Parallel c-axis, anhydrous. 
Calibrated between 1000-1500C.From Hoshino et al. (1985) 
(https://doi.org/10.1016/0022-3697(85)90079-4).
"""
function Rt_Ti_Hoshino1985_para_c()
    data = DiffusionData(
        Name = "Ti self-diffusion in Rutile (parallel to c-axis) | Hoshino et al. (1985)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "Ti",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Fluid = "anhydrous",
        D0 = 6.5cm^2 / s,
        log_D0_1σ = log(1.33e-4)NoUnits,
        Ea = (66.11)u"kcal / mol",
        Ea_1σ = (0.56)u"kcal / mol",
        Charge = -2,
        T_range_min = 1000C,
        T_range_max = 1500C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",        
        BibTex_Reference = "
          @article{Hoshino1985,
            abstract = {The self-diffusion of 44Ti has been measured both parallel to and perpendicular to the c axis in rutile single crystals by a serial-sectioning technique as a function of temperature (1000-1500°C) and oxygen partial pressure (10-14 - 1 atm). The oxygen-partial-pressure dependence of. D*Ti indicates that cation selfdiffusion occurs by an interstitial-type mechanism and that both trivalent and tetravalent interstitial titanium ions may contribute to cation self-diffusion. At po2 = 1.50 × 10-7 atm where impurity-induced defects are unimportant,D*Ti(∥c)=6.50 +1.33 -1.11exp- (66.11±0.56 kcal mole RT cm2 S and D*Ti(⊥c)= 4.55 +1.78 -1.28exp- (64.08±0.99) kcal mole RT cm2 S. In the intrinsic region, the ratio D*Ti (⊥c)/D*Ti(∥c) was found to increase from 1.2 to 1.6 as the temperature decreased from 1500 to 1000°C. Computations based upon the defect model of Kofstad (involving the atomic defects Ti...iTi....iand V..o), of Marucco et al. (Ti....i and V..o), and of Blumenthal et al. (Ti...i and Ti....i) are compared with the experimental data on deviation from stoichiometry, electrical conductivity, cation self-diffusion and chemical diffusion in TiO2-x. These comparisons provide values of the defect concentrations, cation-defect diffusivities, electron mobility and reasonable values of the correlation factor for cation diffusion by the interstitialcy mechanism. Only the model of Kofstad is inconsistent with the data. © 1985.},
            author = {K. Hoshino and N. L. Peterson and C. L. Wiley},
            doi = {10.1016/0022-3697(85)90079-4},
            issn = {00223697},
            issue = {12},
            journal = {Journal of Physics and Chemistry of Solids},
            keywords = {TiO2-x,chemical diffusion,diffusion mechanisms,electrical conductivity,point defects,self-diffusion (cation)},
            pages = {1397-1411},
            title = {Diffusion and point defects in TiO2-x},
            volume = {46},
            year = {1985},
            }
          ",
    )

    return data, info
end


"""
    Rt_Ti_Hoshino1985_perp_c

Diffusion data of Ti self-diffusion in rutile. Perpendicular c-axis, anhydrous. 
Calibrated between 1000-1500C.From Hoshino et al. (1985) 
(https://doi.org/10.1016/0022-3697(85)90079-4).
"""
function Rt_Ti_Hoshino1985_perp_c()
    data = DiffusionData(
        Name = "Ti self-diffusion in Rutile (Perpendicular to c-axis) | Hoshino et al. (1985)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "Ti",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "anhydrous",
        D0 = 4.55cm^2 / s,
        log_D0_1σ = log(1.78e-4)NoUnits,
        Ea = (64.08)u"kcal / mol",
        Ea_1σ = (0.99)u"kcal / mol",
        Charge = -2,
        T_range_min = 1000C,
        T_range_max = 1500C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",        
        BibTex_Reference = "
          @article{Hoshino1985,
            abstract = {The self-diffusion of 44Ti has been measured both parallel to and perpendicular to the c axis in rutile single crystals by a serial-sectioning technique as a function of temperature (1000-1500°C) and oxygen partial pressure (10-14 - 1 atm). The oxygen-partial-pressure dependence of. D*Ti indicates that cation selfdiffusion occurs by an interstitial-type mechanism and that both trivalent and tetravalent interstitial titanium ions may contribute to cation self-diffusion. At po2 = 1.50 × 10-7 atm where impurity-induced defects are unimportant,D*Ti(∥c)=6.50 +1.33 -1.11exp- (66.11±0.56 kcal mole RT cm2 S and D*Ti(⊥c)= 4.55 +1.78 -1.28exp- (64.08±0.99) kcal mole RT cm2 S. In the intrinsic region, the ratio D*Ti (⊥c)/D*Ti(∥c) was found to increase from 1.2 to 1.6 as the temperature decreased from 1500 to 1000°C. Computations based upon the defect model of Kofstad (involving the atomic defects Ti...iTi....iand V..o), of Marucco et al. (Ti....i and V..o), and of Blumenthal et al. (Ti...i and Ti....i) are compared with the experimental data on deviation from stoichiometry, electrical conductivity, cation self-diffusion and chemical diffusion in TiO2-x. These comparisons provide values of the defect concentrations, cation-defect diffusivities, electron mobility and reasonable values of the correlation factor for cation diffusion by the interstitialcy mechanism. Only the model of Kofstad is inconsistent with the data. © 1985.},
            author = {K. Hoshino and N. L. Peterson and C. L. Wiley},
            doi = {10.1016/0022-3697(85)90079-4},
            issn = {00223697},
            issue = {12},
            journal = {Journal of Physics and Chemistry of Solids},
            keywords = {TiO2-x,chemical diffusion,diffusion mechanisms,electrical conductivity,point defects,self-diffusion (cation)},
            pages = {1397-1411},
            title = {Diffusion and point defects in TiO2-x},
            volume = {46},
            year = {1985},
            }
          ",
    )

    return data, info
end