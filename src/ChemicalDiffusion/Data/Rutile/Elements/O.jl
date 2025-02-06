"""
    Rt_O_Arita1979_para_c

Diffusion data of O in rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 1150-1450K.
From Arita et al. (1979) (https://doi.org/10.1111/j.1151-2916.1979.tb19101.x).
"""
function Rt_O_Arita1979_para_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Arita et al. (1979)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        D0 = 3.4e-7m^2 / s,
        Ea = 251kJ / mol,
        Ea_1σ = 75.3kJ / mol,
        Charge = -2,
        T_range_min = 1150K,
        T_range_max = 1450K,
        P0 = 6000Pa
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Arita1979,
            abstract = {Secondary ion mass spectrumètry has been applied for measuring the tracer diffusivity of oxygen in the c direction of single‐crystal rutile for a temperature range of 1150 to 1450 K at 6000 Pa pressure of oxygen gas. Specimens diffusion‐annealed in oxygen gas containing 18O were subsequently continuously sputtered and analyzed for 16O and 18O. The tracer diffusivity was determined from the depth profile of 18O, taking into account a surface exchange reaction of oxygen. The tracer diffusivity in Cr2O3‐doped rutile was 3 to 8 times larger than that in pure rutile. For pure rutile, the diffusivity is expressed by D(m2/s)=3.4×10−7, exp [‐251(kJ/mol)/RT], and for 0.08 mol% Cr2O3‐doped rutile, by D(m2/s)= 2.0×10−8 exp[‐204(KJ/mol)/RT]. The Cr2O3 doping had a catalytic effect on the rate constant of the surface exchange reaction on the c surface. The rate constant is represented, for pure rutile, by K(m/s)= 2.4×10−1 exp[‐246(KJ/mol)/RT], and for 0.08 mol% Cr2O3‐doped rutile, k(m/s)= 3.5×10−5 exp[‐131(KJ/mol)/RT]. Copyright © 1979, Wiley Blackwell. All rights reserved},
            author = {M. Arita and M. Hosoya and M. Kobayashi and M. Someno},
            doi = {10.1111/j.1151-2916.1979.tb19101.x},
            issn = {15512916},
            issue = {9-10},
            journal = {Journal of the American Ceramic Society},
            pages = {443-446},
            title = {Depth Profile Measurement by Secondary Ion Mass Spectrometry for Determining the Tracer Diffusivity of Oxygen in Rutile},
            volume = {62},
            year = {1979},
        }

          ",
    )

    return data, info
end

"""
    Rt_O_Arita1979_para_c_Cr

Diffusion data of O in Cr-doped rutile. With anhydrous conditions and parallel to c-axis. Calibrated between 1150-1450K.
From Arita et al. (1979) (https://doi.org/10.1111/j.1151-2916.1979.tb19101.x).
"""
function Rt_O_Arita1979_para_c_Cr()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Arita et al. (1979)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Doping = "Cr",
        D0 = 2.0e-8m^2 / s,
        Ea = 204kJ / mol,
        Ea_1σ = 61.2kJ / mol,
        Charge = -2,
        T_range_min = 1150K,
        T_range_max = 1450K,
        P0 = 6000Pa
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Arita1979,
            abstract = {Secondary ion mass spectrumètry has been applied for measuring the tracer diffusivity of oxygen in the c direction of single‐crystal rutile for a temperature range of 1150 to 1450 K at 6000 Pa pressure of oxygen gas. Specimens diffusion‐annealed in oxygen gas containing 18O were subsequently continuously sputtered and analyzed for 16O and 18O. The tracer diffusivity was determined from the depth profile of 18O, taking into account a surface exchange reaction of oxygen. The tracer diffusivity in Cr2O3‐doped rutile was 3 to 8 times larger than that in pure rutile. For pure rutile, the diffusivity is expressed by D(m2/s)=3.4×10−7, exp [‐251(kJ/mol)/RT], and for 0.08 mol% Cr2O3‐doped rutile, by D(m2/s)= 2.0×10−8 exp[‐204(KJ/mol)/RT]. The Cr2O3 doping had a catalytic effect on the rate constant of the surface exchange reaction on the c surface. The rate constant is represented, for pure rutile, by K(m/s)= 2.4×10−1 exp[‐246(KJ/mol)/RT], and for 0.08 mol% Cr2O3‐doped rutile, k(m/s)= 3.5×10−5 exp[‐131(KJ/mol)/RT]. Copyright © 1979, Wiley Blackwell. All rights reserved},
            author = {M. Arita and M. Hosoya and M. Kobayashi and M. Someno},
            doi = {10.1111/j.1151-2916.1979.tb19101.x},
            issn = {15512916},
            issue = {9-10},
            journal = {Journal of the American Ceramic Society},
            pages = {443-446},
            title = {Depth Profile Measurement by Secondary Ion Mass Spectrometry for Determining the Tracer Diffusivity of Oxygen in Rutile},
            volume = {62},
            year = {1979},
        }

          ",
    )

    return data, info
end


"""
    Rt_O_Dennis1993_perp_c

Diffusion data of O in rutile. With hydrothermal conditions and perpendicular to c-axis. 
Calibrated between 873-1373K. From Dennis et al. (1993) (https://doi.org/10.1007/BF00414275).
"""
function Rt_O_Dennis1993_perp_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (perpendicular to c-axis and hydrothermal conditions) | Dennis et al. (1993)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "hydrothermal",
        D0 = 2.41e-12m^2 / s,
        log_D0_1σ = 25.4NoUnits,
        Ea = 172.5kJ / mol,
        Ea_1σ = 23.6kJ / mol,
        Charge = -2,
        T_range_min = 873K,
        T_range_max = 1373K,
        P0 = 100MPa
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Dennis1993,
            abstract = {Oxygen self-diffusion coefficients have been determined for synthetic and natural rutile single crystals under hydrothermal conditions at 100 MPa total pressure and in the temperature range 873-1373 K. The diffusion coefficients are lower than the results from dry gas exchange studies would predict. Between 973 and 1373 K the results can be characterized by two linear Arrhenius relationships. D=1.14×10-11 exp(-168.8 kJ mol-1/RT) m-2s-1 for the natural rutile, and D=2.41×10-12 exp(-172.5 kJ mol-1/RT) m2s-1 for the synthetic crystal. The results have been interpreted in terms of a defect model involving the dissolution of water in rutile as substitutional hydroxyl defects on oxygen lattice sites, with a solution enthalpy in the range 81-106 kJmol-1. © 1993 Chapman & Hall.},
            author = {P. F. Dennis and R. Freer},
            doi = {10.1007/BF00414275},
            issn = {00222461},
            issue = {17},
            journal = {Journal of Materials Science},
            pages = {4804-4810},
            title = {Oxygen self-diffusion in rutile under hydrothermal conditions},
            volume = {28},
            year = {1993},
        }
          ",
    )

    return data, info
end


"""
    Rt_O_Dennis1993_perp_c_nat

Diffusion data of O in natural rutile. With hydrothermal conditions and perpendicular to c-axis. 
Calibrated between 873-1373K. From Dennis et al. (1993) (https://doi.org/10.1007/BF00414275).
"""
function Rt_O_Dennis1993_perp_c_nat()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (perpendicular to c-axis and hydrothermal conditions) | Dennis et al. (1993)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        Fluid = "hydrothermal",
        Doping = "natural",
        D0 = 1.14e-11m^2 / s,
        log_D0_1σ = 23.1NoUnits,
        Ea = 168.8kJ / mol,
        Ea_1σ = 32.9kJ / mol,
        Charge = -2,
        T_range_min = 873K,
        T_range_max = 1373K,
        P0 = 100MPa
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Dennis1993,
            abstract = {Oxygen self-diffusion coefficients have been determined for synthetic and natural rutile single crystals under hydrothermal conditions at 100 MPa total pressure and in the temperature range 873-1373 K. The diffusion coefficients are lower than the results from dry gas exchange studies would predict. Between 973 and 1373 K the results can be characterized by two linear Arrhenius relationships. D=1.14×10-11 exp(-168.8 kJ mol-1/RT) m-2s-1 for the natural rutile, and D=2.41×10-12 exp(-172.5 kJ mol-1/RT) m2s-1 for the synthetic crystal. The results have been interpreted in terms of a defect model involving the dissolution of water in rutile as substitutional hydroxyl defects on oxygen lattice sites, with a solution enthalpy in the range 81-106 kJmol-1. © 1993 Chapman & Hall.},
            author = {P. F. Dennis and R. Freer},
            doi = {10.1007/BF00414275},
            issn = {00222461},
            issue = {17},
            journal = {Journal of Materials Science},
            pages = {4804-4810},
            title = {Oxygen self-diffusion in rutile under hydrothermal conditions},
            volume = {28},
            year = {1993},
        }
          ",
    )

    return data, info
end


"""
    Rt_O_Derry1981_para_c

Diffusion data of O in natural rutile. With parallel to c-axis. 
Calibrated between 1173-1673K. From Derry et al. (1981) (https://doi.org/10.1016/0022-3697(81)90011-1).
"""
function Rt_O_Derry1981_para_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis) | Derry et al. (1981)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        D0 = 2.4e-6m^2 / s,
        Ea = 282.6kJ / mol,
        Ea_1σ = 5kJ / mol,
        Charge = -2,
        T_range_min = 1173K,
        T_range_max = 1673K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Derry1981,
            author = {D. J. Derry and D. G. Lees and J. M. Calvert},
            doi = {https://doi.org/10.1016/0022-3697(81)90011-1},
            journal = {J. Phys. Chem. Solids},
            pages = {57-64},
            title = {A study of oxygen self-diffusion in the c-direction of rutile using a nuclear technique},
            volume = {42},
            year = {1981},
            }
          ",
    )

    return data, info
end

"""
    Rt_O_Haul1965_unor

Diffusion data of O in rutile. Unoriented, anhydrous. 
Calibrated between 710-1300C. From Haul et al. (1965) (https://doi.org/10.1016/0022-3697(65)90066-1).
"""
function Rt_O_Haul1965_unor()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (Unoriented) | Haul et al. (1965)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Unoriented",
        Crystallography = "Tetragonal",
        Fluid = "anhydrous",
        D0 = 2.0e-3cm^2 / s,
        Ea = (60.0e3)u"J / mol",
        Ea_1σ = (1.5e3)u"J / mol",
        Charge = -2,
        T_range_min = 710C,
        T_range_max = 1300C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
          @article{Haul1965,
            abstract = {Zusammenfassung-Die Sauerstoff-Selbstdiffusion in Rutilkristallen wurde durch gas/fest Isoto-penaustausch mit 180 markiertem Sauerstoff im Temperaturbereich von 710-1300°C untersucht. D* = 2,0.10-s exp[-(60 _L l,S)lOs/Rn. Aus den Diffusionsmessungen ergibt sich fiir die Fehlord-nung in den untersuchten Rutilkristallen das Vorliegen von Sauerstoff-Leerstellen. Eine Abhtigig-keit vom Sauerstoffdruck (760-10-s Torr) wurde nicht beobachtet, da die Konzentration der durch den Fremdoxidgehalt (AlaOs) bedingten Leerstellen diejenige der thermisch erzeugten iibersteigt. Die Sauerstoff-Diffusion ist senkrecht zur c-Achse griisser als parallel c und wird durch Anwesen-heit von Ha0 erniedrigt. Es werden die Bedingungen diskutiert, unter denen sich bei der Auswertung der Diffusions-versuche die aerlagerung einer Phasengrenzreaktion bemerkbar macht. Abstract-Oxygen diffusion in rutile crystals has been measured in the temperature range 710-13OO'C by means of gas/solid isotope exchange with IsO labelled oxygen. D* = 2,0-10-s exp[-(60 f 1,s) x 10S/RT]. With respect to disorder in the studied rutile crystals the results are clear evidence of the presence of oxygen vacancies. No dependence on oxygen pressure (760-10-s Torr) has been observed since the concentration of vacancies due to the impurity oxide content (AlsOs) is in excess to the concentration of those formed thermally. Oxygen diffusion is greater _L than 11 to the c-axis and is decreased by the presence of HsO. Conditions are discussed under which the superposition of a phase boundary exchange reaction becomes noticeable in the evaluation of diffusion experiments.},
            author = {R Haul and G. Dümbgen},
            doi = {https://doi.org/10.1016/0022-3697(65)90066-1},
            journal = {J. Phys. Chem. Solids},
            pages = {1-10},
            title = {Sauerstoff-Selbstdiffusion in Rutilkirstallen},
            volume = {26},
            url = {https://ac.els-cdn.com/0022369765900661/1-s2.0-0022369765900661-main.pdf?_tid=c8da796d-6617-4d05-b688-98d52f35e0dd&acdnat=1552934512_56530708d8184c76bdf6bd271d85bf76},
            year = {1965},
            }
          ",
    )

    return data, info
end


"""
    Rt_O_Lundy1973_para_c

Diffusion data of O in natural rutile. With parallel to c-axis. 
Calibrated between 1200-1500K. From Lundy et al. (1973) (https://doi.org/10.1051/jphyscol:1973953).
"""
function Rt_O_Lundy1973_para_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis) | Lundy et al. (1973)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        D0 = 0.046cm^2 / s,
        Ea = 59.9u"cal/mol",
        Charge = -2,
        T_range_min = 1200K,
        T_range_max = 1500K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
        @article{Lundy1973,
        author = {T Lundy and W Coghlan and T Lundy and W Coghlan Processus and D E Transport Dans and L E S Oxydescation},
        doi = {10.1051/jphyscol:1973953},
        journal = {Journal de Physique Colloques},
        pages = {C9-299-C9-302},
        title = {PROCESSUS DE TRANSPORT DANS LES OXYDESCATION SELF DIFFUSION IN RUTILE},
        volume = {34},
        year = {1973},
        }
          ",
    )

    return data, info
end


"""
    Rt_O_Lundy1973_perp_c

Diffusion data of O in natural rutile. With perpendicular to c-axis. 
Calibrated between 1200-1500K. From Lundy et al. (1973) (https://doi.org/10.1051/jphyscol:1973953).
"""
function Rt_O_Lundy1973_perp_c()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (perpendicular to c-axis) | Lundy et al. (1973)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "⊥c",
        Crystallography = "Tetragonal",
        D0 = 0.0024cm^2 / s,
        Ea = 48.5u"cal/mol",
        Charge = -2,
        T_range_min = 1200K,
        T_range_max = 1500K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
        @article{Lundy1973,
        author = {T Lundy and W Coghlan and T Lundy and W Coghlan Processus and D E Transport Dans and L E S Oxydescation},
        doi = {10.1051/jphyscol:1973953},
        journal = {Journal de Physique Colloques},
        pages = {C9-299-C9-302},
        title = {PROCESSUS DE TRANSPORT DANS LES OXYDESCATION SELF DIFFUSION IN RUTILE},
        volume = {34},
        year = {1973},
        }
          ",
    )

    return data, info
end


"""
    Rt_O_Moore1998_para_c_fast

Diffusion data of O in natural rutile. With hydrous and anhydrous conditions, parallel to c-axis, at 0.1 to 1000 MPa. 
Calibrated between 1200-1500K. From Moore et al. (1998) (https://doi.org/10.2138/am-1998-7-803).
"""
function Rt_O_Moore1998_para_c_fast()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis) | Moore et al. (1998)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Buffer = "1atm to NNO",
        D0 = 4.7e-7m^2 / s,
        Ea = 258.0e3J / mol,
        Ea_1σ = 22.0e3J / mol,
        Charge = -2,
        T_range_min = 1200K,
        T_range_max = 1500K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{Moore1998,
            abstract = {Oxygen self-diffusion in rutile was studied in synthetic and natural samples over the temperature range 750 to 1000 °C and the pressure range 0.1 to 1000 MPa using an 180-enriched source. Most experiments investigated the dependence of D on temperature, water pressure, crystallographic direction, and experiment duration. A few experiments investigated the dependence of D on fo2 and confining pressure. The uptake profiles of 18O in experimental products were measured by nuclear reaction analysis using the reaction 18O(p,α)15N. Two mechanisms are responsible for O diffusion in rutile, and one is faster than the other by about an order of magnitude. O that diffuses by the faster mechanism is described by the diffusion law: D(∥c) = 4.7 × 10-7 exp(-258 ±22 × 103/RT) D0 in mVs; £EA in J/mol; Tin K. Diffusion by the slower mechanism is described by this law: D(∥c) = 5.9 × 10-5 exp(-330 ±15 × 103VRT) D0 in mVs; EA in J/mol; T in K. Oxygen fugacity in itself does not affect D at fugacities between 1 atm and Ni-NiO. However, the presence or absence of water during reduction does affect the diffusion behavior. When water is absent during rutile growth and/or subsequent reduction, only the faster mechanism operates, and when water is present during growth or reduction, both mechanisms operate simultaneously, though the contribution from the slow mechanism dominates that of the fast mechanism. Because few geologic environments are truly dry, the slower law should generally be used for modeling O diffusion for rutile in nature. Comparison with other studies of rutile suggests that migration of O vacancies is the mechanism responsible for the faster diffusion law whereas migration of Ti interstitials is responsible for the slower diffusion law. Oxygen diffusion in rutile, is slower perpendicular to the c axis than parallel to that axis by about half an order of magnitude. There is no perceptible effect of confining pressure on D below 100 MPa, or between 600 and 1000 MPa. However, between 100 and 600 MPa, D decreases by nearly an order of magnitude. Closure temperatures for O diffusion in rutile are high - 650 °C for a crystal with a 100 Um radius and a 10 °C/Ma cooling rate. Rutile is retentive of its O isotopic composition. A crystal with a 100 μ radius will retain its initial core composition for just over 10 million years at 600 °C.},
            author = {D. K. Moore and Daniele J. Cherniak and E. B. Watson},
            doi = {10.2138/am-1998-7-803},
            issn = {0003004X},
            journal = {American Mineralogist},
            pages = {700-711},
            title = {Oxygen diffusion in rutile from 750 to 1000 °C and 0.1 to 1000 MPa},
            volume = {83},
            year = {1998},
            }
          ",
    )

    return data, info
end


"""
    Rt_O_Moore1998_para_c_slow

Diffusion data of O in natural rutile. With hydrous and anhydrous conditions, parallel to c-axis, at 0.1 to 1000 MPa. 
Calibrated between 1200-1500K. From Moore et al. (1998) (https://doi.org/10.2138/am-1998-7-803).
"""
function Rt_O_Moore1998_para_c_slow()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile (parallel to c-axis) | Moore et al. (1998)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Orientation = "Ξc",
        Crystallography = "Tetragonal",
        Buffer = "1atm to NNO",
        D0 = 5.9e-5m^2 / s,
        Ea = 330.0e3J / mol,
        Ea_1σ = 15.0e3J / mol,
        Charge = -2,
        T_range_min = 1200K,
        T_range_max = 1500K
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{Moore1998,
            abstract = {Oxygen self-diffusion in rutile was studied in synthetic and natural samples over the temperature range 750 to 1000 °C and the pressure range 0.1 to 1000 MPa using an 180-enriched source. Most experiments investigated the dependence of D on temperature, water pressure, crystallographic direction, and experiment duration. A few experiments investigated the dependence of D on fo2 and confining pressure. The uptake profiles of 18O in experimental products were measured by nuclear reaction analysis using the reaction 18O(p,α)15N. Two mechanisms are responsible for O diffusion in rutile, and one is faster than the other by about an order of magnitude. O that diffuses by the faster mechanism is described by the diffusion law: D(∥c) = 4.7 × 10-7 exp(-258 ±22 × 103/RT) D0 in mVs; £EA in J/mol; Tin K. Diffusion by the slower mechanism is described by this law: D(∥c) = 5.9 × 10-5 exp(-330 ±15 × 103VRT) D0 in mVs; EA in J/mol; T in K. Oxygen fugacity in itself does not affect D at fugacities between 1 atm and Ni-NiO. However, the presence or absence of water during reduction does affect the diffusion behavior. When water is absent during rutile growth and/or subsequent reduction, only the faster mechanism operates, and when water is present during growth or reduction, both mechanisms operate simultaneously, though the contribution from the slow mechanism dominates that of the fast mechanism. Because few geologic environments are truly dry, the slower law should generally be used for modeling O diffusion for rutile in nature. Comparison with other studies of rutile suggests that migration of O vacancies is the mechanism responsible for the faster diffusion law whereas migration of Ti interstitials is responsible for the slower diffusion law. Oxygen diffusion in rutile, is slower perpendicular to the c axis than parallel to that axis by about half an order of magnitude. There is no perceptible effect of confining pressure on D below 100 MPa, or between 600 and 1000 MPa. However, between 100 and 600 MPa, D decreases by nearly an order of magnitude. Closure temperatures for O diffusion in rutile are high - 650 °C for a crystal with a 100 Um radius and a 10 °C/Ma cooling rate. Rutile is retentive of its O isotopic composition. A crystal with a 100 μ radius will retain its initial core composition for just over 10 million years at 600 °C.},
            author = {D. K. Moore and Daniele J. Cherniak and E. B. Watson},
            doi = {10.2138/am-1998-7-803},
            issn = {0003004X},
            journal = {American Mineralogist},
            pages = {700-711},
            title = {Oxygen diffusion in rutile from 750 to 1000 °C and 0.1 to 1000 MPa},
            volume = {83},
            year = {1998},
            }
          ",
    )

    return data, info
end

"""
    Rt_O_Venkatu1970

Diffusion data of O in natural rutile.  
Calibrated between 900-1300C. From Venkatu et al. (1970).
"""
function Rt_O_Venkatu1970()
    data = DiffusionData(
        Name = "O self-diffusion in Rutile | Venkatu et al. (1970)",
        Phase = "Rutile",
        Formula = "TiO2",
        Species = "O",
        Crystallography = "Tetragonal",
        D0 = 6.4e-2cm^2 / s,
        Ea = 61400u"cal/mol",
        Ea_1σ = 22.0e3J / mol,
        Charge = -2,
        T_range_min = 900C,
        T_range_max = 1300C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by ML (05.02.25)",
        BibTex_Reference = "
            @article{Venkatu1970,
            author = {D.A. Venkatu and L.E. Poteat},
            journal = {Material Science and Engineering},
            pages = {258-262},
            title = {Diffusion of Titanium in Single Crystal Rutile},
            volume = {5},
            year = {1970},
            }
          ",
    )

    return data, info
end
