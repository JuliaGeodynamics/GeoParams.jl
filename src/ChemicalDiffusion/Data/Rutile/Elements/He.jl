"""
    Rt_He_Cherniak2011_para_c

Diffusion data of He in rutile. With anhydrous conditions and parallel to c-axis. 
Calibrated between 250-500°C. From Cherniak et al. (2011) 
(https://doi.org/10.1016/j.chemgeo.2011.07.015).
"""
function Rt_He_Cherniak2011_para_c()
    data = DiffusionData(
        Name = "He diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Cherniak et al. (2011)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "He",  # element or species being diffused
        Orientation = "Ξc",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 1.75e-8m^2 / s,  # pre-exponential factor
        log_D0_1σ = 1.3NoUnits,
        Ea = 120kJ / mol,  # activation energy
        Ea_1σ = 7kJ / mol,  # uncertainty at 1σ of the activation energy
        Charge = 4,  # charge of the cation
        T_range_min = 250C,  # temperature min of the experiment
        T_range_max = 500C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (15.01.25)",
        BibTex_Reference = "
          @article{Cherniak2011,
            abstract = {Diffusion of helium has been characterized in natural sphene (titanite) and synthetic rutile. Polished slabs of sphene and rutile were implanted with 100keV 3He at a dose of 5×1015 3He/cm2 and annealed in 1-atm furnaces. 3He distributions following experiments were measured with nuclear reaction analysis using the reaction 3He(d,p)4He. For diffusion in rutile we obtain the following Arrhenius relations:D||a=2.26 × 10-10exp(-126± 11 kJmol-1/RT) m2sec-1.D||c=1.75 × 10-8exp(-120± 7 kJmol-1/RT) m2sec-1.Although activation energies for diffusion parallel to the c- and a-axes are comparable, there is marked diffusional anisotropy, with diffusion parallel to the c-axis about 2 orders of magnitude faster than transport parallel to the a-axis. These diffusivities bracket the values determined for He diffusion in rutile in bulk release experiments (Stockli et al., 2005, 2007; Wolfe, 2009), although the role of anisotropy could not be directly evaluated in those measurements.In titanite, the following Arrhenius relation was obtained over the temperature range 252-550 °C for diffusion parallel to the a-axis:. D=2.14 × 10-6exp(-148± 8 kJmol-1/RT) m2sec-1.In contrast to rutile and zircon (Cherniak et al., 2009), titanite shows little evidence of anisotropy, as diffusivities parallel to the a- and c-axes are similar, and diffusivities for titanites from two different localities are similar. He diffusion coefficients obtained in this study are similar to those measured through bulk release of He by step heating (Shuster et al., 2003). Over the investigated temperature range, diffusion of He in titanite is similar to that of He diffusion in rutile parallel to the c-axis, but much faster than diffusion parallel to the a-axis.Since the diffusion of He in rutile exhibits such pronounced anisotropy, we model diffusional loss of He with a recently developed finite-element code (Watson et al., 2010) created to simulate diffusion in cylindrical geometry with differing radial and axial diffusion coefficients. We present example applications to evaluate helium losses from rutile grains as a function of grain size and length to diameter ratios.In efforts to better understand the occurrence of pronounced anisotropy for He diffusion in some crystals (e.g., rutile, zircon) but not in others (e.g., apatite, titanite), we consider the density and distribution of interstitial apertures in the crystal structure that might permit He migration. These determinations, an extension of the concept of ionic porosity, are consistent with observations of relative diffusivities and the existence or absence of significant anisotropy for the minerals rutile, zircon, apatite and titanite. © 2011 Elsevier B.V.},
            author = {D. J. Cherniak and E. B. Watson},
            doi = {10.1016/j.chemgeo.2011.07.015},
            issn = {00092541},
            issue = {3-4},
            journal = {Chemical Geology},
            keywords = {Diffusion,Helium,Nuclear reaction analysis,Rutile,Thermochronology,Titanite (sphene)},
            pages = {149-161},
            publisher = {Elsevier B.V.},
            title = {Helium diffusion in rutile and titanite, and consideration of the origin and implications of diffusional anisotropy},
            volume = {288},
            url = {http://dx.doi.org/10.1016/j.chemgeo.2011.07.015},
            year = {2011},
            }
          ",
    )

    return data, info
end

"""
    Rt_He_Cherniak2011_perp_c

Diffusion data of He in rutile. With anhydrous conditions and perpendicular to c-axis. 
Calibrated between 300-600°C. From Cherniak et al. (2011) 
(https://doi.org/10.1016/j.chemgeo.2011.07.015).
"""
function Rt_He_Cherniak2011_perp_c()
    data = DiffusionData(
        Name = "He diffusion in Rutile (perpendicular to c-axis and anhydrous conditions) | Cherniak et al. (2011)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "He",  # element or species being diffused
        Orientation = "⊥c",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 2.48e-10m^2 / s,  # pre-exponential factor
        log_D0_1σ = 1.75NoUnits,
        Ea = 126kJ / mol,  # activation energy
        Ea_1σ = 11kJ / mol,  # uncertainty at 1σ of the activation energy
        Charge = 4,  # charge of the cation
        T_range_min = 300C,  # temperature min of the experiment
        T_range_max = 600C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (15.01.25)",
        BibTex_Reference = "
          @article{Cherniak2011,
            abstract = {Diffusion of helium has been characterized in natural sphene (titanite) and synthetic rutile. Polished slabs of sphene and rutile were implanted with 100keV 3He at a dose of 5×1015 3He/cm2 and annealed in 1-atm furnaces. 3He distributions following experiments were measured with nuclear reaction analysis using the reaction 3He(d,p)4He. For diffusion in rutile we obtain the following Arrhenius relations:D||a=2.26 × 10-10exp(-126± 11 kJmol-1/RT) m2sec-1.D||c=1.75 × 10-8exp(-120± 7 kJmol-1/RT) m2sec-1.Although activation energies for diffusion parallel to the c- and a-axes are comparable, there is marked diffusional anisotropy, with diffusion parallel to the c-axis about 2 orders of magnitude faster than transport parallel to the a-axis. These diffusivities bracket the values determined for He diffusion in rutile in bulk release experiments (Stockli et al., 2005, 2007; Wolfe, 2009), although the role of anisotropy could not be directly evaluated in those measurements.In titanite, the following Arrhenius relation was obtained over the temperature range 252-550 °C for diffusion parallel to the a-axis:. D=2.14 × 10-6exp(-148± 8 kJmol-1/RT) m2sec-1.In contrast to rutile and zircon (Cherniak et al., 2009), titanite shows little evidence of anisotropy, as diffusivities parallel to the a- and c-axes are similar, and diffusivities for titanites from two different localities are similar. He diffusion coefficients obtained in this study are similar to those measured through bulk release of He by step heating (Shuster et al., 2003). Over the investigated temperature range, diffusion of He in titanite is similar to that of He diffusion in rutile parallel to the c-axis, but much faster than diffusion parallel to the a-axis.Since the diffusion of He in rutile exhibits such pronounced anisotropy, we model diffusional loss of He with a recently developed finite-element code (Watson et al., 2010) created to simulate diffusion in cylindrical geometry with differing radial and axial diffusion coefficients. We present example applications to evaluate helium losses from rutile grains as a function of grain size and length to diameter ratios.In efforts to better understand the occurrence of pronounced anisotropy for He diffusion in some crystals (e.g., rutile, zircon) but not in others (e.g., apatite, titanite), we consider the density and distribution of interstitial apertures in the crystal structure that might permit He migration. These determinations, an extension of the concept of ionic porosity, are consistent with observations of relative diffusivities and the existence or absence of significant anisotropy for the minerals rutile, zircon, apatite and titanite. © 2011 Elsevier B.V.},
            author = {D. J. Cherniak and E. B. Watson},
            doi = {10.1016/j.chemgeo.2011.07.015},
            issn = {00092541},
            issue = {3-4},
            journal = {Chemical Geology},
            keywords = {Diffusion,Helium,Nuclear reaction analysis,Rutile,Thermochronology,Titanite (sphene)},
            pages = {149-161},
            publisher = {Elsevier B.V.},
            title = {Helium diffusion in rutile and titanite, and consideration of the origin and implications of diffusional anisotropy},
            volume = {288},
            url = {http://dx.doi.org/10.1016/j.chemgeo.2011.07.015},
            year = {2011},
            }
          ",
    )

    return data, info
end
