"""
    Rt_Sc_Sasaki1985_para_c

Diffusion data of Sc in rutile. With anhydrous conditions, in air and parallel to c-axis.
Calibrated between 900-1500°C. From Sasaki et al. (1985) 
(https://doi.org/10.1016/0022-3697(85)90129-5).
"""
function Rt_Sc_Sasaki1985_para_c()
    data = DiffusionData(
        Name = "Sc diffusion in Rutile (parallel to c-axis and anhydrous conditions) | Sasaki et al. (1985)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Sc",  # element or species being diffused
        Orientation = "Ξc",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "air",
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 1.98e-8m^2 / s,  # pre-exponential factor
        Ea = 188kJ / mol,  # activation energy
        Charge = 3,  # charge of the cation
        T_range_min = 900C,  # temperature min of the experiment
        T_range_max = 1500C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Values re-fitted by ML (10.02.25)",
        BibTex_Reference = "
            @article{Sasaki1985,
            abstract = {By means of the radioactive-tracer sectioning technique, the tracer diffusion of the impurity ions, 46Sc, 51Cr, 54Mn, 59Fe, 60Co, 63Ni and 95Zr, in rutile single crystals was measured as functions of crystal orientation, temperature, oxygen partial pressure and Al impurity content. The diffusion coefficients are very sensitive to the electric charge of the impurity ions. Divalent impurities (e.g., Co and Ni) diffuse extremely rapidly in TiO2, compared to cation self-diffusion, and exhibit an extreme anisotropy in diffusion behavior, divalent-impurity diffusion parallel to the c-axis is much larger than it is perpendicular to the c-axis. Trivalent impurity ions (Sc and Cr) and tetravalent impurity ions (Zr) diffuse similar to cation self-diffusion, both as functions of temperature and oxygen partial pressure. The divalent impurity ions Co and Ni apparently diffuse as interstitial ions along open channels parallel to the c-axis. The results suggest that Sc, Cr and Zr ions diffuse by an interstitialcy mechanism involving the simultaneous and cooperative migration of tetravalent interstitial titanium ions and the tracer-impurity ions. Iron ions diffuse both as divalent and as trivalent ions. The impurity diffusion as functions of oxygen partial pressure and Al-impurity content are consistent with calculations of point-defect concentrations in rutile. © 1985.},
            author = {Jun Sasaki and N. L. Peterson and K. Hoshino},
            doi = {10.1016/0022-3697(85)90129-5},
            issn = {00223697},
            issue = {11},
            journal = {Journal of Physics and Chemistry of Solids},
            keywords = {diffusion mechanisms,interstitialcy diffusion mechanism,point defects,rutile (TiO2),tracer impurity diffusion},
            pages = {1267-1283},
            title = {Tracer impurity diffusion in single-crystal rutile (TiO2-x)},
            volume = {46},
            year = {1985},
            }
          ",
    )

    return data, info
end

"""
    Rt_Sc_Sasaki1985_perp_c

Diffusion data of Sc in rutile. With anhydrous conditions, in air and perpendicular to c-axis.
Calibrated between 900-1100°C. From Sasaki et al. (1985) 
(https://doi.org/10.1016/0022-3697(85)90129-5).
"""
function Rt_Sc_Sasaki1985_perp_c()
    data = DiffusionData(
        Name = "Sc diffusion in Rutile (perpendicular to c-axis and anhydrous conditions) | Sasaki et al. (1985)",
        Phase = "Rutile",  # name of the mineral
        Formula = "TiO2",  # chemical formula of the mineral
        Species = "Sc",  # element or species being diffused
        Orientation = "⊥c",  # Crystal orientation from the diffusion experiment
        Crystallography = "Tetragonal",  # Crystallographic system of the mineral
        Buffer = "air",
        Fluid = "Anhydrous",  # Fluid condition (e.g., anhydrous) during the experiment
        D0 = 1.57e-9m^2 / s,  # pre-exponential factor
        Ea = 169kJ / mol,  # activation energy
        Charge = 3,  # charge of the cation
        T_range_min = 900C,  # temperature min of the experiment
        T_range_max = 1100C  # temperature max of the experiment
    )
    info = MaterialParamsInfo(;
        Comment = "Values re-fitted by ML (10.02.25)",
        BibTex_Reference = "
            @article{Sasaki1985,
            abstract = {By means of the radioactive-tracer sectioning technique, the tracer diffusion of the impurity ions, 46Sc, 51Cr, 54Mn, 59Fe, 60Co, 63Ni and 95Zr, in rutile single crystals was measured as functions of crystal orientation, temperature, oxygen partial pressure and Al impurity content. The diffusion coefficients are very sensitive to the electric charge of the impurity ions. Divalent impurities (e.g., Co and Ni) diffuse extremely rapidly in TiO2, compared to cation self-diffusion, and exhibit an extreme anisotropy in diffusion behavior, divalent-impurity diffusion parallel to the c-axis is much larger than it is perpendicular to the c-axis. Trivalent impurity ions (Sc and Cr) and tetravalent impurity ions (Zr) diffuse similar to cation self-diffusion, both as functions of temperature and oxygen partial pressure. The divalent impurity ions Co and Ni apparently diffuse as interstitial ions along open channels parallel to the c-axis. The results suggest that Sc, Cr and Zr ions diffuse by an interstitialcy mechanism involving the simultaneous and cooperative migration of tetravalent interstitial titanium ions and the tracer-impurity ions. Iron ions diffuse both as divalent and as trivalent ions. The impurity diffusion as functions of oxygen partial pressure and Al-impurity content are consistent with calculations of point-defect concentrations in rutile. © 1985.},
            author = {Jun Sasaki and N. L. Peterson and K. Hoshino},
            doi = {10.1016/0022-3697(85)90129-5},
            issn = {00223697},
            issue = {11},
            journal = {Journal of Physics and Chemistry of Solids},
            keywords = {diffusion mechanisms,interstitialcy diffusion mechanism,point defects,rutile (TiO2),tracer impurity diffusion},
            pages = {1267-1283},
            title = {Tracer impurity diffusion in single-crystal rutile (TiO2-x)},
            volume = {46},
            year = {1985},
            }
          ",
    )

    return data, info
end