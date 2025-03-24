"""
    Ol_Fe_Mg_Dohmen2007_perp_c()

Interdiffusion data of Fe and Mg in natural gem quality olivine for all oxygen fugacities. Calibrated with experiments conducted between 700-1200°C at atmospheric pressure with oxygen fugacity between 1e-12 to 1e-5 and parallel to c axis from Dohmen et al. (2007) and data from the literature. Note that the molar fraction of Fe is required to calculate the interdiffusion coefficient. To calculate the interdiffusion coefficient parallel to c, multiply the interdiffusion coefficient by 1e6.
From Dohmen and Chakraborty (2007) (https://doi.org/10.1007/s00269-007-0158-6).
"""
function Ol_Fe_Mg_Dohmen2007_perp_c()
    data = DiffusionData(
        Name = "General Fe-Mg interdiffusion in olivine | Dohmen and Chakraborty (2007)",
        Phase = "Olivine",
        Formula = "(Mg,Fe)2[SiO4]",
        Species = "Mg and Fe",
        Orientation = "⊥c",
        Crystallography = "Orthorhombic",
        Buffer = "CO-CO2 mixture",
        D0 = (5.37 * 1.0e-9 * 2.303)m^2 / s,
        Ea = (226000)u"J / mol",
        ΔV = (7 * 1.0e-6)m^3 / mol,
        aX = (3 * 2.303)NoUnits,
        bX = -0.14NoUnits,
        P0 = 1bar,
        Charge = 2,  # charge of the cation
        T_range_min = 700C,
        T_range_max = 1550C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (21.03.25)",
        BibTex_Reference = "
            @article{dohmen2007fe,
            title={Fe--Mg diffusion in olivine II: point defect chemistry, change of diffusion mechanisms and a model for calculation of diffusion coefficients in natural olivine},
            author={Dohmen, Ralf and Chakraborty, Sumit},
            journal={Physics and Chemistry of Minerals},
            volume={34},
            number={6},
            pages={409--430},
            year={2007},
            publisher={Springer}
            }
          ",
    )

    return data, info
end

"""
    Ol_Fe_Mg_Dohmen2007_TaMED_perp_c()

Interdiffusion data of Fe and Mg for the transition metal extrinsic (TaMED) mechanism. Calibrated with experiments conducted between 900-1200°C at atmospheric pressure with oxygen fugacity between 1e-10 to 1e-5 and parallel to c axis from Dohmen et al., 2007 and data from the literature. Note that the molar fraction of Fe and the oxygen fugacity are required to calculate the interdiffusion coefficient. To calculate the interdiffusion coefficient parallel to c, multiply the interdiffusion coefficient by 1e6.
From Dohmen and Chakraborty (2007) (https://doi.org/10.1007/s00269-007-0158-6).
"""
function Ol_Fe_Mg_Dohmen2007_TaMED_perp_c()
    data = DiffusionData(
        Name = "Fe-Mg interdiffusion in olivine of the TaMED mechanism | Dohmen and Chakraborty (2007)",
        Phase = "Olivine",
        Formula = "(Mg,Fe)2[SiO4]",
        Species = "Mg and Fe",
        Orientation = "⊥c",
        Crystallography = "Orthorhombic",
        Buffer = "CO-CO2 mixture",
        D0 = (6.17 * 1.0e-10 * 2.303)m^2 / s,
        Ea = (201000)u"J / mol",
        ΔV = (7 * 1.0e-6)m^3 / mol,
        nfO2 = (1 / 6)NoUnits,
        dfO2 = 1.0e-7NoUnits,
        aX = (3 * 2.303)NoUnits,
        bX = -0.1NoUnits,
        P0 = 1bar,
        Charge = 2,  # charge of the cation
        T_range_min = 900C,
        T_range_max = 1550C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (21.03.25)",
        BibTex_Reference = "
            @article{dohmen2007fe,
            title={Fe--Mg diffusion in olivine II: point defect chemistry, change of diffusion mechanisms and a model for calculation of diffusion coefficients in natural olivine},
            author={Dohmen, Ralf and Chakraborty, Sumit},
            journal={Physics and Chemistry of Minerals},
            volume={34},
            number={6},
            pages={409--430},
            year={2007},
            publisher={Springer}
            }
          ",
    )

    return data, info
end


"""
    Ol_Fe_Mg_Dohmen2007_PED_perp_c()

Interdiffusion data of Fe and Mg in natural gem quality olivine for the purely extrinsic (PED) mechanism. Calibrated with experiments conducted between 700-900°C at atmospheric pressure with oxygen fugacity between 1e-12 to 1e-10 and parallel to c axis from Dohmen et al. (2007) and data from the literature. Note that the molar fraction of Fe is required to calculate the interdiffusion coefficient. To calculate the interdiffusion coefficient parallel to c, multiply the interdiffusion coefficient by 1e6.
From Dohmen and Chakraborty (2007) (https://doi.org/10.1007/s00269-007-0158-6).
"""
function Ol_Fe_Mg_Dohmen2007_PED_perp_c()
    data = DiffusionData(
        Name = "Fe-Mg interdiffusion in olivine of the PED mechanism | Dohmen and Chakraborty (2007)",
        Phase = "Olivine",
        Formula = "(Mg,Fe)2[SiO4]",
        Species = "Mg and Fe",
        Orientation = "⊥c",
        Crystallography = "Orthorhombic",
        Buffer = "CO-CO2 mixture",
        D0 = (1.23 * 1.0e-9 * 2.303)m^2 / s,
        Ea = (220000)u"J / mol",
        ΔV = (7 * 1.0e-6)m^3 / mol,
        aX = (3 * 2.303)NoUnits,
        bX = -0.1NoUnits,
        P0 = 1bar,
        Charge = 2,  # charge of the cation
        T_range_min = 700C,
        T_range_max = 900C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (21.03.25)",
        BibTex_Reference = "
            @article{dohmen2007fe,
            title={Fe--Mg diffusion in olivine II: point defect chemistry, change of diffusion mechanisms and a model for calculation of diffusion coefficients in natural olivine},
            author={Dohmen, Ralf and Chakraborty, Sumit},
            journal={Physics and Chemistry of Minerals},
            volume={34},
            number={6},
            pages={409--430},
            year={2007},
            publisher={Springer}
            }
          ",
    )

    return data, info
end
