"""
    Ol_Mg_Chakraborty1994_forsterite()

Trace diffusion data of Mg in synthetic pure forsterite. Calibrated with experiments conducted between 1000-1300°C at atmospheric pressure and up to 10 GPa in multianvil apparatus.
From Chakraborty and al. (1994) (https://doi.org/10.1007/BF00203923).
"""
function Ol_Mg_Chakraborty1994_forsterite()
    data = DiffusionData(
        Name = "Mg diffusion in synthetic pure forsteriste | Chakraborty and al. (1994)",
        Phase = "Olivine",
        Formula = "(Mg,Fe)2[SiO4]",
        Species = "Mg",
        Orientation = "⊥c",
        Crystallography = "Orthorhombic",
        Buffer = "Air and CO-CO2 mixture",
        D0 = 9.6*1e-4m^2 / s,
        Ea = 400u"kJ / mol",
        Ea_1σ = (60/2)u"kJ / mol",
        Charge = 2,  # charge of the cation
        T_range_min = 1000C,
        T_range_max = 1300C
    )
    info = MaterialParamsInfo(;
        Comment = "Checked values by HD (19.03.25)",
        BibTex_Reference = "
            @article{chakraborty1994mg,
            title={Mg tracer diffusion in synthetic forsterite and San Carlos olivine as a function of P, T and fO 2},
            author={Chakraborty, Sumit and Farver, John R and Yund, Richard A and Rubie, David C},
            journal={Physics and Chemistry of Minerals},
            volume={21},
            pages={489--500},
            year={1994},
            publisher={Springer}
            }
          ",
    )

    return data, info
end

